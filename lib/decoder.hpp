#pragma once

#include "libONIAK/oniakDataStructure/oupq.h"
#include "libONIAK/oniakHash/ohash.h"
#include "libONIAK/oniakMath/orange.h"

#include <doro/doro.hpp>
#include <doro/iblt.h>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <format>
#include <functional>
#include <iostream>
#include <map>
#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Doro {
    enum class PursuitChoice {
        L1,
        L2,
    };
    struct DecodeConfig {
        bool verbose, debug;
        // lower and upper bound of value range
        double lb, ub;
        int max_num_peels = -1;
        PursuitChoice pursuit_choice;
    };

    template <typename IndexType>
    IndexType from_binary_vector(const std::vector<uint8_t>& vec) {
        assert(vec.size() == sizeof(IndexType));
        IndexType result;
        std::memcpy(&result, vec.data(), vec.size());
        return result;
    }


    template <typename DoroCodeT, ONIAK::UpdatePQBackend backend = ONIAK::UpdatePQBackend::PriorityQueue>
    class DoroDecoder {
    public:
        using ArrType = typename DoroCodeT::ArrayT;
        using IndexType = typename DoroCodeT::IndexT;
        using TwoDimVector = std::vector<std::vector<IndexType>>;
        using UpdatePQ = ONIAK::UpdatePQAdapter<IndexType, ArrType, backend>::type;
        using IBLTResultType = std::unordered_map<IndexType, std::vector<uint8_t>>;

        DoroDecoder(DoroCodeT &code, const std::unordered_set<IndexType> &setA,
                    DecodeConfig &config, int max_recenter_rounds = 10,
                    ONIAK::WYHash *finger_hash = nullptr,
                    ONIAK::WYHash *resolving_hash = nullptr)
            : code_(&code),
              setA_(&setA), k_is_five_(code.k() == 5), valid_neighbors_(false),
              max_recenter_rounds_(max_recenter_rounds),
              finger_hash_(finger_hash), resolving_hash_(resolving_hash),
              neighbors_(code.size()), neighbors2_(code.size()), result_(),
              fingerprints_(), my_fingerprints_(), new_elements_(),
              pairs_5in3_(), affected_neighbors_(),
              priority_queue_(
                  // function from sensed strength to priority
                  [](ArrType a) { return std::abs(a); }),
              config_(&config), collision_resolving_(false),
              unresolved_elements_(), unresolved_ids_() {}

        void peel_until_empty() {
            while (!priority_queue_.empty()) {
                auto [strength, element] = priority_queue_.top();
                new_elements_.insert(element);

                ArrType delta = strength_to_delta(element, strength);
                code_->peel(element, delta);
                auto cur_iter = result_.find(element);
                ArrType cur_element_value = (cur_iter != result_.end()) ? cur_iter->second : 0;
                result_[element] = cur_element_value + delta;

                affected_neighbors_.clear();
                auto all_hashes = code_->hash_all(element);
                for (auto [index, sign] : all_hashes) {
                    notify_neighbors(index, -sign * delta);  // minus because the counter is reduced by delta
                }
                update_neighbor_strengths(element, /*append53*/ true);

                if (config_->verbose) {
                    std::cout << std::format("power: {}, element: {}, cur_element_value: {}, delta: {}\n",
                        strength, element, cur_element_value, delta);
                    code_->show_result();
                }
            }
        }

        int decode() {
            priority_queue_.clear();
            bool is_l2 = config_->pursuit_choice == PursuitChoice::L2; // else is l1

            scan_setA(is_l2);
            // peel the counters until nothing can be done.
            new_elements_.clear();
            for (int rnd = 0; rnd < max_recenter_rounds_; ++rnd) {
                peel_until_empty();

                ArrType delta = 0;
                if (collision_resolving_) {
                    delta = recenter_code();
                    if (config_->pursuit_choice == PursuitChoice::L2) {
                        config_->pursuit_choice = PursuitChoice::L1;
                        scan_setA(false);
                    }
                }
                if (config_->max_num_peels > 0 && code_->num_peels() >= config_->max_num_peels) break;
                if (delta == 0 && priority_queue_.empty()) break;
            }
            if (collision_resolving_ && !pairs_5in3_.empty()) {
                resolve_5in3();
            }

            unresolved_elements_.clear();
            unresolved_ids_.clear();
            if (collision_resolving_ && finger_hash_ != nullptr && resolving_hash_ != nullptr) {
                add_unresolved_elements();
            }
            return code_->num_peels();
        }

        // fingerprints are add-only.
        // In other words, if an element is added and then removed, its fingerprint will not be removed.
        // This is because it costs more to transmit deletions.
        std::vector<int8_t> update_fingerprints() {
            assert(finger_hash_ != nullptr && resolving_hash_ != nullptr);
            // i.e., finger_s in two_party
            size_t capacity = finger_hash_->mod();
            std::vector<int8_t> fingerprints(capacity, 0);
            for (auto [key, value] : result_) {
                uint32_t finger = finger_hash_->hash_in_range(key);
                if (std::abs(value) < 1e-6)  continue;
                if (!my_fingerprints_.contains(finger)) {
                    fingerprints[finger] = 1;
                }
                my_fingerprints_.insert({ finger, key });
            }
            return fingerprints;
        }

        void load_fingerprints(const std::vector<int8_t> fingerprints) {
            assert(finger_hash_ != nullptr);
            assert(fingerprints.size() == finger_hash_->mod());
            for (size_t i : ONIAK::nonzero_indices(fingerprints)) {
                if (!fingerprints_.contains(i)) {
                    fingerprints_.insert(i);
                }
            }
        }

        std::unordered_map<IndexType, ArrType>& result() {
            return result_;
        }

        // returns number of detected collisions
        int resolve_collision(DoroDecoder& other) {
            assert(finger_hash_ != nullptr && *other.finger_hash_ == *finger_hash_);
            assert(resolving_hash_ != nullptr && *other.resolving_hash_ == *resolving_hash_);
            int num_collisions = 0;

            for (auto [i, finger_i] : std::views::enumerate(other.unresolved_elements_)) {
                auto [finger1, finger2] = finger_i;
                auto [finger_iter, iter_end] = my_fingerprints_.equal_range(finger1);
                for (; finger_iter != iter_end; ++finger_iter) {
                    IndexType element = finger_iter->second;
                    IndexType other_element = other.unresolved_ids_[i];
                    uint64_t f2_element = resolving_hash_->hash64(element);
                    if (f2_element == finger2) {  // collision detected, reverting
                        // CAUTION: if there is any collision with finger2, then set reconciliation would fail.
                        assert(element == other_element);
                        ++num_collisions;
                        ArrType cur_value = result_.at(element);
                        ArrType other_value = other.result_.at(other_element);
                        result_[element] = 0;
                        other.result_[element] = 0;
                        code_->peel(element, -cur_value);
                        other.code_->peel(other_element, -other_value);
                        break;
                    } // else accept the change, nothing is needed
                }
            }
            other.unresolved_elements_.clear();
            return num_collisions;
        }

        void resolve_5in3() {
            if constexpr (!std::is_integral_v<ArrType>) return;
            for (auto [e1, e2] : pairs_5in3_) {
                ArrType e1_val = result_.contains(e1) ? result_.at(e1) : 0;
                ArrType e2_val = result_.contains(e2) ? result_.at(e2) : 0;
                if (!e1_val && !e2_val) continue;
                assert(!e1_val || !e2_val);
                ArrType val = e1_val ? e1_val : e2_val;
                IndexType cur = e1_val ? e1 : e2;
                IndexType other = e1_val ? e2 : e1;
                ArrType strength = (code_->sense(cur) - code_->sense(other)) * (-val);
                if (strength >= 3) {
                    code_->peel(cur, -val);
                    code_->peel(other, val);
                    result_[cur] = 0;
                    result_[other] = val;
                    new_elements_.insert(other);  // what if other is already decoded by the other party?

                    affected_neighbors_.clear();
                    auto all_hashes = code_->hash_all(cur);
                    for (auto [index, sign] : all_hashes) {
                        notify_neighbors(index, sign * val);  // each counter value increases by val.
                    }
                    all_hashes = code_->hash_all(other);
                    for (auto [index, sign] : all_hashes) {
                        notify_neighbors(index, -sign * val);
                    }
                    update_neighbor_strengths(-1);
                }
            }
        }

        void encode_iblt(IBLT<IndexType>& iblt) const {
            for (const IndexType& element : *setA_) {
                // elements that are assumed not to be in result
                if (!(result_.contains(element) && result_.at(element) != 0)) {
                    iblt.insert(static_cast<uint64_t>(element), to_binary_vector(element));
                }
            }
        }

        std::vector<IndexType> decode_iblt(IBLT<IndexType> iblt) {
            // compute difference iblt
            for (const IndexType& element : *setA_) {
                if (!(result_.contains(element) && result_.at(element) != 0)) {
                    iblt.erase(element, to_binary_vector(element));
                }
            }
            IBLTResultType positive, negative;
            // return if failure
            if (!iblt.listEntries(positive, negative)) return {};

            std::vector<IndexType> other_result;
            ArrType my_value = (config_->ub > 0) ? config_->ub : config_->lb;
            // negative entries are my extra ones
            for (const auto& [key, value] : negative) {
                IndexType element = from_binary_vector<IndexType>(value);
                // special case, if we have no error, it must be from the other side.
                if (my_value == 0) {
                    other_result.push_back(element);
                } else {
                    // here, doro is not updated, since we no longer attempt to peel them.
                    result_[element] = my_value;
                    new_elements_.insert(element);
                }
            }  // these new elements still may collide with those in the other party.
            add_unresolved_elements();

            // positive entries are the other's extra ones

            for (const auto& [key, value] : positive) {
                IndexType element = from_binary_vector<IndexType>(value);
                if (setA_->contains(element)) {  // false positive is at my side
                    result_[element] = 0;
                } else
                    other_result.push_back(element);
            }
            return other_result;
        }

        void add_iblt_elements(std::vector<IndexType> list) {
            ArrType my_value = (config_->ub > 0) ? config_->ub : config_->lb;
            for (const IndexType& element : list) {
                // considering the special case above, we flip our result
                if (result_.contains(element) && result_.at(element) != 0) result_[element] = 0;
                else
                    result_[element] = my_value;
            }
        }

        void enter_resolving() {
            collision_resolving_ = true;
        }

        std::vector<std::pair<int, uint64_t>>& unresolved_elements() {
            return unresolved_elements_;
        }

        bool has_new_elements() const {
            return !new_elements_.empty();
        }

        void add_to_priority_queue(IndexType element, ArrType strength, ArrType delta) {
            // if fingerprint mechanism is active, then this element must not collide with any known fingerprint
            // this is set to avoid the case in which one party adds an element, and another party subtracts it.
            bool fingerprint_condition = (finger_hash_ == nullptr) || !fingerprints_.contains(finger_hash_->hash_in_range(element));
            bool delta_condition = (std::abs(delta) > 1e-6);
            if ((collision_resolving_ || fingerprint_condition) && delta_condition)
                priority_queue_.push(element, strength);
            else
                priority_queue_.set(element, strength);
        }

        std::unordered_set<int>& fingerprints() {
            return fingerprints_;
        }

        const std::unordered_set<int>& fingerprints() const {
            return fingerprints_;
        }

    private:
        ArrType get_delta(double delta, ArrType cur_value) const {
            double new_value = cur_value + delta;
            new_value = std::max(config_->lb, new_value);
            new_value = std::min(config_->ub, new_value);
            delta = new_value - cur_value;
            if constexpr (std::is_integral_v<ArrType>) {
                if (std::abs(delta) < 0.50001)
                    delta = 0; // if delta is 0.5, default round to 1
                else         // but we let it be 0 to avoid oscillation
                    delta = std::round(delta);
            }
            return delta;
        }

        ArrType strength_to_delta(IndexType element, ArrType strength) const {
            bool is_l2 = config_->pursuit_choice == PursuitChoice::L2;
            auto cur_iter = result_.find(element);
            ArrType cur_element_value = (cur_iter != result_.end()) ? cur_iter->second : 0;
            // In L2, delta is computed from average instead of sum
            double normalized_strength = is_l2 ? static_cast<double>(strength) / static_cast<double>(code_->k()) : code_->get_median(element);
            return get_delta(normalized_strength, cur_element_value);
        }

        ArrType new_strength2(IndexType element, ArrType strength, ArrType delta) const {
            bool is_l2 = config_->pursuit_choice == PursuitChoice::L2;
            ArrType new_str = is_l2 ? strength + delta : code_->sense_l1(element);
            assert(!is_l2 || new_str == code_->sense(element));
            // the bookkeeping of strengths should always be correct. Assertion is only checked in debug build.
            return new_str;
        }

        void notify_neighbors(int counter_idx, ArrType strength_change) {
            for (IndexType neighbor : neighbors_[counter_idx]) {
                if (!affected_neighbors_.contains(neighbor)) affected_neighbors_[neighbor] = 0;
                affected_neighbors_[neighbor] += strength_change;
            }
            for (IndexType neighbor : neighbors2_[counter_idx]) {
                if (!affected_neighbors_.contains(neighbor)) affected_neighbors_[neighbor] = 0;
                affected_neighbors_[neighbor] -= strength_change;
            }
        }

        void update_neighbor_strengths(IndexType element, bool append_53 = false) {
            for (auto [neighbor, strength_change] : affected_neighbors_) {
                if (append_53 && k_is_five_ && element != neighbor && std::abs(strength_change) >= 3) {
                    std::pair<IndexType, IndexType> ele_nei = (element < neighbor) ? std::make_pair(element, neighbor) : std::make_pair(neighbor, element);
                    if (!pairs_5in3_.contains(ele_nei)) {
                        pairs_5in3_.insert(ele_nei);
                    }
                }
                ArrType neighbor_strength = new_strength2(neighbor, priority_queue_[neighbor], strength_change);
                ArrType neighbor_delta = strength_to_delta(neighbor, neighbor_strength);
                add_to_priority_queue(neighbor, neighbor_strength, neighbor_delta);
            }
            affected_neighbors_.clear();
        }

        // in case of quantization error, recenter all counters to within the range [lb, ub)
        ArrType recenter_code() {
            affected_neighbors_.clear();
            ArrType max_delta = std::ranges::max(code_->code(), {},
                std::bind(&DoroCodeT::deviation, code_, std::placeholders::_1)); // calls member function
            max_delta = code_->deviation(max_delta);
            if (max_delta <= 1) return 0;
            for (auto [i, cur_value] : code_->code() | std::views::enumerate) {
                if (code_->deviation(cur_value) > 0)
                    if (code_->deviation(cur_value) == max_delta) {
                        ArrType new_value = code_->recenter(cur_value);
                        ArrType delta = new_value - cur_value;
                        code_->code()[i] = new_value;
                        ++code_->num_recenters();
                        notify_neighbors(i, delta);
                    }
            }
            update_neighbor_strengths(-1);
            return max_delta;
        }

        void scan_setA(bool is_l2) {
            for (IndexType element : *setA_) {
                if (!valid_neighbors_) {
                    auto all_hashes = code_->hash_all(element);
                    for (auto [index, sign] : all_hashes) {
                        if (sign > 0)
                            neighbors_[index].push_back(element);
                        else neighbors2_[index].push_back(element);
                    }
                }
                ArrType strength = is_l2 ? code_->sense(element) : code_->sense_l1(element);
                ArrType delta = strength_to_delta(element, strength);
                add_to_priority_queue(element, strength, delta);
            }
            valid_neighbors_ = true;
        }

        void add_unresolved_elements() {
            for (auto key : new_elements_) {
                if (std::abs(result_.at(key)) > 1e-6) {
                    int finger1 = finger_hash_->hash_in_range(key);
                    if (fingerprints_.contains(finger1)) {
                        uint64_t finger2 = resolving_hash_->hash64(key);
                        unresolved_elements_.push_back({ finger1, finger2 });
                        unresolved_ids_.push_back(key);
                    }
                }
            }
        }

        DoroCodeT* code_;
        const std::unordered_set<IndexType>* setA_;
        bool k_is_five_;  // a special case for handling rarely colliding element pairs when k = 5
        bool valid_neighbors_; // the neighbors set are already established
        int max_recenter_rounds_;
        ONIAK::WYHash* finger_hash_, * resolving_hash_;
        // neighbors_[index] stores elements that have a plus sign at index
        // neighbors2_[index] stores those with a minus sign
        TwoDimVector neighbors_, neighbors2_;
        // decoded result
        std::unordered_map<IndexType, ArrType> result_;
        // fingerprint values from another party
        std::unordered_set<int> fingerprints_;
        // multimap from fingerprint to peeled value/key
        std::unordered_multimap<int, IndexType> my_fingerprints_;
        // new elements decoded in this round
        std::unordered_set<IndexType> new_elements_;
        // if k = 5, this stores pairs of elements with 3 collisions.
        std::set<std::pair<IndexType, IndexType>> pairs_5in3_;
        // a buffer that stores all affected neighbors in case of code counter updates.
        std::unordered_map<IndexType, ArrType> affected_neighbors_;   // <id, delta>
        UpdatePQ priority_queue_;
        DecodeConfig* config_;
        bool collision_resolving_;
        std::vector<std::pair<int, uint64_t>> unresolved_elements_;
        std::vector<IndexType> unresolved_ids_;
    };

}
