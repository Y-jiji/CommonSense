#ifndef __DORO_DECODER_HPP__
#define __DORO_DECODER_HPP__

#include "doro.hpp"
#include "oniakDataStructure/oupq.h"
#include "oniakHash/ohash.h"

#include <algorithm>
#include <cassert>
#include <format>
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
  L2
};
struct DecodeConfig {
  int tk, max_stage; // used for stage decode.
  // tk is stage size
  int ta; // terminate if any element thrashes for ta times.
  bool verbose, debug;
  // lower and upper bound of value range
  double lb, ub;
  PursuitChoice pursuit_choice;
};

template <typename ArrType = int32_t, ONIAK::UpdatePQBackend backend = ONIAK::UpdatePQBackend::PriorityQueue>
class DoroDecoder {
public:
  using DoroCodeT = DoroCode<ArrType>;
  using TwoDimVector = std::vector<std::vector<int>>;
  using UpdatePQ = ONIAK::UpdatePQAdapter<int, ArrType, backend>::type;

  DoroDecoder(ONIAK::WYHash* finger_hash = nullptr) : finger_hash_(finger_hash), priority_queue_(
    // function from sensed strength to priority
    [](ArrType a) {return std::abs(a);}
  ), collision_resolving_(false) {
  }

  int stage_decode(
    DoroCodeT* code,
    const std::vector<int>& setA, /*possible candidates of code*/
    const DecodeConfig* config,
    std::unordered_map<int, ArrType>*& result) {
    code_ = code;
    config_ = config;
    thrashing_.clear();
    result = &result_;
    bool is_l2 = config_->pursuit_choice == PursuitChoice::L2; // else is l1
    for (int element : setA) {
      thrashing_[element] = 0;
    }

    for (int rnd : std::views::iota(0, config_->max_stage)) {
      std::vector<std::tuple<ArrType, int, ArrType>> signals;
      for (auto element_thrash : thrashing_) {
        int element = element_thrash.first;
        ArrType strength = is_l2 ? code_->sense(element) : code_->sense_l1(element);
        ArrType delta = strength_to_delta(element, strength);
        if (std::abs(delta) > 1e-6)
          signals.emplace_back(-std::abs(strength), element, delta);
      }
      std::nth_element(signals.begin(), signals.begin() + config_->tk, signals.end());
      if (signals.size() > static_cast<size_t>(config_->tk))
        signals.resize(config_->tk);
      bool finished = true;
      for (auto [power, element, delta] : signals) {
        finished = false;
        int thrash = thrashing_.find(element)->second;
        if (thrash > config_->ta) {
          finished = true;
          break;
        }
        thrashing_[element] = thrash + 1;
        code_->peel(element, delta);
        auto cur_iter = result_.find(element);
        ArrType cur_element_value = (cur_iter != result_.end()) ? cur_iter->second : 0;
        result_[element] = cur_element_value + delta;

        if (config_->verbose) {
          std::cout << std::format("thrash: {}, power: {}, element: {}, cur_element_value: {}, delta: {}\n",
            thrash, power, element, cur_element_value, delta);
        }
      }
      if (config_->verbose) {
        std::cout << std::format("Round: {}\n", rnd);
        code_->show_result();
      }
      if (finished)
        return rnd;
    }
    return config_->max_stage;
  }

  int decode_resense(
    DoroCodeT* code,
    const std::vector<int>& setA, /*possible candidates of code*/
    const DecodeConfig* config,
    std::unordered_map<int, ArrType>*& result) {
    code_ = code;
    config_ = config;
    thrashing_.clear();
    result = &result_;
    neighbors_.assign(code_->size(), {});
    priority_queue_.clear();
    bool is_l2 = config_->pursuit_choice == PursuitChoice::L2; // else is l1

    for (int element : setA) {
      for (int i : std::views::iota(0, code_->k())) {
        auto [index, sign] = code_->hash(i, element);
        neighbors_[index].push_back(element);
      }
      ArrType strength = is_l2 ? code_->sense(element) : code_->sense_l1(element);
      ArrType delta = strength_to_delta(element, strength);

      if (std::abs(delta) > 1e-6) {
        priority_queue_.push(element, strength);
      }
      else
        priority_queue_.set(element, strength);
      // priority_queue_.push(element, code_->sense(element));
    }

    while (!priority_queue_.empty()) {
      auto [strength, element] = priority_queue_.top();
      auto thrash_iter = thrashing_.find(element);
      int thrash = (thrash_iter != thrashing_.end()) ? thrash_iter->second : 0;
      if (thrash > config_->ta) break;
      thrashing_[element] = thrash + 1;

      ArrType delta = strength_to_delta(element, strength);
      code_->peel(element, delta);
      auto cur_iter = result_.find(element);
      ArrType cur_element_value = (cur_iter != result_.end()) ? cur_iter->second : 0;
      result_[element] = cur_element_value + delta;
      priority_queue_.pop();
      priority_queue_.set(element, new_strength(element, strength, delta, code_->k()));

      std::vector<int> affected_neighbors;
      for (int i : std::views::iota(0, code_->k())) {
        auto [index, sign] = code_->hash(i, element);
        affected_neighbors.insert(affected_neighbors.end(), neighbors_[index].begin(), neighbors_[index].end());
      }
      // There may be duplicates.
      std::sort(affected_neighbors.begin(), affected_neighbors.end());
      for (size_t i : std::views::iota(0u, affected_neighbors.size())) {
        int neighbor = affected_neighbors[i];
        int neighbor_cnt = 1;
        if (neighbor == element)  continue;
        while (i < affected_neighbors.size() - 1 && neighbor == affected_neighbors[i + 1]) {
          ++neighbor_cnt;
          ++i;
        }
        ArrType neighbor_strength = new_strength(neighbor, priority_queue_[neighbor], delta, neighbor_cnt);
        ArrType neighbor_delta = strength_to_delta(neighbor, neighbor_strength);
        if (std::abs(neighbor_delta) > 1e-6)
          priority_queue_.update(neighbor, neighbor_strength);
        else
          priority_queue_.set(neighbor, neighbor_strength);
      }
      if (config_->verbose) {

        std::cout << std::format("thrash: {}, power: {}, element: {}, cur_element_value: {}, delta: {}, sl: {}\n",
          thrash, strength, element, cur_element_value, delta, priority_queue_.size());
        code_->show_result();
      }
    }
    return code_->num_peels();
  }

  int decode(
    DoroCodeT* code,
    const std::vector<int>& setA, /*possible candidates of code*/
    const DecodeConfig* config,
    std::unordered_map<int, ArrType>*& result) {
    code_ = code;
    config_ = config;
    thrashing_.clear();
    result = &result_;
    neighbors_.assign(code_->size(), {});
    neighbors2_.assign(code_->size(), {});
    priority_queue_.clear();
    bool is_l2 = config_->pursuit_choice == PursuitChoice::L2; // else is l1

    for (int element : setA) {
      for (int i : std::views::iota(0, code_->k())) {
        auto [index, sign] = code_->hash(i, element);
        if (sign > 0)
          neighbors_[index].push_back(element);
        else neighbors2_[index].push_back(element);
      }
      ArrType strength = is_l2 ? code_->sense(element) : code_->sense_l1(element);
      ArrType delta = strength_to_delta(element, strength);

      // if fingerprint mechanism is active, then this element must not collide with any known fingerprint
      bool fingerprint_condition = (finger_hash_ == nullptr) || !fingerprints_.contains((*finger_hash_)(element));
      if (std::abs(delta) > 1e-6 && fingerprint_condition) {
        priority_queue_.push(element, strength);
      }
      else
        priority_queue_.set(element, strength);
      // priority_queue_.push(element, code_->sense(element));
    }

    while (!priority_queue_.empty()) {
      auto [strength, element] = priority_queue_.top();
      auto thrash_iter = thrashing_.find(element);
      int thrash = (thrash_iter != thrashing_.end()) ? thrash_iter->second : 0;
      if (thrash > config_->ta) break;
      thrashing_[element] = thrash + 1;

      ArrType delta = strength_to_delta(element, strength);
      code_->peel(element, delta);
      auto cur_iter = result_.find(element);
      ArrType cur_element_value = (cur_iter != result_.end()) ? cur_iter->second : 0;
      result_[element] = cur_element_value + delta;
      priority_queue_.pop();
      priority_queue_.set(element, new_strength2(element, strength, delta, code_->k()));

      std::unordered_map<int, int> affected_neighbors;
      for (int i : std::views::iota(0, code_->k())) {
        auto [index, sign] = code_->hash(i, element);
        std::vector<int>& pos_neighbors = (sign > 0) ? neighbors_[index] : neighbors2_[index];
        std::vector<int>& neg_neighbors = (sign < 0) ? neighbors_[index] : neighbors2_[index];
        for (int neighbor : pos_neighbors) {
          if (neighbor == element) continue;
          auto neighbor_iter = affected_neighbors.find(neighbor);
          affected_neighbors[neighbor] = (neighbor_iter != affected_neighbors.end()) ? neighbor_iter->second + 1 : 1;
        }
        for (int neighbor : neg_neighbors) {
          if (neighbor == element) continue;
          auto neighbor_iter = affected_neighbors.find(neighbor);
          affected_neighbors[neighbor] = (neighbor_iter != affected_neighbors.end()) ? neighbor_iter->second - 1 : -1;
        }
      }
      for (auto [neighbor, neighbor_cnt] : affected_neighbors) {
        ArrType neighbor_strength = new_strength2(neighbor, priority_queue_[neighbor], delta, neighbor_cnt);
        ArrType neighbor_delta = strength_to_delta(neighbor, neighbor_strength);
        bool fingerprint_condition = (finger_hash_ == nullptr) || !fingerprints_.contains((*finger_hash_)(element));
        if (std::abs(neighbor_delta) > 1e-6 && fingerprint_condition)
          priority_queue_.update(neighbor, neighbor_strength);
        else
          priority_queue_.set(neighbor, neighbor_strength);
      }
      if (config_->verbose) {
        std::cout << std::format("thrash: {}, power: {}, element: {}, cur_element_value: {}, delta: {}\n",
          thrash, strength, element, cur_element_value, delta);
        code_->show_result();
      }
    }
    return code_->num_peels();
  }

  void update_fingerprints(DoroDecoder& other) const {
    assert(finger_hash_ != nullptr && *other.finger_hash_ == *finger_hash_);
    for (auto [key, value] : result_) {
      if (std::abs(value) > 1e-6)
        other.fingerprints_.insert((*finger_hash_)(key));
    }
  }

  void reset() {
    result_.clear();
    fingerprints_.clear();
    thrashing_.clear();
    priority_queue_.clear();
    collision_resolving_ = false;
  }

  void enter_resolving() {
    collision_resolving_ = true;
  }

  std::unordered_set<int>& fingerprints() {
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

  ArrType strength_to_delta(int element, ArrType strength) const {
    bool is_l2 = config_->pursuit_choice == PursuitChoice::L2;
    auto cur_iter = result_.find(element);
    ArrType cur_element_value = (cur_iter != result_.end()) ? cur_iter->second : 0;
    // In L2, delta is computed from average instead of sum
    double normalized_strength = is_l2 ? static_cast<double>(strength) / static_cast<double>(code_->k()) : strength;
    return get_delta(normalized_strength, cur_element_value);
  }

  ArrType new_strength(int element, ArrType strength, ArrType delta, int k) const {
    bool is_l2 = config_->pursuit_choice == PursuitChoice::L2;
    return is_l2 ? code_->sense(element) : code_->sense_l1(element);
  }

  ArrType new_strength2(int element, ArrType strength, ArrType delta, int k) const {
    bool is_l2 = config_->pursuit_choice == PursuitChoice::L2;
    return is_l2 ? strength - delta * k : code_->sense_l1(element);
  }

  DoroCodeT* code_;
  ONIAK::WYHash* finger_hash_;
  TwoDimVector neighbors_, neighbors2_;
  std::unordered_map<int, ArrType> result_;
  std::unordered_set<int> fingerprints_;
  std::unordered_map<int, int> thrashing_;
  UpdatePQ priority_queue_;
  const DecodeConfig* config_;
  bool collision_resolving_;
};
}
#endif