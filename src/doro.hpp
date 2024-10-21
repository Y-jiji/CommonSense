#ifndef __DORO_HPP__
#define __DORO_HPP__

#include "oniakHash/ohash.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <inttypes.h>
#include <iostream>
#include <numeric>
#include <ranges>
#include <unordered_map>
#include <vector>

namespace Doro {
template <typename ArrType = int32_t>
class DoroCode {
public:
  using IndexType = int;
  using SparseVector = std::unordered_map<IndexType, ArrType>;
  using PairType = std::pair<IndexType, ArrType>;

  template <typename RandomDevice>
  DoroCode(int size, int k, bool is_cbf, RandomDevice& rng, int lb = 0, int ub=0) :
   arr_(size), is_cbf_(is_cbf), k_(k), num_peels_(0), num_correct_peels_(0),
    num_recenters_(0), lb_(lb), ub_(ub), interval_(ub-lb) {
    hash_funcs_.reserve(k_);
    for ([[maybe_unused]] int i : std::views::iota(0, k_)) {
      hash_funcs_.emplace_back(/*mask*/ 0, /*mod*/ 2 * size, /*seed*/ rng());
    }
  }

  void encode(const SparseVector& kvpairs) {
    assert(values_.empty());
    for (auto [x, val] : kvpairs) {
      for (int i : std::views::iota(0, k_)) {
        auto [index, sign] = hash(i, x);
        arr_[index] += sign * val;
      }
      values_[x] = val;
    }
  }

  // if delta is 1, all counters will be REDUCED by 1 (if positive sign)
  // or INCREASED by 1 (if negative sign)
  void peel(int element, ArrType delta) {
    num_peels_ += 1;
    for (int i : std::views::iota(0, k_)) {
      auto [index, sign] = hash(i, element);
      arr_[index] -= sign * delta;
    }
    auto element_iter = values_.find(element);
    int ori_value = (element_iter != values_.end()) ? element_iter->second : 0;
    int new_value = ori_value - delta;
    values_[element] = new_value;
    if (std::abs(new_value) < std::abs(ori_value))
      num_correct_peels_ += 1;
  }

  // senses the signal of an element by inner product
  // performance bottleneck
  ArrType sense(int element) const {
    ArrType signal = 0;
    for (int i : std::views::iota(0, k_)) {
      auto [index, sign] = hash(i, element);
      signal += sign * arr_[index];
    }
    return signal;
  }

  void print_key(int element) const {
    for (int i : std::views::iota(0, k_)) {
      auto [index, sign] = hash(i, element);
      std::cout << index << "\\" << sign * arr_[index] << " ";
    }
    std::cout << std::endl;
  }

  // recenter a counter so that it lies in [lb, ub)
  ArrType recenter(ArrType x) {
    if (interval_ <= 0) return x;   // never recenters
    x -= std::floor(static_cast<float>(x - lb_) / interval_) * interval_;
    assert(x >= lb_ && x < ub_);
    return x;
  }

  // Returns median signal. If k_ is even, returns the smaller one in absolute value.
  ArrType sense_l1(int element) const {
    std::vector<ArrType> signals(k_);
    for (int i : std::views::iota(0, k_)) {
      auto [index, sign] = hash(i, element);
      signals[i] = sign * arr_[index];
    }
    auto med_iter = signals.begin() + k_ / 2;
    std::nth_element(signals.begin(), med_iter, signals.end());
    if (k_ % 2 == 0) {
      auto med_iter2 = med_iter - 1;
      std::nth_element(signals.begin(), med_iter2, med_iter);
      if (*med_iter * (*med_iter2) < 0)
        return 0;
      else if (std::abs(*med_iter) > std::abs(*med_iter2))
        return *med_iter2;
    }
    return *med_iter;
  }

  int nonzero_num() const {
    return std::count_if(values_.begin(), values_.end(), [](PairType pair) { return pair.second != 0; });
  }

  ArrType mae() const {
    return std::accumulate(values_.begin(), values_.end(), 0, [](ArrType sum, PairType pair) { return sum + std::abs(pair.second); });
  }

  // if any two hashes do not collide, return 0
  // if there are two hashes colliding under the same sign, return 1
  // if colliding under opposite signs, return 2
  int is_colliding(int element) const {
    std::unordered_map<int, int> hash_to_sign;
    for (auto i : std::views::iota(0, k_)) {
      auto [index, sign] = hash(i, element);
      auto iter = hash_to_sign.find(index);
      if (iter == hash_to_sign.end()) {
        hash_to_sign[index] = sign;
      } else if (iter->second != sign) {
        return 2;
      } else {
        return 1;
      }
    }
    return 0;
  }

  void reset() {
    arr_.assign(size(), 0);
    num_peels_ = 0;
    num_correct_peels_ = 0;
    num_recenters_ = 0;
    values_ = {};
  }

  int size() const {
    return arr_.size();
  }

  std::vector<ArrType>& code() { return arr_; }
  const std::vector<ArrType>& code() const { return arr_; }

  int k() const { return k_; }
  int num_peels() const { return num_peels_; }
  int& num_peels() { return num_peels_; }
  int& num_correct_peels() { return num_correct_peels_; }
  int& num_recenters() { return num_recenters_; }
  ArrType interval() const { return interval_; }
  double mid_point() const { return (lb_ + ub_ - 1.0) / 2.0; }

  void show_result() const {
    std::cout << "Number of peels: " << num_peels_ << std::endl;
    std::cout << "Number of correct peels: " << num_correct_peels_ << std::endl;
    std::cout << "Number of wrong peels: " << num_peels_ - num_correct_peels_ << std::endl;
    std::cout << "Size of nonzero indices: " << nonzero_num() << std::endl;
    std::cout << "L1 Norm of indices: " << mae() << std::endl;
  }

  std::pair<IndexType, int> hash(int i, IndexType x) const {
    int hash_value = hash_funcs_[i].hash_in_range(x);
    int idx = (hash_value >= size()) ? hash_value - size() : hash_value;
    int sign = (is_cbf_ || hash_value < size()) ? 1 : -1;
    return { idx, sign };
  }

  bool is_cbf() const { return is_cbf_; }

  bool empty() const {
    return std::all_of(arr_.begin(), arr_.end(), [](ArrType val) { return val == 0; });
  }

private:
  std::vector<ArrType> arr_;
  SparseVector values_;

  // If is_cbf_ is true, then all hashed signs are positive.
  bool is_cbf_; // cbf: counting Bloom filter
  std::vector<ONIAK::WYHash> hash_funcs_;
  int k_, num_peels_, num_correct_peels_, num_recenters_;
  // [lower bounds, and upper bounds) in quantization
  // if lb >= ub, then quantization is not in effect.
  ArrType lb_, ub_;
  float interval_; // ub - lb
};

} // namespace Doro

#endif