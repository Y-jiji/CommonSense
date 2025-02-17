#ifndef __DORO_HPP__
#define __DORO_HPP__

#include "libONIAK/oniakHash/ohash.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <inttypes.h>
#include <iostream>
#include <numeric>
#include <ranges>
#include <unordered_map>
#include <vector>

namespace Doro {

constexpr int kDoroNumHashFuncs = 30;
constexpr int kDoroNumHashFuncsMargin = 10;

template <typename IndexType, typename ArrType = int32_t>
class DoroCode {
public:
  using SparseVector = std::unordered_map<IndexType, ArrType>;
  using PairType = std::pair<IndexType, ArrType>;
  using IndexT = IndexType;
  using ArrayT = ArrType;

  template <typename RandomDevice>
  DoroCode(int size, int k, bool is_cbf, RandomDevice& rng, int lb = 0, int ub = 0) :
    arr_(size), is_cbf_(is_cbf), k_(k), num_peels_(0), num_correct_peels_(0),
    num_recenters_(0), lb_(lb), ub_(ub), interval_(ub - lb) {
    int num_hash_funcs = std::max(kDoroNumHashFuncs, k + kDoroNumHashFuncsMargin);
    hash_funcs_.reserve(num_hash_funcs);
    std::ranges::for_each(std::views::iota(0, num_hash_funcs),
      [&](int) {hash_funcs_.emplace_back(/*mask*/ 0, /*mod*/ 2 * size, /*seed*/ rng());});
  }

  void encode(const SparseVector& kvpairs) {
    assert(values_.empty());
    for (auto [x, val] : kvpairs) {
      auto all_hashes = hash_all(x);
      for (auto [index, sign] : all_hashes) {
        arr_[index] += sign * val;
      }
      values_[x] = val;
    }
  }

  // if delta is 1, all counters will be REDUCED by 1 (if positive sign)
  // or INCREASED by 1 (if negative sign)
  void peel(IndexType element, ArrType delta) {
    num_peels_ += 1;
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      arr_[index] -= sign * delta;
    }
    auto element_iter = values_.find(element);
    ArrType ori_value = (element_iter != values_.end()) ? element_iter->second : 0;
    ArrType new_value = ori_value - delta;
    values_[element] = new_value;
    if (std::abs(new_value) < std::abs(ori_value))
      num_correct_peels_ += 1;
  }

  // senses the signal of an element by inner product
  // was performance bottleneck
  ArrType sense(IndexType element) const {
    ArrType signal = 0;
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      signal += sign * arr_[index];
    }
    return signal;
  }

  void print_key(IndexType element) const {
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      std::print(std::cout, "{}\\{} ", index, sign * arr_[index]);
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

  ArrType deviation(ArrType x) const {
    if (interval_ <= 0 || (x >= lb_ && x < ub_)) return 0;
    return std::max(lb_ - x, x - ub_ + 1);
  }

  // Returns median signal. If k_ is even, returns the smaller one in absolute value.
  ArrType sense_l1(IndexType element) const {
    std::vector<ArrType> signals;
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      signals.push_back(sign * arr_[index]);
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

  SparseVector& values() { return values_; }

  // each value is an (index, sign) pair
  std::vector<std::pair<int, int>> hash_all(IndexType x) const {
    std::vector<std::pair<int, int>> all_hashes;
    all_hashes.reserve(k_);
    for (auto& hfunc : hash_funcs_) {
      int hash_value = hfunc.hash_in_range(x);
      int idx = (hash_value >= size()) ? hash_value - size() : hash_value;
      int sign = (is_cbf_ || hash_value < size()) ? 1 : -1;
      if (!std::ranges::contains(all_hashes, idx, [](auto x) {return x.first;}))
        all_hashes.emplace_back(idx, sign);
      // otherwise reject this sample.
      if (all_hashes.size() == static_cast<size_t>(k_)) return all_hashes;
    }
    throw std::runtime_error("Error: Run out of hash functions.");
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

void foo();

#endif