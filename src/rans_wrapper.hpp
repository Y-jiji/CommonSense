#ifndef __RANS_WRAPPER_HPP__
#define __RANS_WRAPPER_HPP__

#include <cassert>
#include <cstdint>
#include <ranges>
#include <unordered_map>
#include <utility>
#include <vector>

#include <iostream>

#include "rans-tricks-master/rans_byte.h"

namespace Doro {

static const uint32_t prob_bits = 15, prob_scale = 1 << prob_bits;

// Only the RansWrapper can edit RansCode.
// Otherwise RansCode is read-only.
class RansCode {
public:
  RansCode(int size) : code_(size, 0), offset_(size - 2) {}
  size_t size() const { return code_.size() - offset_; }
  auto code_view() const { return std::views::all(code_) | std::views::drop(offset_); }

private:
  uint8_t* resize_if_full() {
    if (offset_ < 4) {  // result buffer is full
      size_t result_size = code_.size();
      size_t new_offset = offset_ + result_size;
      result_size *= 2;
      std::vector<uint8_t> result2(result_size, 0);
      std::copy(code_.data() + offset_, code_.data() + code_.size(), result2.data() + new_offset);
      code_ = std::move(result2);
      offset_ = new_offset;
    }
    return code_.data() + offset_;
  }
  uint8_t* data() { return code_.data() + offset_; }
  const uint8_t* data() const { return code_.data() + offset_; }
  const uint8_t* end() const { return code_.data() + code_.size(); }
  const uint8_t* begin() const { return code_.data(); }

  std::vector<uint8_t> code_;
  size_t offset_;
  template <typename T>
  friend class RansWrapper;
};

template <typename T> // type to be compressed
class RansWrapper {
public:
  explicit RansWrapper(const std::unordered_map<T, int>& frequencies) {
    std::vector<uint32_t> freqs, cum_freqs;
    uint32_t cum = 0;
    cum_freqs.push_back(0);
    for (auto& [sym, freq] : frequencies) {
      freqs.push_back(freq);
      cum += freq;
      cum_freqs.push_back(cum);
    }
    if (cum > prob_scale) {
      throw std::invalid_argument("frequency sum greater than 2^15.");
    }
    else if (frequencies.size() < 1) {
      throw std::invalid_argument("empty frequency map.");
    }

    for (auto& cfeq : cum_freqs) {
      cfeq = (cfeq * prob_scale) / cum;
    }
    for (size_t i = 0; i < freqs.size(); i++) {
      freqs[i] = cum_freqs[i + 1] - cum_freqs[i];
      assert(freqs[i] > 0);
    }
    assert(cum_freqs.back() == prob_scale);

    size_t i = 0, j = 0;
    // Build inverse cdf table
    inverse_cum_.resize(prob_scale);
    for (auto& [sym, freq] : frequencies) {
      RansEncSymbol symbol;
      RansEncSymbolInit(&symbol, cum_freqs[i], freqs[i], prob_bits);
      esyms_[sym] = symbol;
      RansDecSymbol dsymbol;
      RansDecSymbolInit(&dsymbol, cum_freqs[i], freqs[i]);
      dsyms_[sym] = dsymbol;
      ++i;
      while (j < cum_freqs[i]) {
        inverse_cum_[j] = sym;
        ++j;
      }
    }
  }

  RansCode encode(const std::vector<T>& data) const {
    RansState rans;
    RansEncInit(&rans);
    size_t result_size = data.size() * sizeof(T) + 2;
    RansCode result(result_size);
    uint8_t* ptr = result.data();
    for (auto& sym : data | std::views::reverse) {
      // out_of_range error if sym not in the frequency table in construction.
      RansEncPutSymbol(&rans, &ptr, &esyms_.at(sym));
      result.offset_ = ptr - result.begin();
      ptr = result.resize_if_full();
    }
    RansEncFlush(&rans, &ptr);
    result.offset_ = ptr - result.begin();
    return result;
  }

  std::vector<T> decode(const RansCode& code) const {
    RansState rans;
    const uint8_t* ptr = code.data();
    RansDecInit(&rans, &ptr);
    std::vector<T> result;
    while (true) {
      uint32_t value = RansDecGet(&rans, prob_bits);
      T sym = inverse_cum_.at(value);
      RansDecAdvanceSymbol(&rans, &ptr, &dsyms_.at(sym), prob_bits);
      if (ptr >= code.end()) {
        break;
      }
      result.push_back(sym);
    }
    return result;
  }

private:
  std::unordered_map<T, RansEncSymbol> esyms_;
  std::unordered_map<T, RansDecSymbol> dsyms_;
  std::vector<T> inverse_cum_;
};


} // namespace Doro
#endif