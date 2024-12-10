#ifndef __RANS_WRAPPER_HPP__
#define __RANS_WRAPPER_HPP__

#include <cassert>
#include <cstdint>
#include <ranges>
#include <unordered_map>
#include <utility>
#include <vector>

#include <iostream>

#include "rans-tricks/rans_byte.h"

namespace Doro {
// set to avoid errors in rANS for non-enough precision.
static const uint32_t max_symbol_freq = 1u << 15;
static const uint32_t prob_bits = 15;
static const uint32_t prob_scale = 1 << prob_bits;

// Only the RansWrapper can edit RansCode.
// Otherwise RansCode is read-only.
//
// In the default mode, the code_ only contains one symbol: the default one.
// All other symbols (presumbed to be rare) are in extra_.
class RansCode {
public:
  RansCode(int size) : code_(size, 0), offset_(size - 2), default_mode_(false), data_length_(0) {}
  // size is in terms of bits.
  size_t size() const {
    size_t code_size = code_.size() - offset_ - 2; // last two bytes are never used 
    return code_size * 8 + extra_index_.size() * 32 + extra_.size() * 8 + 1  // default_marker
      + (default_mode_ ? 32 : 0); // four bytes for data length in default mode.
  }
  // raw codes not containing extra data.
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

  std::vector<uint8_t> code_, extra_;
  // compressed data start from offset_ to code_.size()
  // extra data is the storage for data not in frequency table
  size_t offset_;
  std::vector<uint32_t> extra_index_;
  bool default_mode_;
  uint32_t data_length_; // length of data in default mode.
  // used only for debugging
  mutable std::vector<uint32_t> rans_states_;

  template <typename T, typename VAL>
  friend class RansWrapper;
};

template <typename T, typename VAL> // type to be compressed
class RansWrapper {
public:
  explicit RansWrapper(const std::unordered_map<T, VAL>& frequencies)
    : default_sym_(0), default_mode_(false), new_scale_bits_(prob_bits) {
    if (frequencies.size() < 1) {
      throw std::invalid_argument("empty frequency map.");
    }

    std::vector<uint32_t> freqs, cum_freqs;
    std::vector<double> real_cum;
    double cum = 0.0;
    real_cum.push_back(0.0);
    for (auto& [sym, freq] : frequencies) {
      cum += freq;
      real_cum.push_back(cum);
    }

    for (auto cfeq : real_cum) {
      cum_freqs.push_back((cfeq * prob_scale) / cum);
    }

    for (size_t i = 0; i < frequencies.size(); i++) {
      freqs.push_back(cum_freqs[i + 1] - cum_freqs[i]);
    }

    uint32_t scale_down = 1;
    uint32_t max_freqs = *std::max_element(freqs.begin(), freqs.end());
    if (max_freqs > max_symbol_freq) {
      while (max_freqs / scale_down > max_symbol_freq) {
        scale_down *= 2;
        --new_scale_bits_;
      }
      for (auto& cfreq : cum_freqs) {
        cfreq /= scale_down;
      }
      freqs.clear();
      for (size_t i = 0; i < frequencies.size(); i++) {
        freqs.push_back(cum_freqs[i + 1] - cum_freqs[i]);
      }
    }
    uint32_t new_prob_scale = prob_scale / scale_down;
    assert(cum_freqs.back() == new_prob_scale);
    size_t i = 0, j = 0;
    // Build inverse cdf table
    inverse_cum_.resize(new_prob_scale);
    for (auto& [sym, freq] : frequencies) {
      if (freqs[i] == 0) {
        ++i;
        continue;
      }
      default_sym_ = sym;  // in default mode, default_sym_ will be the only symbol with nonzero frequency.
      RansEncSymbol symbol;
      RansEncSymbolInit(&symbol, cum_freqs[i], freqs[i], new_scale_bits_);
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
    if (esyms_.size() == 1) default_mode_ = true;
  }

  RansCode encode(const std::vector<T>& data) const {
    RansState rans;
    RansEncInit(&rans);
    size_t result_size = default_mode_ ? sizeof(T) : data.size() * sizeof(T) + 2;
    RansCode result(result_size);
    if (default_mode_) {
      auto result_head = reinterpret_cast<T*>(result.data());
      *result_head = default_sym_;
      result.offset_ = 0;
      result.default_mode_ = true;
      result.data_length_ = data.size();
    }

    uint8_t* ptr = result.data();
    for (int i = data.size() - 1; i >= 0; --i) {
      const auto& sym = data[i];
      // if symbol is not in frequency table
      // In default mode, the frequency table only contains the default symbol.
      if (!esyms_.contains(sym)) {
        result.extra_index_.push_back(i);
        auto sym_raw = reinterpret_cast<const uint8_t*>(&sym);
        result.extra_.insert(result.extra_.end(), sym_raw, sym_raw + sizeof(T));
        // space to store this data and its index
        continue;
      }
      // In default mode, the default symbol is not encoded.
      if (default_mode_) continue;
#ifdef ONIAK_DEBUG
      result.rans_states_.push_back(rans);
#endif

      RansEncPutSymbol(&rans, &ptr, &esyms_.at(sym));
      result.offset_ = ptr - result.begin();
      ptr = result.resize_if_full();
    }
    if (!default_mode_) {
      RansEncFlush(&rans, &ptr);
      result.offset_ = ptr - result.begin();
    }

    return result;
  }

  std::vector<T> decode(const RansCode& code) const {
    RansState rans;
    const uint8_t* ptr = code.data();
    RansDecInit(&rans, const_cast<uint8_t**>(&ptr));
    std::vector<T> result;
    auto extra_index = code.extra_index_.end();
    const uint8_t* extra_ptr = code.extra_.data() + code.extra_.size();
    while (true) {
      while (extra_index > code.extra_index_.begin() && *(extra_index - 1) == result.size()) {
        extra_index--;
        extra_ptr -= sizeof(T);
        // symbol in extra list.
        result.push_back(*reinterpret_cast<const T*>(extra_ptr));
      }
      T sym;
      if (code.default_mode_) {
        sym = *reinterpret_cast<const T*>(code.data());
        if (result.size() >= code.data_length_) {
          break;
        }
      }
      else {
        uint32_t value = RansDecGet(&rans, new_scale_bits_);
        sym = inverse_cum_.at(value);
        RansDecAdvanceSymbol(&rans, const_cast<uint8_t**>(&ptr), &dsyms_.at(sym), new_scale_bits_);

#ifdef ONIAK_DEBUG
        if (!code.rans_states_.empty()) {
          assert(rans == code.rans_states_.back());
          code.rans_states_.pop_back();
        }
#endif

        if (ptr >= code.end()) {
          break;
        }
      }
      result.push_back(sym);
    }
    return result;
  }

private:
  std::unordered_map<T, RansEncSymbol> esyms_;
  std::unordered_map<T, RansDecSymbol> dsyms_;
  std::vector<T> inverse_cum_;
  T default_sym_;
  bool default_mode_;
  int new_scale_bits_;
};


} // namespace Doro
#endif