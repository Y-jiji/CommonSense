#ifndef __BCH_WRAPPER_HPP
#define __BCH_WRAPPER_HPP

#include "bch_codec/bch_codec.h"

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace Doro {

class BCHWrapper {
public:
  // if primitive_polynomial is 0, the default one will be used.
  BCHWrapper(int order, int capability, unsigned primitive_polynomial = 0)
    : bch_(init_bch(order, capability, primitive_polynomial)) {
  }
  // forbid copy and move
  BCHWrapper(const BCHWrapper&) = delete;
  BCHWrapper& operator=(const BCHWrapper&) = delete;
  BCHWrapper(BCHWrapper&&) = delete;
  BCHWrapper& operator=(BCHWrapper&&) = delete;
  // Destructor to handle c-ctyle ownership
  ~BCHWrapper() { free_bch(bch_); }

  int codeword_size() const { return bch_->n; }
  int order() const { return bch_->m; }
  int capability() const { return bch_->t; }
  int parity_bits_size() const { return bch_->ecc_bits; }
  size_t message_size() const { return bch_->n - bch_->ecc_bits; }
  size_t message_size_bytes() const { return message_size() / 8; }

  // if bit_by_bit, each element of input data is one bit. Only the LSB is used.
  std::vector<uint8_t> encode(const std::vector<uint8_t>& data, bool bit_by_bit = false) const {
    std::vector<uint8_t> parity_bits(bch_->ecc_bytes, 0);
    if (bit_by_bit) {
      assert(data.size() == message_size());
      encodebits_bch(bch_, data.data(), parity_bits.data());
    }
    else {
      assert(data.size() == message_size_bytes());
      encode_bch(bch_, data.data(), data.size(), parity_bits.data());
    }
    return parity_bits;
  }

  std::vector<uint8_t> decode(const std::vector<uint8_t>& data,
    const std::vector<uint8_t>& parity_bits, bool bit_by_bit = false) const {
    std::vector<unsigned int> errLocOut(bch_->t);
    int nerrFound;
    assert(parity_bits.size() == bch_->ecc_bytes);
    if (bit_by_bit) {
      assert(data.size() == message_size());
      nerrFound = decodebits_bch(bch_, data.data(), parity_bits.data(), errLocOut.data());
    } else {
      assert(data.size() == message_size_bytes());
      nerrFound = decode_bch(bch_, data.data(), data.size(), parity_bits.data(), NULL, NULL, errLocOut.data());
    }

    if (nerrFound < 0) {
      // error correction failed
      return data;
    }

    std::vector<uint8_t> decoded_data(data);
    if (bit_by_bit) {
      correctbits_bch(bch_, decoded_data.data(), errLocOut.data(), nerrFound);
    } else {
      correct_bch(bch_, decoded_data.data(), data.size(), errLocOut.data(), nerrFound);
    }
    return decoded_data;
  }

private:
  bch_control* bch_;
};

}

#endif