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
    if (!bch_) {
      throw std::runtime_error("Failed to initialize BCH codec.");
    }
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
    std::vector<uint8_t> parity_bits, data_copy = data;
    if (bit_by_bit) {
      size_t number_blocks = (data.size() + message_size() - 1) / message_size();
      data_copy.resize(number_blocks * message_size(), 0);
      parity_bits.assign(number_blocks * parity_bits_size(), 0);
      for (size_t b = 0; b < number_blocks; ++b) {
        encodebits_bch(bch_, data_copy.data() + b * message_size(), parity_bits.data() + b * parity_bits_size());
      }
    }
    else {
      parity_bits.assign(bch_->ecc_bytes, 0);
      assert(data.size() == message_size_bytes());
      encode_bch(bch_, data.data(), data.size(), parity_bits.data());
    }
    return parity_bits;
  }

  int encoded_parity_size(const std::vector<uint8_t>& data) const {
    size_t number_blocks = (data.size() + message_size() - 1) / message_size();
    return number_blocks * parity_bits_size();
  }

  std::vector<uint8_t> decode(const std::vector<uint8_t>& data,
    const std::vector<uint8_t>& parity_bits, bool bit_by_bit = false) const {
    std::vector<unsigned int> errLocOut(bch_->t);
    std::vector<uint8_t> decoded_data(data);
    int nerrFound;
    if (bit_by_bit) {
      size_t number_blocks = (data.size() + message_size() - 1) / message_size();
      decoded_data.resize(number_blocks * message_size(), 0);
      assert(parity_bits.size() == parity_bits_size() * number_blocks);
      for (size_t b = 0; b < number_blocks; ++b) {
        nerrFound = decodebits_bch(bch_, decoded_data.data() + b * message_size(),
          parity_bits.data() + b * parity_bits_size(), errLocOut.data());
        if (nerrFound < 0) { // error correction failed
          continue;
        }
        else {
          correctbits_bch(bch_, decoded_data.data() + b * message_size(), errLocOut.data(), nerrFound);
        }
      }
      decoded_data.resize(data.size());
    }
    else {
      assert(parity_bits.size() == bch_->ecc_bytes);
      assert(data.size() == message_size_bytes());
      nerrFound = decode_bch(bch_, data.data(), data.size(), parity_bits.data(), NULL, NULL, errLocOut.data());
      if (nerrFound < 0) { // error correction failed
        return data;
      }
      else {
        correct_bch(bch_, decoded_data.data(), data.size(), errLocOut.data(), nerrFound);
      }
    }
    return decoded_data;
  }

  // num_failed_blocks and num_corrected_errors are incremented.
  std::vector<size_t> get_error_positions(const std::vector<uint8_t>& data,
    const std::vector<uint8_t>& parity_bits, bool bit_by_bit=false, int* num_failed_blocks=nullptr,
    int* num_corrected_errors=nullptr) const {
    assert(bit_by_bit == true); // Otherwise, not implemented.
    std::vector<unsigned int> errLocOut(bch_->t);
    std::vector<uint8_t> decoded_data = data;
    std::vector<size_t> error_pos;
    int nerrFound;

    if (bit_by_bit) {
      size_t number_blocks = (data.size() + message_size() - 1) / message_size();
      decoded_data.resize(number_blocks * message_size(), 0);
      assert(parity_bits.size() == parity_bits_size() * number_blocks);
      for (size_t b = 0; b < number_blocks; ++b) {
        nerrFound = decodebits_bch(bch_, decoded_data.data() + b * message_size(),
          parity_bits.data() + b * parity_bits_size(), errLocOut.data());
        if (nerrFound < 0) { // error correction failed
          if (num_failed_blocks) {
            (*num_failed_blocks)++;
          }
          continue;
        }
        else {
          if (num_corrected_errors) {
            *num_corrected_errors += nerrFound;
          }
          for (int i = 0; i < nerrFound; ++i) {
            error_pos.push_back(b * message_size() + errLocOut[i]);
          }
        }
      }
    }
    return error_pos;
  }

private:
  bch_control* bch_;
};

}

#endif