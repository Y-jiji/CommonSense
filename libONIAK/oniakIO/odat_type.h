#pragma once

#include <cstdint>

namespace ONIAK {

template <typename T>
struct ODAT_code {};

template <>
struct ODAT_code<float> {
  static constexpr unsigned char value = 0x52;
};

template <>
struct ODAT_code<double> {
  static constexpr unsigned char value = 0x62;
};
template <>
struct ODAT_code<int32_t> {
  static constexpr unsigned char value = 0x51;
};
template <>
struct ODAT_code<int64_t> {
  static constexpr unsigned char value = 0x61;
};
template <>
struct ODAT_code<int8_t> {
  static constexpr unsigned char value = 0x31;
};
template <>
struct ODAT_code<int16_t> {
  static constexpr unsigned char value = 0x41;
};
template <>
struct ODAT_code<uint8_t> {
  static constexpr unsigned char value = 0x33;
};

};


 