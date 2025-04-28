#pragma once
#include "libONIAK/oniakMath/overylarge.h"
#include <bitset>
#include <cstdint>
#include <functional>
#include <cstring>
#include <doro/sha256.hpp>

using IndexType = ONIAK::VeryLargeInt<256>;
using IndexTypeI64 = int64_t;

template <typename K> inline std::vector<uint8_t> to_binary_vector(const K &x) {
  constexpr size_t size = sizeof(K);
  std::vector<uint8_t> result(size);
  std::memcpy(result.data(), &x, size);
  return result;
}