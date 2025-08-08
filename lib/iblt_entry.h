#pragma once

#include <cassert>
#include <cstring>
#include <doro/key_types.hpp>
#include <doro/sha256.hpp>
#include <inttypes.h>
#include <vector>

template <typename K>
class IBLTHashTableEntry {
public:
  int64_t count = 0;
  K keySum;
  uint64_t keyCheck = 0;
  std::vector<uint8_t> valueSum{};

  bool isPure() const;
  bool empty() const;
  void addValue(const std::vector<uint8_t> v);
};

template <typename K> bool IBLTHashTableEntry<K>::isPure() const {
  if (count == 1 || count == -1) {
    uint64_t check = HashWithSalt(to_binary_vector(keySum), N_HASHCHECK);
    return (keyCheck == check);
  }
  return false;
}

template <typename K> bool IBLTHashTableEntry<K>::empty() const {
  return (count == 0ull && keySum == K() && keyCheck == 0ull);
}

template <typename K>
void IBLTHashTableEntry<K>::addValue(const std::vector<uint8_t> v) {
  valueSum.resize(v.size());
  for (size_t i = 0; i < v.size(); i++) {
    valueSum[i] ^= v[i];
  }
}
