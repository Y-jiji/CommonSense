#pragma once

#include <cassert>
#include <cstring>
#include <doro/iblt_entry.h>
#include <doro/key_types.hpp>
#include <doro/sha256.hpp>
#include <inttypes.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

template<typename K>
class IBLT {
public:

  IBLT(size_t _expectedNumEntries, size_t _ValueSize, float hedge = 4.0, size_t _numHashes = 4);
  IBLT(const IBLT<K>& other);
  virtual ~IBLT();

  void insert(const K& k, const std::vector<uint8_t> v);
  void erase(const K& k, const std::vector<uint8_t> v);

  bool listEntries(
    std::unordered_map<K, std::vector<uint8_t>>& positive,
    std::unordered_map<K, std::vector<uint8_t>>& negative,
    const std::unordered_set<K>* positive_supset = nullptr,
    const std::unordered_set<K>* negative_supset = nullptr
  ) const;

  // Subtract two IBLTs
  IBLT<K> operator-(const IBLT<K>& other) const;

  int hashTableSize();

  // these need to be public for python integration 
  size_t valueSize;
  size_t numHashes;

private:
  void insertImpl(int plusOrMinus, const K k, const std::vector<uint8_t> v);
  std::vector<IBLTHashTableEntry<K>> hashTable;
};

template <typename K>
int IBLT<K>::hashTableSize() { return (hashTable.size()); }

template <typename K>
IBLT<K>::IBLT(size_t _expectedNumEntries, size_t _valueSize, float hedge,
  size_t _numHashes)
  : valueSize(_valueSize), numHashes(_numHashes) {
  assert(numHashes > 1);
  size_t nEntries = (int)(((float)_expectedNumEntries) * hedge);
  nEntries = (nEntries + numHashes - 1) / numHashes * numHashes;
  hashTable.resize(nEntries);
}

template <typename K>
IBLT<K>::IBLT(const IBLT<K>& other) {
  valueSize = other.valueSize;
  hashTable = other.hashTable;
  numHashes = other.numHashes;
}

template <typename K>
IBLT<K>::~IBLT() {}

template <typename K>
void IBLT<K>::insertImpl(int plusOrMinus, const K key, const std::vector<uint8_t> v) {
  assert(v.size() == valueSize);
  assert(key != K());
  assert(plusOrMinus == 1 || plusOrMinus == -1);
  std::unordered_set<uint64_t> hashPositions;
  for (size_t i = 0; i < numHashes; i++) {
    uint64_t h = HashWithSalt(to_binary_vector(key), i) % hashTable.size();
    if (hashPositions.contains(h))
      continue;
    hashPositions.insert(h);
    IBLTHashTableEntry<K>& entry = hashTable.at(h);
    entry.count += plusOrMinus;
    entry.keySum ^= key;
    entry.keyCheck ^= HashWithSalt(to_binary_vector(key), N_HASHCHECK);
    if (entry.empty()) {
      entry.valueSum.clear();
    } else {
      entry.addValue(v);
    }
  }
}

template <typename K>
void IBLT<K>::insert(const K& k, const std::vector<uint8_t> v) {
  insertImpl(1, k, v);
}

template <typename K>
void IBLT<K>::erase(const K& k, const std::vector<uint8_t> v) {
  insertImpl(-1, k, v);
}

template <typename K>
bool IBLT<K>::listEntries(std::unordered_map<K, std::vector<uint8_t>>& positive,
  std::unordered_map<K, std::vector<uint8_t>>& negative,
  const std::unordered_set<K>* positive_supset,
  const std::unordered_set<K>* negative_supset) const {
  IBLT<K> peeled = *this;
  while (true) {
    bool any = false;
    for (size_t i = 0; i < peeled.hashTable.size(); i++) {
      IBLTHashTableEntry<K>& entry = peeled.hashTable.at(i);
      if (!entry.isPure()) {
        continue;
      }
      if (entry.count == 1) {
        positive[entry.keySum] = entry.valueSum;
      } else if (entry.count == -1) {
        positive[entry.keySum] = entry.valueSum;
      } else {
        exit(-1);
      }
      peeled.insertImpl(-entry.count, entry.keySum, entry.valueSum);
      any = true;
    }
    if (!any) break;
  }
  return true;
}

template <typename K>
IBLT<K> IBLT<K>::operator-(const IBLT<K>& other) const {
  // IBLT's must be same params/size:
  assert(valueSize == other.valueSize);
  assert(hashTable.size() == other.hashTable.size());
  IBLT<K> result(*this);
  for (size_t i = 0; i < hashTable.size(); i++) {
    IBLTHashTableEntry<K>& e1 = result.hashTable.at(i);
    const IBLTHashTableEntry<K>& e2 = other.hashTable.at(i);
    e1.count -= e2.count;
    e1.keySum ^= e2.keySum;
    e1.keyCheck ^= e2.keyCheck;
    if (e1.empty()) {
      e1.valueSum.clear();
    } else {
      e1.addValue(e2.valueSum);
    }
  }
  return result;
}