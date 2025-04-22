#include <cassert>
#include <iostream>
#include <list>
#include <sstream>
#include <utility>
#include <iomanip>      // std::setbase
#include <doro/iblt.h>
#include <doro/key_types.hpp>
#include <unordered_set>
#include <doro/sha256.hpp>
#include <cstring>
#define DEBUGPRINT

template <typename IndexType>
std::vector<uint8_t> to_binary_vector(const IndexType& x) {
    constexpr size_t size = sizeof(IndexType);
    std::vector<uint8_t> result(size);
    std::memcpy(result.data(), &x, size);
    return result;
}

bool IBLTHashTableEntry::isPure() const {
    if (count == 1 || count == -1) {
        uint64_t check = HashWithSalt(to_binary_vector(keySum), N_HASHCHECK);
        return (keyCheck == check);
    }
    return false;
}

bool IBLTHashTableEntry::empty() const {
    return (count == 0ull && keySum == IndexType() && keyCheck == 0ull);
}

void IBLTHashTableEntry::addValue(const std::vector<uint8_t> v) {
    valueSum.resize(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        valueSum[i] ^= v[i];
    }
}

int IBLT::hashTableSize() {
    return(hashTable.size());
}

IBLT::IBLT(size_t _expectedNumEntries, size_t _valueSize, float hedge, size_t _numHashes) :
    valueSize(_valueSize),
    numHashes(_numHashes)
{
    assert(numHashes > 1);
    size_t nEntries = (int)(((float)_expectedNumEntries) * hedge);
    nEntries = (nEntries + numHashes - 1) / numHashes * numHashes;
    hashTable.resize(nEntries);
}

IBLT::IBLT(const IBLT& other) {
    valueSize = other.valueSize;
    hashTable = other.hashTable;
    numHashes = other.numHashes;
}

IBLT::~IBLT() {}

void IBLT::_insert(int plusOrMinus, const IndexType key, const std::vector<uint8_t> v) {
    assert(v.size() == valueSize);
    assert(key != IndexType());
    assert(plusOrMinus == 1 || plusOrMinus == -1);
    std::unordered_set<uint64_t> hashPositions;
    for (size_t i = 0; i < numHashes; i++) {
        uint64_t h = HashWithSalt(to_binary_vector(key), i) % hashTable.size();
        if (hashPositions.contains(h)) continue;
        hashPositions.insert(h);
        IBLTHashTableEntry& entry = hashTable.at(h);
        // DEBUGPRINT std::cout << "POSITION: " << h << std::endl;
        // DEBUGPRINT std::cout << "SIGN: " << plusOrMinus << std::endl;
        // DEBUGPRINT std::cout << "INSERT KEY: " << key << std::endl;
        // DEBUGPRINT std::cout << "INSERT KEY CHECK: " << HashWithSalt(to_binary_vector(key), N_HASHCHECK) << std::endl;
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

void IBLT::insert(const IndexType& k, const std::vector<uint8_t> v) {
    _insert(1, k, v);
}

void IBLT::erase(const IndexType& k, const std::vector<uint8_t> v) {
    _insert(-1, k, v);
}

bool IBLT::listEntries(
    std::set<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>& positive,
    std::set<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>& negative,
    const std::unordered_set<IndexType>* positive_supset,
    const std::unordered_set<IndexType>* negative_supset
) const {
    IBLT peeled = *this;
    while (true) {
        bool any = false;
        for (size_t i = 0; i < peeled.hashTable.size(); i++) {
            IBLTHashTableEntry& entry = peeled.hashTable.at(i);
            if (!entry.isPure()) { continue; }
            // DEBUGPRINT std::cout << "POSIITON " << i << " IS PURE" << std::endl;
            auto keySumVec = to_binary_vector(entry.keySum);
            if (entry.count == 1) {
                positive.insert(std::make_pair(keySumVec, entry.valueSum));
                // DEBUGPRINT if (positive_supset && !positive_supset->contains(entry.keySum)) {
                    // DEBUGPRINT std::cout << "DECODED KEY NOT IN VALID RANGE (POSITIVE)" << std::endl;
                    // DEBUGPRINT std::cout << entry.keySum << std::endl;
                    // DEBUGPRINT std::cout << "CHECK " << entry.keyCheck << std::endl;
                    // DEBUGPRINT if (negative_supset->contains(entry.keySum)) {
                        // DEBUGPRINT std::cout << "IT SHOULD BE IN NEGATIVE" << std::endl;
                    // DEBUGPRINT }
                    // DEBUGPRINT if (negative.contains(std::make_pair(keySumVec, entry.valueSum))) {
                        // DEBUGPRINT std::cout << "IT IS IN NEGATIVE" << std::endl;
                    // DEBUGPRINT }
                    // DEBUGPRINT exit(-1);
                // DEBUGPRINT }
            } else if (entry.count == -1) {
                negative.insert(std::make_pair(keySumVec, entry.valueSum));
                // DEBUGPRINT if (negative_supset && !negative_supset->contains(entry.keySum)) {
                    // DEBUGPRINT std::cout << "DECODED KEY NOT IN VALID RANGE (NEGATIVE)" << std::endl;
                    // DEBUGPRINT std::cout << entry.keySum << std::endl;
                    // DEBUGPRINT std::cout << "CHECK " << entry.keyCheck << std::endl;
                    // DEBUGPRINT if (positive_supset->contains(entry.keySum)) {
                    // DEBUGPRINT   std::cout << "IT SHOULD BE IN POSITIVE" << std::endl;
                    // DEBUGPRINT }
                    // DEBUGPRINT if (positive.contains(std::make_pair(keySumVec, entry.valueSum))) {
                    // DEBUGPRINT   std::cout << "IT IS IN POSITIVE" << std::endl;
                    // DEBUGPRINT }
                    // DEBUGPRINT exit(-1);
                // DEBUGPRINT }
            } else {
                // DEBUGPRINT std::cout << "WHAT?" << std::endl;
                exit(-1);
            }
            peeled._insert(-entry.count, entry.keySum, entry.valueSum);
            any = true;
        }
        if (!any) break;
    }
    // for (size_t i = 0; i < peeled.hashTable.size() / numHashes; i++) {
    //     if (peeled.hashTable.at(i).count != 0) return false;
    // }
    return true;
}

IBLT IBLT::operator-(const IBLT& other) const {
    // IBLT's must be same params/size:
    assert(valueSize == other.valueSize);
    assert(hashTable.size() == other.hashTable.size());

    IBLT result(*this);
    for (size_t i = 0; i < hashTable.size(); i++) {
        IBLTHashTableEntry& e1 = result.hashTable.at(i);
        const IBLTHashTableEntry& e2 = other.hashTable.at(i);
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