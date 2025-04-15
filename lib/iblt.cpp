#include <cassert>
#include <iostream>
#include <list>
#include <sstream>
#include <utility>
#include <iomanip>      // std::setbase
#include <doro/iblt_hash.hpp>
#include <doro/iblt.h>
#include <doro/utilstrencodings.h>

static const size_t N_HASHCHECK = 11;

bool IBLTHashTableEntry::isPure() const {
    if (count == 1 || count == -1) {
        uint64_t check = MurmurHashTwice(N_HASHCHECK, ToVec(keySum));
        return (keyCheck == check);
    }
    return false;
}

bool IBLTHashTableEntry::empty() const {
    return (count == 0 && keySum == 0 && keyCheck == 0);
}

void IBLTHashTableEntry::addValue(const std::vector<uint8_t> v) {
    if (v.empty()) {
        return;
    }
    if (valueSum.size() < v.size()) {
        valueSum.resize(v.size());
    }
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
    size_t nEntries = (int)(((float)_expectedNumEntries) * (hedge));
    while (numHashes * (nEntries / numHashes) != nEntries) ++nEntries;
    hashTable.resize(nEntries);
}

IBLT::IBLT(const IBLT& other) {
    valueSize = other.valueSize;
    hashTable = other.hashTable;
    numHashes = other.numHashes;
}

IBLT::~IBLT() {}

void IBLT::_insert(int plusOrMinus, uint64_t k, const std::vector<uint8_t> v) {
    assert(v.size() == valueSize);

    std::vector<uint8_t> kvec = ToVec(k);

    size_t bucketsPerHash = hashTable.size() / numHashes;
    for (size_t i = 0; i < numHashes; i++) {
        size_t startEntry = i * bucketsPerHash;
        uint64_t h = MurmurHashTwice(i, kvec);
        IBLTHashTableEntry& entry = hashTable.at(startEntry + (h % bucketsPerHash));
        entry.count += plusOrMinus;
        entry.keySum ^= k;
        entry.keyCheck ^= MurmurHashTwice(N_HASHCHECK, kvec);
        if (entry.empty()) {
            entry.valueSum.clear();
        } else {
            entry.addValue(v);
        }
    }
}

void IBLT::insert(uint64_t k, const std::vector<uint8_t> v) {
    _insert(1, k, v);
}

void IBLT::erase(uint64_t k, const std::vector<uint8_t> v) {
    _insert(-1, k, v);
}

bool IBLT::get(uint64_t k, std::vector<uint8_t>& result) const {
    result.clear();

    std::vector<uint8_t> kvec = ToVec(k);

    size_t bucketsPerHash = hashTable.size() / numHashes;
    for (size_t i = 0; i < numHashes; i++) {
        size_t startEntry = i * bucketsPerHash;
        uint64_t h = MurmurHashTwice(i, kvec);
        std::cout << "murmur hash " << h << std::endl;
        const IBLTHashTableEntry& entry = hashTable.at(startEntry + (h % bucketsPerHash));
        if (entry.empty()) {
            // Definitely not in table. Leave
            // result empty, return true.
            return true;
        } else if (entry.isPure()) {
            if (entry.keySum == k) {
                // Found!
                result.assign(entry.valueSum.begin(), entry.valueSum.end());
                return true;
            } else {
                // Definitely not in table.
                return true;
            }
        }
    }

    // Don't know if k is in table or not; "peel" the IBLT to try to find it:
    IBLT peeled = *this;
    size_t nErased = 0;
    for (size_t i = 0; i < peeled.hashTable.size(); i++) {
        IBLTHashTableEntry& entry = peeled.hashTable.at(i);
        if (entry.isPure()) {
            if (entry.keySum == k) {
                // Found!
                result.assign(entry.valueSum.begin(), entry.valueSum.end());
                return true;
            }
            ++nErased;
            // std::cout<<" peeling\n";
            peeled._insert(-entry.count, entry.keySum, entry.valueSum);
        }
    }
    if (nErased > 0) {
        // Recurse with smaller IBLT
        return peeled.get(k, result);
    }

    return false;
}

bool IBLT::listEntries(std::set<std::pair<uint64_t, std::vector<uint8_t> > >& positive,
    std::set<std::pair<uint64_t, std::vector<uint8_t> > >& negative) const {
    IBLT peeled = *this;
    size_t nErased = 0;
    do {
        nErased = 0;
        for (size_t i = 0; i < peeled.hashTable.size(); i++) {
            IBLTHashTableEntry& entry = peeled.hashTable.at(i);
            if (entry.isPure()) {
                if (entry.count == 1) {
                    positive.insert(std::make_pair(entry.keySum, entry.valueSum));
                } else {
                    negative.insert(std::make_pair(entry.keySum, entry.valueSum));
                }
                peeled._insert(-entry.count, entry.keySum, entry.valueSum);
                ++nErased;
            }
        }
    } while (nErased > 0);

    // If any buckets for one of the hash functions is not empty,
    // then we didn't peel them all:
    for (size_t i = 0; i < peeled.hashTable.size() / numHashes; i++) {
        if (peeled.hashTable.at(i).empty() != true) return false;
    }
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

// For debugging during development:
std::string IBLT::DumpEntry(size_t i) const {
    std::ostringstream result;
    const IBLTHashTableEntry& entry = hashTable.at(i);
    result << std::dec << entry.count << "\t" << entry.keySum << "\t";
    result << ((MurmurHashTwice(N_HASHCHECK, ToVec(entry.keySum)) == entry.keyCheck) || (entry.count == 0) ? "true" : "false");
    result << "\t" << std::setw(8) << std::setbase(16) << unsigned(entry.keyCheck) << "\t[";
    for (uint8_t i : entry.valueSum)
        result << std::setw(1) << std::setbase(16) << unsigned(i);
    result << "]\n";
    return result.str();
}

// For debugging during development:
std::string IBLT::DumpTable() const {
    std::ostringstream result;
    for (size_t i = 0; i < hashTable.size(); i++) {
        const IBLTHashTableEntry& entry = hashTable.at(i);
        if (entry.count != 0 || entry.keySum != 0 || ((MurmurHashTwice(N_HASHCHECK, ToVec(entry.keySum)) == entry.keyCheck))) {
            result << std::dec << i << "\t";
            result << std::dec << entry.count << "\t" << entry.keySum << "\t";
            result << ((MurmurHashTwice(N_HASHCHECK, ToVec(entry.keySum)) == entry.keyCheck) || (entry.count == 0) ? "true" : "false");
            result << "\t" << std::setw(8) << std::setbase(16) << unsigned(entry.keyCheck) << "\t[";
            for (uint8_t i : entry.valueSum)
                result << std::setw(1) << std::setbase(16) << unsigned(i);
            result << "]\n";
        }
    }
    return result.str();
}