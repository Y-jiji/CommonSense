#pragma once

#include <inttypes.h>
#include <set>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <optional>
#include "libONIAK/oniakMath/overylarge.h"
#include <doro/key_types.hpp>

//
// Invertible Bloom Lookup Table implementation
// References:
//
// "What's the Difference? Efficient Set Reconciliation
// without Prior Context" by Eppstein, Goodrich, Uyeda and
// Varghese
//
// "Invertible Bloom Lookup Tables" by Goodrich and
// Mitzenmacher
//

class IBLTHashTableEntry
{
public:
    int64_t  count = 0;
    IndexType keySum;
    uint64_t keyCheck = 0;
    std::vector<uint8_t> valueSum{};

    bool isPure() const;
    bool empty() const;
    void addValue(const std::vector<uint8_t> v);
};

class IBLT
{
public:

    IBLT(size_t _expectedNumEntries, size_t _ValueSize, float hedge = 4.0, size_t _numHashes = 4);
    IBLT(const IBLT& other);
    virtual ~IBLT();

    void insert(const IndexType& k, const std::vector<uint8_t> v);
    void erase(const IndexType& k, const std::vector<uint8_t> v);

    bool listEntries(
        std::set<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>& positive,
        std::set<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>& negative,
        const std::unordered_set<IndexType>*   positive_supset = nullptr,
        const std::unordered_set<IndexType>*   negative_supset = nullptr
    ) const;

    // Subtract two IBLTs
    IBLT operator-(const IBLT& other) const;

    int hashTableSize();

    // these need to be public for python integration 
    size_t valueSize;
    size_t numHashes;

private:
    void _insert(int plusOrMinus, const IndexType k, const std::vector<uint8_t> v);

    std::vector<IBLTHashTableEntry> hashTable;
};