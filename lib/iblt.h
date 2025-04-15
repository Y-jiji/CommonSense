#ifndef IBLT_H
#define IBLT_H

#include <inttypes.h>
#include <set>
#include <vector>
#include <fstream>
#include <unordered_map>
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
    int32_t count;
    uint64_t keySum;
    uint32_t keyCheck;
    std::vector<uint8_t> valueSum;

    bool isPure() const;
    bool empty() const;
    void addValue(const std::vector<uint8_t> v);
};

class IBLT
{
public:

    IBLT(size_t _expectedNumEntries, size_t _ValueSize, float hedge = 1.36, size_t _numHashes = 4);
    IBLT(const IBLT& other);
    virtual ~IBLT();

    void insert(uint64_t k, const std::vector<uint8_t> v);
    void erase(uint64_t k, const std::vector<uint8_t> v);

    bool get(uint64_t k, std::vector<uint8_t>& result) const;

    bool listEntries(std::set<std::pair<uint64_t,std::vector<uint8_t> > >& positive,
        std::set<std::pair<uint64_t,std::vector<uint8_t> > >& negative) const;

    // Subtract two IBLTs
    IBLT operator-(const IBLT& other) const;

    // For debugging:
    std::string DumpTable() const;

    std::string DumpEntry(size_t i) const;

    int hashTableSize();

    // these need to be public for python integration 
    size_t valueSize;
    size_t numHashes;  

private:
    void _insert(int plusOrMinus, uint64_t k, const std::vector<uint8_t> v);
    
    static std::string parameter_file; 

    std::vector<IBLTHashTableEntry> hashTable;
};

#endif /* IBLT_H */
