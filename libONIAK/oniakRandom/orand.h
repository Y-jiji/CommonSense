#pragma once

#include "../global.h"
#include "../oniakMath/overylarge.h"

#include <bit>
#include <concepts>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ONIAK
{ 
    // specialized overload for consistency test with Python-based Graphene
    uint64_t random_value(std::mt19937 &rand) {
        uint64_t res = rand();
        return (res << 32) | rand();
    }

    template <std::integral T, typename RandomGenerator>
    T random_value(RandomGenerator &rand) {
        std::uniform_int_distribution<T> dist;
        return dist(rand);
    }

    template <isVeryLargeInt T, typename RandomGenerator>
    T random_value(RandomGenerator &rand) {
        constexpr int rg_width = std::bit_width(RandomGenerator::max());
        T res;
        int rem_width = res.size;
        while (rem_width > 0) {
            int width = std::min(rem_width, rg_width);
            T rand_val(rand() & (-1ull >> (64 - width)));
            res <<= width;
            res |= rand_val;
            rem_width -= width;
        }
        return res;
    }

    // Generate n random numbers without replacement in [0, U)
    template <typename T, typename RandomGenerator>
    std::vector<T> random_shuffle(int U, int n, RandomGenerator &rand) {
        std::vector<T> pool(U);
        std::iota(pool.begin(), pool.end(), 0);
        std::shuffle(pool.begin(), pool.end(), rand);
        pool.resize(n);
        return pool;
    }

    template <typename T, typename RandomGenerator>
    std::vector<T> random_nonrepetitive(int n, RandomGenerator &rand) {
        std::vector<T> pool(n);
        std::unordered_set<T> set;
        for (T& val : pool) {
          do {
            val = random_value<T>(rand);
          } while (!set.insert(val).second);
        }
        return pool;
    }
}

