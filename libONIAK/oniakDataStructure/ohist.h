#pragma once

#include <unordered_map>
#include <vector>

namespace ONIAK {

template <typename T>
std::unordered_map<T, int> frequency_count(const std::vector<T>& vec) {
    std::unordered_map<T, int> result;
    for (const T& item : vec) {
        auto iter = result.find(item);
        int count = (iter != result.end()) ? iter->second : 0;
        ++count;
        result[item] = count;
    }
    return result;
}

// hashable tuple to be used as keys for hash containers
// hash values are simply XOR's of individual hashes

struct variadic_hash {
template <typename T>
size_t operator()(const T& t) {
    return std::hash<T>()(t);
}

template <typename T, typename ... Targs>
size_t operator()(const T& t, const Targs&... fargs) {
    size_t result = std::hash<T>()(t);
    return result ^ variadic_hash()(fargs...);
}
};

template <typename ... Targs>
struct TupleHash {
    std::size_t operator() (const std::tuple<Targs...>& tuple) const {
        if constexpr (sizeof...(Targs) == 0) return 0;
        return std::apply(
            variadic_hash(),
            tuple);
    }
};

}


