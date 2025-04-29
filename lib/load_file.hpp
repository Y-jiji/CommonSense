#pragma once
#include "libONIAK/oniakMath/overylarge.h"
#include "libONIAK/oniakRandom/orand.h"
#include "nlohmann/json.hpp"
#include <array>
#include <bitset>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <unordered_set>
#include <random>
#include <doro/key_types.hpp>

using namespace std;

std::bitset<256> stokey(const std::string &str) {
  std::array<unsigned char, 32> key{0};
  // skip '0x'
  for (size_t i = 2; i < 34; i += 2) {
    unsigned char a = str[i];
    unsigned char b = str[i];
    if (a >= '0' && a <= '9') {
      a -= '0';
    } else {
      a -= 'a';
    }
    if (b >= '0' && b <= '9') {
      b -= '0';
    } else {
      b -= 'a';
    }
    key[i - 2] = (a << 4) | b;
  }
  std::bitset<256> key_bitset;
  for (size_t i = 0; i < 32; ++i) {
    for (size_t j = 0; j < 8; ++j) {
      key_bitset.set(i * 8 + j, (key[i] >> (7 - j)) & 1);
    }
  }
  return key_bitset;
}

std::vector<std::bitset<256>> load(const std::string &path) {
  std::ifstream file(path);
  std::string line;
  std::vector<std::bitset<256>> values;
  while (std::getline(file, line)) {
    auto value = stokey(line);
    values.push_back(value);
  }
  return values;
}

template<typename K>
auto load_dataset_k256_or_k32(const nlohmann::json& config) {
  if (config.contains("A path") && config.contains("B path") &&
      config.contains("Intersect path")) {
    if constexpr(std::is_same<K, std::bitset<256>>()) {
      auto vecA = load(config.at("A path"));
      auto setA = std::unordered_set<K>(vecA.begin(), vecA.end());
      auto vecB = load(config.at("B path"));
      auto setB = std::unordered_set<K>(vecB.begin(), vecB.end());
      auto vecIntersect = load(config.at("Intersect path"));
      auto setIntersect =
          std::unordered_set<K>(vecIntersect.begin(), vecIntersect.end());
      auto A_size = setA.size();
      auto B_size = setB.size();
      auto A_minus_B_size = A_size - setIntersect.size();
      auto B_minus_A_size = B_size - setIntersect.size();
      auto setA_minus_B = setA;
      auto setB_minus_A = setB;
      for (const auto &v : setIntersect) {
        setA_minus_B.erase(v);
        setB_minus_A.erase(v);
      }
      return std::make_tuple(setA, setB, setA_minus_B, setB_minus_A,
                            A_minus_B_size, B_minus_A_size, setIntersect.size());
    } else {
      std::cout << "Bad data type, loading from file only supports u256" << std::endl;
      exit(-1);
    }
  } else {
    size_t seed = 0;
    if (config.contains("seed"))
      seed = config.at("seed");
    int A_size = config.at("A size");
    int B_size = config.at("B size");
    int A_minus_B_size = A_size - B_size;
    if (config.contains("A minus B size"))
      A_minus_B_size = config.at("A minus B size");
    int A_minus_B_size_minimum = std::max(0, A_size - B_size);
    if (A_minus_B_size < A_minus_B_size_minimum) {
      cout << "Warning: A_minus_B_size is too small. Reset to minimum "
              "possible value."
           << endl;
      A_minus_B_size = A_minus_B_size_minimum;
    }
    if (A_minus_B_size > A_size) {
      cout << "Warning: A_minus_B_size is larger than A_size. Reset to A_size."
           << endl;
      A_minus_B_size = A_size;
    }
    size_t A_union_B_size = B_size + A_minus_B_size;
    size_t A_intersect_B_size = A_size - A_minus_B_size;
    size_t B_minus_A_size = B_size - A_intersect_B_size;

    mt19937 rng(seed);
    auto rand_vec = ONIAK::random_nonrepetitive<K>(A_union_B_size, rng);

    // Extract A, B, A - B, B - A
    unordered_set<K> setA(rand_vec.begin(), rand_vec.begin() + A_size);
    unordered_set<K> setB(rand_vec.begin() + A_minus_B_size,
                                  rand_vec.end());
    unordered_set<K> setA_minus_B(rand_vec.begin(),
                                          rand_vec.begin() + A_minus_B_size);
    unordered_set<K> setB_minus_A(rand_vec.begin() + A_size,
                                          rand_vec.end());
    return std::make_tuple(setA, setB, setA_minus_B, setB_minus_A,
                           A_minus_B_size, B_minus_A_size, A_intersect_B_size);
  }
}