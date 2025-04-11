#include <doro/bch_wrapper.hpp>
#include <doro/doro.hpp>
#include <doro/decoder.hpp>
#include <doro/probability.hpp>
#include <doro/rans_wrapper.hpp>

#include "libONIAK/oniakDataStructure/ohist.h"
#include "libONIAK/oniakDebug/odebug.h"
#include "libONIAK/oniakMath/overylarge.h"
#include "libONIAK/oniakRandom/orand.h"
#include "libONIAK/oniakTimer/otime.h"
#include "nlohmann/json.hpp"
#include "IBLT-opt/iblt.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <format>
#include <iostream>
#include <print>
#include <random>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <boost/math/distributions/binomial.hpp>

using namespace std;
using namespace Doro;
using namespace ONIAK;
using namespace nlohmann;
using namespace boost::math;

enum class Party { Alis = 0, Bela = 1 };
enum class Status { CollisionAvoiding = 0, CollisionResolving = 1, Finished = 2 };
using CounterType = int16_t;
using IndexType = VeryLargeInt<256>;
using DoroCodeType = DoroCode<IndexType, CounterType>;

// The following constants are used for parameter tuning
constexpr double diff_coding_error_min = 1e-9;
constexpr double diff_coding_error_max = 0.1;
constexpr double diff_coding_error_final = 1e-6;
constexpr double bch_block_error_rate = 0.01;
constexpr int iblt_element_margin = 100;
constexpr double iblt_element_factor = 1.25;
constexpr int iblt_cell_size = sizeof(IndexType) * 8 + 32; // 32 key-check bits.

int main(int argc, char* argv[]) {
  if (argc < 2) {
    println("Error: please specify config json file.");
    exit(1);
  }
  // configure experiement configuration input
  std::ifstream fin(argv[1]);
  json config = json::parse(fin);

  // configure experiment output
  std::string save_path = config["result filename"];
  std::ofstream fout(save_path);
  if (!fout.is_open()) {
    cout << "Error: cannot open file " << save_path << endl;
    return 1;
  }

  size_t seed = 0;
  if (config.contains("seed")) seed = config.at("seed");

  // a helper function to convert value from bitvector to bitset
  auto to_bit_set = [](const std::vector<uint8_t>& bv) {
    IndexType result;
    std::memcpy(&result, bv.data(), bv.size());
    return result;
  };

  // Generate a vector
  // [  A (] B  )
  int A_size = config.at("A size");
  int B_size = config.at("B size");
  int A_minus_B_size = A_size - B_size;
  if (config.contains("A minus B size")) A_minus_B_size = config.at("A minus B size");
  int A_minus_B_size_minimum = std::max(0, A_size - B_size);
  if (A_minus_B_size < A_minus_B_size_minimum) {
    cout << "Warning: A_minus_B_size is too small. Reset to minimum possible value." << endl;
    A_minus_B_size = A_minus_B_size_minimum;
  }
  if (A_minus_B_size > A_size) {
    cout << "Warning: A_minus_B_size is larger than A_size. Reset to A_size." << endl;
    A_minus_B_size = A_size;
  }
  int A_union_B_size = B_size + A_minus_B_size;
  int A_intersect_B_size = A_size - A_minus_B_size;
  int B_minus_A_size = B_size - A_intersect_B_size;

  mt19937 rng(seed);
  auto rand_vec = random_nonrepetitive<IndexType>(A_union_B_size, rng);

  // Extract A, B, A - B, B - A
  unordered_set<IndexType> setA(rand_vec.begin(), rand_vec.begin() + A_size);
  unordered_set<IndexType> setB(rand_vec.begin() + A_minus_B_size, rand_vec.end());
  unordered_set<IndexType> setA_minus_B(rand_vec.begin(), rand_vec.begin() + A_minus_B_size);
  unordered_set<IndexType> setB_minus_A(rand_vec.begin() + A_size, rand_vec.end());

  //                   //
  // Set Reconcilation //
  //                   //
  int round = 0;
  // compute IBLT for A and B
  IBLT setA_IBLT((B_minus_A_size + A_minus_B_size), sizeof(IndexType));
  for (auto a: setA) {
    setA_IBLT.insert(static_cast<uint64_t>(a), to_binary_vector(a));
  }
  IBLT setB_IBLT((B_minus_A_size + A_minus_B_size), sizeof(IndexType));
  for (auto b: setB) {
    setB_IBLT.insert(static_cast<uint64_t>(b), to_binary_vector(b));
  }
  // patch B with A's iblt
  auto setA_minus_B_IBLT = setA_IBLT - setB_IBLT;
  auto A_minus_B_IBLT_estimated = std::set<std::pair<uint64_t, std::vector<uint8_t>>>();
  auto B_minus_A_IBLT_estimated = std::set<std::pair<uint64_t, std::vector<uint8_t>>>();
  setA_minus_B_IBLT.listEntries(A_minus_B_IBLT_estimated, B_minus_A_IBLT_estimated);
  for (auto [key, value]: A_minus_B_IBLT_estimated) {
    setA.erase(to_bit_set(value));
  }
  for (auto [key, value]: B_minus_A_IBLT_estimated) {
    setB.erase(to_bit_set(value));
  }
  // add to transmitted size
  int transmitted_size_A_to_B = 
    setA_IBLT.hashTableSize() * iblt_cell_size;
  int transmitted_size_B_to_A = 
    accumulate(B_minus_A_IBLT_estimated.begin(), B_minus_A_IBLT_estimated.end(), 0, 
    [](size_t sum, const auto& vec) { 
      // reducing on iblt transmitted size
      return sum + sizeof(std::pair<uint64_t, std::vector<uint8_t>>) + vec.second.size(); 
    });

  //                 //
  // Validate Result //
  //                 //
  
  // a helper function for counting elements of x in y
  auto count_x_in_y = [](
    const std::unordered_set<IndexType>& x, 
    const std::unordered_set<IndexType>& y
  ) {
    size_t count = 0;
    for (auto i: x) { count += y.contains(i) ? 1 : 0; }
    return count;
  };
  // validate setA and setB are just the intersection
  bool success = true;
  if (
    count_x_in_y(setA, setA_minus_B) != 0 ||
    count_x_in_y(setB, setB_minus_A) != 0 ||
    count_x_in_y(setA, setB_minus_A) != 0 ||
    count_x_in_y(setB, setA_minus_B) != 0
  ) {
    success = false;
    std::cout << "A in A\\B: " << count_x_in_y(setA, setA_minus_B) << std::endl;
    std::cout << "A in B\\A: " << count_x_in_y(setA, setB_minus_A) << std::endl;
    std::cout << "B in A\\B: " << count_x_in_y(setB, setA_minus_B) << std::endl;
    std::cout << "B in B\\A: " << count_x_in_y(setB, setB_minus_A) << std::endl;
    println("failure");
  }
  cout << "iblt total cost a to b " << transmitted_size_A_to_B << std::endl;
  cout << "iblt total cost b to a " << transmitted_size_B_to_A << std::endl;
  cout << "round " << round << std::endl;
  config["iblt total cost a to b"] = transmitted_size_A_to_B;
  config["iblt total cost b to a"] = transmitted_size_B_to_A;
  config["success"] = success;
  config["round"] = round;
  fout << std::setw(4) << config << endl;  // indent 4
}