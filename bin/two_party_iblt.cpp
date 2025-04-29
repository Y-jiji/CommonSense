#include <doro/load_file.hpp>
#include <doro/bch_wrapper.hpp>
#include <doro/doro_decoder.hpp>
#include <doro/doro.hpp>
#include <doro/probability.hpp>
#include <doro/rans_wrapper.hpp>

#include "doro/key_types.hpp"
#include "libONIAK/oniakDataStructure/ohist.h"
#include "libONIAK/oniakDebug/odebug.h"
#include "libONIAK/oniakMath/overylarge.h"
#include "libONIAK/oniakRandom/orand.h"
#include "libONIAK/oniakTimer/otime.h"
#include "nlohmann/json.hpp"

#include <doro/iblt.h>
#include <boost/math/distributions/binomial.hpp>
#include <cmath>
#include <iostream>
#include <print>
#include <stdint.h>

using namespace std;
using namespace Doro;
using namespace ONIAK;
using namespace nlohmann;
using namespace boost::math;

enum class Party { Alis = 0, Bela = 1 };
enum class Status { CollisionAvoiding = 0, CollisionResolving = 1, Finished = 2 };
using CounterType = int16_t;
using IndexType = VeryLargeInt<256>;

// The following constants are used for parameter tuning
constexpr double diff_coding_error_min = 1e-9;
constexpr double diff_coding_error_max = 0.1;
constexpr double diff_coding_error_final = 1e-6;
constexpr double bch_block_error_rate = 0.01;
constexpr int iblt_element_margin = 100;
constexpr double iblt_element_factor = 1.25;
constexpr int iblt_cell_size = 32 * 8 + 32; // 32 key-check bits.

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

  // a helper function to convert value from bitvector to bitset
  auto to_bit_set = [](const std::vector<uint8_t>& bv) {
    IndexType result((std::bitset<256>()));
    std::memcpy(&result, bv.data(), bv.size());
    return result;
  };

  // Generate a vector
  // [  A (] B  )
  auto [setA, setB, setA_minus_B, setB_minus_A, A_minus_B_size, B_minus_A_size,
        A_intersect_B_size] = load_dataset_k256_or_k32<IndexType>(config);
  std::cout << "A minus B: " << A_minus_B_size << std::endl;
  std::cout << "B minus A: " << B_minus_A_size << std::endl;

  //                   //
  // Set Reconcilation //
  //                   //
  int round = 0;
  int transmitted_size_A_to_B = 0;
  int transmitted_size_B_to_A = 0;

  while (true) {
    // stop flag
    bool stop = true;
    // compute IBLT for A and B
    IBLT<IndexType> setA_IBLT((B_minus_A_size + A_minus_B_size), sizeof(IndexType));
    std::cout << "BUILD IBLT FOR A" << std::endl;
    for (auto a: setA) {
      setA_IBLT.insert(a, to_binary_vector(a));
    }
    std::cout << "BUILD IBLT FOR B" << std::endl;
    IBLT<IndexType> setB_IBLT((B_minus_A_size + A_minus_B_size),
                              sizeof(IndexType));
    for (auto b: setB) {
      setB_IBLT.insert(b, to_binary_vector(b));
    }
    // patch B with A's iblt
    auto setA_minus_B_IBLT = setA_IBLT - setB_IBLT;
    std::cout << "BUILD DIFFERENCE" << std::endl;
    auto A_minus_B_IBLT_estimated =
        std::unordered_map<VeryLargeInt<256>, std::vector<uint8_t>>();
    auto B_minus_A_IBLT_estimated =
        std::unordered_map<VeryLargeInt<256>, std::vector<uint8_t>>();
    std::cout << "LIST ENTRIES" << std::endl;
    auto ok = setA_minus_B_IBLT.listEntries(A_minus_B_IBLT_estimated, B_minus_A_IBLT_estimated, &setA_minus_B, &setB_minus_A);
    if (!ok) {
      std::cout << "DECODING FAILURE!" << std::endl;
      exit(-1);
    }
    std::cout << "REMOVING DELTA" << std::endl;
    std::cout << "A minus B decoded: " << A_minus_B_IBLT_estimated.size() << std::endl;
    for (auto [key, value]: A_minus_B_IBLT_estimated) {
      if (!setA.contains(to_bit_set(value))) {
        std::cout << "KEY in decoded A minus B not in setA" << std::endl;
        if (setB.contains(to_bit_set(value))) {
          std::cout << "It is in B" << std::endl;
        }
      }
      setA.erase(to_bit_set(value));
      stop = false;
    }
    std::cout << "REST: " << std::endl;
    for (auto x: setA) {
      if (!setA_minus_B.contains(x)) continue;
      std::cout << x << std::endl;
    }
    std::cout << "B minus A decoded: " << B_minus_A_IBLT_estimated.size() << std::endl;
    for (auto [key, value]: B_minus_A_IBLT_estimated) {
      if (!setB.contains(to_bit_set(value))) {
        std::cout << "KEY in B minus A not in setB" << std::endl;
        if (setA.contains(to_bit_set(value))) {
          std::cout << "It is in A" << std::endl;
        }
      }
      setB.erase(to_bit_set(value));
      stop = false;
    }
    std::cout << "REST: " << std::endl;
    for (auto x: setB) {
      if (!setB_minus_A.contains(x)) continue;
      std::cout << x << std::endl;
    }
    if (stop) {
      std::cout << "STOP" << std::endl;
      break;
    }
    // add to transmitted size
    transmitted_size_A_to_B += 
      setA_IBLT.hashTableSize() * iblt_cell_size;
    transmitted_size_B_to_A += 
      accumulate(B_minus_A_IBLT_estimated.begin(), B_minus_A_IBLT_estimated.end(), 0, 
      [](size_t sum, const auto& vec) { 
        // reducing on iblt transmitted size
        return sum + 256 + vec.second.size(); 
      });
    std::cout << "ROUND " << round << std::endl;
    std::cout << "A -> B: " << transmitted_size_A_to_B << std::endl;
    std::cout << "B -> A: " << transmitted_size_B_to_A << std::endl;
    round++;
  }

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
  return 0;
}
