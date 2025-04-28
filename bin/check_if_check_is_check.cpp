// #include <doro/load_file.hpp>
// #include <doro/bch_wrapper.hpp>
// #include <doro/decoder.hpp>
// #include <doro/doro.hpp>
// #include <doro/probability.hpp>
// #include <doro/rans_wrapper.hpp>
// #include <doro/sha256.hpp>

// #include "libONIAK/oniakDataStructure/ohist.h"
// #include "libONIAK/oniakDebug/odebug.h"
// #include "libONIAK/oniakMath/overylarge.h"
// #include "libONIAK/oniakRandom/orand.h"
// #include "libONIAK/oniakTimer/otime.h"
// #include "nlohmann/json.hpp"

// #include <doro/iblt.h>
// #include <algorithm>
// #include <array>
// #include <bitset>
// #include <boost/math/distributions/binomial.hpp>
// #include <cmath>
// #include <filesystem>
// #include <format>
// #include <fstream>
// #include <iostream>
// #include <print>
// #include <random>
// #include <ranges>
// #include <stdint.h>
// #include <string>
// #include <tuple>
// #include <unordered_map>
// #include <unordered_set>
// #include <set>
// #include <vector>

// using namespace std;
// using namespace Doro;
// using namespace ONIAK;
// using namespace nlohmann;
// using namespace boost::math;

// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         println("Error: please specify config json file.");
//         exit(1);
//     }
//     // read configuration file from input
//     std::ifstream fin(argv[1]);
//     json config = json::parse(fin);
//     // load dataset
//     auto [setA, setB, setA_minus_B, setB_minus_A, A_minus_B_size, B_minus_A_size,
//         A_intersect_B_size] = load_dataset_k256(config);
//     // see if hash for key check collide
//     auto setAKeySum = std::set<uint64_t>();
//     for (auto x: setA) {
//         auto checkX = HashWithSalt(Doro::to_binary_vector<VeryLargeInt<256>>(x), N_HASHCHECK);
//         if (setAKeySum.contains(checkX)) {
//             std::cout << "COLLISION DETECTED " << std::endl;
//             std::cout << "KEY: " << x << std::endl;
//             std::cout << "HASHED: " << checkX << std::endl;
//             exit(-1);
//         }
//         setAKeySum.insert(checkX);
//     }
//     // see if hash for key check collide
//     auto setBKeySum = std::set<uint64_t>();
//     for (auto x: setB) {
//         auto checkX = HashWithSalt(to_binary_vector(x), N_HASHCHECK);
//         if (setBKeySum.contains(checkX)) {
//             std::cout << "COLLISION DETECTED " << std::endl;
//             std::cout << "KEY: " << x << std::endl;
//             std::cout << "HASHED: " << checkX << std::endl;
//             exit(-1);
//         }
//         setBKeySum.insert(checkX);
//     }
// }