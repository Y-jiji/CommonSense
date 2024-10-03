#include "doro.hpp"
#include "decoder.hpp"
#include "probability.hpp"

#include "oniakRandom/orand.h"
#include "oniakTimer/otime.h"

#include <format>
#include <iostream>
#include <random>
#include <ranges>
#include <set>
#include <unordered_map>

using namespace std;
using namespace Doro;
using namespace ONIAK;




int main()
{
  mt19937 rng(0);
  DoroCode doro(70000, 8, false, rng);
  auto rand_vec = random_nonrepetitive<int, mt19937>(10000000, 1000000, rng);
  std::shuffle(rand_vec.begin(), rand_vec.end(), rng);
  unordered_map<int, int> ground_truth, *result;
  for (int i : rand_vec | views::take(10000))
  {
    ground_truth[i] = 1;
  }
  doro.encode(std::move(ground_truth));
  DoroDecoder<int32_t, UpdatePQBackend::PriorityQueue> decoder1;
  DoroDecoder<int32_t, UpdatePQBackend::Set> decoder2;
  DoroDecoder<int32_t, UpdatePQBackend::UpdatablePriorityQueue> decoder3;
  DecodeConfig config = {/*tk*/ 5, /*max_round*/ 100, /*ta*/ 10000, /*verbose*/ false, /*debug*/ true, /*lb*/ 0, /*ub*/ 1,
                                              PursuitChoice::L2};
  StopWatch sw; 
  // decoder.stage_decode(&doro, rand_vec, &config, result);
  // cout << format("Time: {}\n", sw.peek());
  // cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  // doro.reset();
  sw.reset_and_start();
  decoder1.decode(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  doro.reset();
  sw.reset_and_start();
  decoder1.decode_resense(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  doro.reset();
  sw.reset_and_start();
  decoder2.decode(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  doro.reset();
  sw.reset_and_start();
  decoder2.decode_resense(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  doro.reset();
  sw.reset_and_start();
  decoder3.decode(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  doro.reset();
  sw.reset_and_start();
  decoder3.decode_resense(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  return 0;
}

/* running result 
A = 1M, diff = 1000, k = 8, d = 10000
Time: 10.011230582    <- PriorityQueue without resense is best
Nonzero: 0
 MAE: 0
Time: 11.636769064
Nonzero: 0


 MAE: 0
Time: 13.661667001
Nonzero: 0
 MAE: 0
Time: 15.37705184
Nonzero: 0
 MAE: 0
Time: 11.050927596
Nonzero: 0
 MAE: 0
Time: 12.832666099
Nonzero: 0
 MAE: 0

 A = 1M, diff = 10000, k = 8, d = 70000
 Time: 13.17794145
Nonzero: 0
 MAE: 0
Time: 15.090222476
Nonzero: 0
 MAE: 0
Time: 18.3383255
Nonzero: 0
 MAE: 0
Time: 20.641058072
Nonzero: 0
 MAE: 0
Time: 14.561330016
Nonzero: 0
 MAE: 0
Time: 16.726342531
Nonzero: 0
 MAE: 0
 */