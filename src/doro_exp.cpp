#include "doro.hpp"
#include "decoder.hpp"

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
  DoroCode doro(10000, 8, false, rng);
  auto rand_vec = random_nonrepetitive<int, mt19937>(10000000, 1000000, rng);
  std::shuffle(rand_vec.begin(), rand_vec.end(), rng);
  unordered_map<int, int> ground_truth, result;
  for (int i : rand_vec | views::take(1000))
  {
    ground_truth[i] = 1;
  }
  doro.encode(std::move(ground_truth));
  using DecoderType = DoroDecoder<int32_t, UpdatePQBackend::Set>;
  DoroDecoder<int32_t, UpdatePQBackend::Set> decoder;
  DecodeConfig config = {/*tk*/ 5, /*max_round*/ 100, /*ta*/ 10000, /*verbose*/ true, /*debug*/ false, /*lb*/ 0, /*ub*/ 1,
                                              PursuitChoice::L2};
  StopWatch sw; 
  // decoder.stage_decode(&doro, rand_vec, &config, result);
  // cout << format("Time: {}\n", sw.peek());
  // cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  // doro.reset();
  sw.reset_and_start();
  decoder.decode(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  doro.reset();
  sw.reset_and_start();
  decoder.decode2(&doro, rand_vec, &config, result);
  cout << format("Time: {}\n", sw.peek());
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  return 0;
}