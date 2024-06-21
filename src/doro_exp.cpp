#include "doro.hpp"
#include "decoder.hpp"

#include "oniakRandom/orand.h"

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
  DoroCode doro(8000, 8, false, rng);
  auto rand_vec = random_nonrepetitive<int, mt19937>(1000000, 100000, rng);
  unordered_map<int, int> ground_truth, result;
  for (int i : rand_vec | views::take(1000))
  {
    ground_truth[i] = 1;
  }
  doro.encode(std::move(ground_truth));
  DoroDecoder<> decoder;
  DoroDecoder<>::DecodeConfig config = {/*tk*/ 100, /*max_round*/ 100, /*ta*/ 10, /*verbose*/ true, /*lb*/ 0, /*ub*/ 1,
                                              DoroDecoder<>::PursuitChoice::L1};
  decoder.stage_decode(&doro, rand_vec, &config, result);
  cout << format("Nonzero: {}\n MAE: {}\n", doro.nonzero_num(), doro.mae());
  return 0;
}