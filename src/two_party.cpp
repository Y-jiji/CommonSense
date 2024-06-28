#include "doro.hpp"
#include "decoder.hpp"
#include "probability.hpp"

#include "oniakRandom/orand.h"
#include "oniakTimer/otime.h"
#include "nlohmann/json.hpp"

#include <format>
#include <iostream>
#include <random>
#include <ranges>
#include <set>
#include <unordered_map>

using namespace std;
using namespace Doro;
using namespace ONIAK;
using namespace nlohmann;

int main()
{
  mt19937 rng(0);

  Skellam skellam = create_corrected_skellam(1, 1000);
  Poisson poisson(1);
  std:: cout << "Skellam: " << entropy(skellam.pmf_map()) << std::endl;
  std:: cout << "Poisson: " << entropy(poisson.pmf_map()) << std::endl;

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

  return 0;
}