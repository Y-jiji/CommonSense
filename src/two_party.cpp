#include "doro.hpp"
#include "decoder.hpp"
#include "probability.hpp"

#include "oniakRandom/orand.h"
#include "oniakTimer/otime.h"
#include "nlohmann/json.hpp"

#include <fstream>
#include <format>
#include <iostream>
#include <print>
#include <random>
#include <ranges>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace Doro;
using namespace ONIAK;
using namespace nlohmann;

enum class Party { Alis = 0, Bela = 1 };
enum class Status {CollisionAvoiding = 0, CollisionResolving = 1, Finished = 2};

int main(int argc, char* argv[]) {
  if (argc < 1) {
    println("Error: please specify config json file.");
    exit(1);
  }
  std::ifstream fin(argv[0]);
  json config = json::parse(fin);

  size_t seed = 0;
  if (config.contains("seed")) seed = config.at("seed");
  int universe = config.at("universe");
  int A_size = config.at("A size");
  int B_size = config.at("B size");
  int A_minus_B_size = A_size - B_size;
  if (config.contains("A minus B size")) A_minus_B_size = config.at("A minus B size");
  int A_union_B_size = B_size + A_minus_B_size;
  int A_intersect_B_size = A_size - A_minus_B_size;
  int B_minus_A_size = B_size - A_intersect_B_size;
  int k = config.at("k");
  int d = config.at("d");
  int s = config.at("s");   // signature size
  int tk = 20;  // stage decode
  if (config.contains("tk")) tk = config.at("tk");
  int ta = 5;
  if (config.contains("ta")) tk = config.at("ta");
  std::string save_path = config.at("result filename");
  int max_rounds = 1'0000'0000;
  if (config.contains("max rounds")) max_rounds = config.at("max rounds");
  int max_comm_rounds = 10;
  if (config.contains("max comm rounds")) max_comm_rounds = config.at("max comm rounds");
  bool counting = config.at("counting");

  mt19937 rng(seed);
  auto rand_vec = random_nonrepetitive<int>(universe, A_union_B_size, rng);
  std::shuffle(rand_vec.begin(), rand_vec.end(), rng);
  // set A is [0, A_size), and set B is [A_minus_B_size, A_union_B_size)
  std::vector<int> setA(rand_vec.begin(), rand_vec.begin() + A_size);
  std::vector<int> setB(rand_vec.begin() + A_minus_B_size, rand_vec.end());
  std::unordered_set<int> setA_minus_B(rand_vec.begin(), rand_vec.begin() + A_minus_B_size);
  std::unordered_set<int> setB_minus_A(rand_vec.begin() + A_size, rand_vec.end());

  DoroCode doro(d, k, counting, rng);
  unordered_map<int, int> ground_truth;
  for (int i : rand_vec | views::take(A_minus_B_size)) {
    ground_truth[i] = -1;  // in Alis but not in Bela
  }
  for (auto i : rand_vec | views::drop(A_size)) {
    ground_truth[i] = 1;  // in Bela but not in Alis
  }
  doro.encode(std::move(ground_truth));
  double lambda = static_cast<double>(A_minus_B_size) * k / d;
  auto A_minus_B_map = get_pmf(lambda, d, counting);
  lambda = static_cast<double>(A_intersect_B_size) * k / d;
  auto A_intersect_B_map = get_pmf(lambda, d, counting);
  lambda = static_cast<double>(B_minus_A_size) * k / d;
  auto B_minus_A_map = get_pmf(lambda, d, counting);
  double first_round_cost = doro_entropy(A_intersect_B_map, A_minus_B_map, B_minus_A_map) * d;
  config["comm costs"] = json::array({ first_round_cost });
  config["time"] = json::array({});
  config["A minus B remaining"] = json::array({});
  config["B minus A remaining"] = json::array({});
  config["A intersect B remaining"] = json::array({});

  Party party = Party::Alis;
  if (s > 31) {
    println("Warning: signature size is too large.");
    s = 31;
  }
  unsigned int mask = (1 << s) - 1;
  Status status = Status::CollisionAvoiding;
  WYHash finger_hash(/*mode*/ 0, mask, /*seed*/ rng());
  DoroDecoder<int32_t> decoder_alis, decoder_bela;
  DecodeConfig dconf_alis(tk, max_rounds, ta, /*verbose*/ false, /*debug*/ false, /*lb*/ 0, /*ub*/ 1, PursuitChoice::L2),
    dconf_bela(tk, max_rounds, ta, /*verbose*/ false, /*debug*/ false, /*lb*/ -1, /*ub*/ 0, PursuitChoice::L2);
  StopWatch sw;
  unordered_map<int, int> result_alis, result_bela;
  unordered_set<int> fingerprints;

 
  
  int comm_rounds = 0;

  while (!doro.isempty() && comm_rounds < max_comm_rounds && status != Status::Finished) {
    // who is current decoder?
    party = (party == Party::Alis) ? Party::Bela : Party::Alis;
    const DecodeConfig& dconf = (party == Party::Alis) ? dconf_alis : dconf_bela;
    DoroDecoder<int32_t>& decoder = (party == Party::Alis) ? decoder_alis : decoder_bela,
      other_decoder = (party == Party::Alis) ? decoder_bela : decoder_alis;
    const std::vector<int>& candidates = (party == Party::Alis) ? setA : setB;
    std::unordered_map<int, int>* result;
    std::unordered_map<int, int>& last_result = (party == Party::Alis) ? result_alis : result_bela;
    decoder.decode(&doro, candidates, &dconf, result);

    double empirical_entropy = entropy(frequency_count(doro.code()));
    int new_extra_size = num_extra_elements(*result, last_result);
    double comm_cost = empirical_entropy * d + new_extra_size * s;
    if (new_extra_size == 0 && status == Status::CollisionAvoiding) {
      status = Status::CollisionResolving;
      decoder_alis.enter_resolving();
      decoder_bela.enter_resolving();
    } else if (new_extra_size == 0 && status == Status::CollisionResolving) {
      status = Status::Finished;
    }

    config["comm costs"].push_back(comm_cost);
    decoder.update_fingerprints(other_decoder);
    last_result = *result;
    config["time"].push_back(sw.peek()); 

    int A_minus_B_remaining_size = 0, B_minus_A_remaining_size = 0, A_intersect_B_remaining_size = 0;
    for (auto [key, value] : *result) {
      if (value == 0) continue;
      if (setA_minus_B.contains(key)) ++A_minus_B_remaining_size;
      else if (setB_minus_A.contains(key)) ++B_minus_A_remaining_size;
      else ++A_intersect_B_remaining_size;
    }
    config["A minus B remaining"].push_back(A_minus_B_remaining_size);
    config["B minus A remaining"].push_back(B_minus_A_remaining_size);
    config["A intersect B remaining"].push_back(A_intersect_B_remaining_size);

    ++comm_rounds;
  }
  config["communication rounds"] = comm_rounds;

  std::ofstream fout(save_path);
  fout << config << endl;
  return 0;
}