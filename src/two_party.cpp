#include "doro.hpp"
#include "decoder.hpp"
#include "probability.hpp"
#include "rans_wrapper.hpp"

#include <oniakDataStructure/ohist.h>
#include "oniakDebug/odebug.h"
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
enum class Status { CollisionAvoiding = 0, CollisionResolving = 1, Finished = 2 };
using CounterType = int8_t;

bool get_sizes(const std::unordered_map<int, CounterType>& result1,
  const std::unordered_map<int, CounterType>& result2, const std::unordered_set<int>& setA_minus_B, const std::unordered_set<int>& setB_minus_A,
  json& config) {
  int A_minus_B_remaining_size = setA_minus_B.size(), B_minus_A_remaining_size = setB_minus_A.size(), A_intersect_B_remaining_size = 0;
  for (auto [key, value] : result1) {
    if (value == 0) continue;
    if (setA_minus_B.contains(key)) --A_minus_B_remaining_size;
    else if (setB_minus_A.contains(key)) --B_minus_A_remaining_size;
    else {
      ++A_intersect_B_remaining_size;
    }
  }
  for (auto [key, value] : result2) {
    if (value == 0) continue;
    if (setA_minus_B.contains(key)) --A_minus_B_remaining_size;
    else if (setB_minus_A.contains(key)) --B_minus_A_remaining_size;
    else ++A_intersect_B_remaining_size;
  }
  config["A minus B remaining"].push_back(A_minus_B_remaining_size);
  config["B minus A remaining"].push_back(B_minus_A_remaining_size);
  config["A intersect B remaining"].push_back(A_intersect_B_remaining_size);
  return A_minus_B_remaining_size == 0 && B_minus_A_remaining_size == 0 && A_intersect_B_remaining_size == 0;
}

template <typename T>
Skellam moment_fit_skellam(const DoroCode<T>& code) {
  auto [mean, variance] = ONIAK::sample_mean_variance(code.code());
  double mu1, mu2;
  if (code.is_cbf()) {
    mu1 = (mean + variance) / 2.0;
    mu2 = (variance - mean) / 2.0;
    if (mu1 < 0) mu1 = 0;
    if (mu2 < 0) mu2 = 0;
  }
  else {
    if (variance < 0) variance = 0;
    mu1 = mu2 = variance / 2.0;
  }
  return { mu1, mu2 };
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    println("Error: please specify config json file.");
    exit(1);
  }
  std::ifstream fin(argv[1]);
  json config = json::parse(fin);

  size_t seed = 0;
  if (config.contains("seed")) seed = config.at("seed");

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
  int universe = config.at("universe");
  if (universe < A_union_B_size) {
    cout << "Warning: universe size is less than A_union_B_size. Reset to A_union_B_size." << endl;
    universe = A_union_B_size;
  }

  int k = config.at("k");
  int d = config.at("d");
  int s = config.at("s");   // signature size
  int s2 = s;               // resolving signature size
  if (config.contains("s2")) s2 = config.at("s2");
  int tk = 20;  // stage decode
  if (config.contains("tk")) tk = config.at("tk");
  int ta = 5;
  if (config.contains("ta")) ta = config.at("ta");
  int to = -1;
  if (config.contains("to")) to = config.at("to");
  int lb = 0;
  if (config.contains("lb")) lb = config.at("lb");
  int ub = 0;
  if (config.contains("ub")) ub = config.at("ub");
  std::string save_path = config.at("result filename");
  int max_rounds = 1'0000'0000;
  if (config.contains("max rounds")) max_rounds = config.at("max rounds");
  int max_comm_rounds = 20;
  if (config.contains("max comm rounds")) max_comm_rounds = config.at("max comm rounds");
  int resolving_round = max_comm_rounds - 5;  // at this round, always start resolving
  if (config.contains("resolving round")) resolving_round = config.at("resolving round");
  bool counting = config.at("counting");

  mt19937 rng(seed);
  auto rand_vec = random_nonrepetitive<int>(universe, A_union_B_size, rng);
  std::shuffle(rand_vec.begin(), rand_vec.end(), rng);
  // set A is [0, A_size), and set B is [A_minus_B_size, A_union_B_size)
  unordered_set<int> setA(rand_vec.begin(), rand_vec.begin() + A_size);
  unordered_set<int> setB(rand_vec.begin() + A_minus_B_size, rand_vec.end());
  unordered_set<int> setA_minus_B(rand_vec.begin(), rand_vec.begin() + A_minus_B_size);
  unordered_set<int> setB_minus_A(rand_vec.begin() + A_size, rand_vec.end());

  DoroCode<CounterType> doro(d, k, counting, rng, lb, ub);
  DoroCode<CounterType> first_round_code(d, k, counting, rng, lb, ub);
  unordered_map<int, CounterType> ground_truth, first_round_gt;
  for (int i : rand_vec | views::take(A_size)) {
    first_round_gt[i] = 1;
  }
  for (int i : rand_vec | views::take(A_minus_B_size)) {
    ground_truth[i] = -1;  // in Alis but not in Bela
  }
  for (auto i : rand_vec | views::drop(A_size)) {
    ground_truth[i] = 1;  // in Bela but not in Alis
  }
  first_round_code.encode(std::move(first_round_gt));
  doro.encode(std::move(ground_truth));

  // recenters all codes to range
  for (auto& val : doro.code()) {
    val = doro.recenter(val);
  }
  for (auto& val : first_round_code.code()) {
    val = first_round_code.recenter(val);
  }

  // rANS encode for uniform distribution in [lb, ub)
  auto uniform_map = uniform_pmf<CounterType>(lb, ub);
  RansWrapper first_round_wrapper(uniform_map);
  RansCode first_round_compressed_code = first_round_wrapper.encode(first_round_code.code());
  auto first_round_decompressed_code = first_round_wrapper.decode(first_round_compressed_code);

  double lambda = static_cast<double>(A_minus_B_size) * k / d;
  auto A_minus_B_map = get_pmf(lambda, d, counting);
  lambda = static_cast<double>(A_intersect_B_size) * k / d;
  auto A_intersect_B_map = get_pmf(lambda, d, counting);
  lambda = static_cast<double>(B_minus_A_size) * k / d;
  auto B_minus_A_map = get_pmf(lambda, d, counting);
  lambda = static_cast<double>(B_size) * k / d;
  auto B_map = get_pmf(lambda, d, counting);
  double B_entropy = entropy(B_map);
  double all_entropy = doro_entropy(A_intersect_B_map, A_minus_B_map, B_minus_A_map);
  // This is conditional entropy H(A | B).
  double first_round_entropy_cost = (all_entropy - B_entropy) * d;
  double first_round_cost = first_round_compressed_code.size() * 8;

  DEBUG_VECTOR_EQUAL(first_round_decompressed_code, first_round_code.code());
  config["comm costs"] = json::array({ first_round_cost });
  config["doro costs"] = json::array({ first_round_cost });
  config["theoretical entropy costs"] = json::array({ first_round_entropy_cost });
  config["finger costs"] = json::array({ 0 });
  config["resolving costs"] = json::array({ 0 });
  config["time"] = json::array({});
  config["num peels"] = json::array({});
  config["A minus B remaining"] = json::array({});
  config["B minus A remaining"] = json::array({});
  config["A intersect B remaining"] = json::array({});
  config["number of recenters"] = 0;

  Party party = Party::Alis;
  if (s > 31) {
    println("Warning: signature size is too large. Reset to 31.");
    s = 31;
  }
  unsigned int mask = (1 << s) - 1;
  WYHash finger_hash(mask, /*mode*/ 0, /*seed*/ rng());
  if (s2 > 31) {
    println("Warning: resolving signature size is too large. Reset to 31.");
    s2 = 31;
  }
  unsigned int mask2 = (1 << s2) - 1;
  WYHash resolving_hash(mask2, /*mode*/ 0, /*seed*/ rng());

  Status status = Status::CollisionAvoiding;
  DoroDecoder<CounterType> decoder_alis(&finger_hash, &resolving_hash), decoder_bela(&finger_hash, &resolving_hash);
  DecodeConfig dconf_alis(tk, max_rounds, to, ta, /*verbose*/ false, /*debug*/ false, /*lb*/ -1, /*ub*/ 0, PursuitChoice::L2),
    dconf_bela(tk, max_rounds, to, ta, /*verbose*/ false, /*debug*/ false, /*lb*/ 0, /*ub*/ 1, PursuitChoice::L2);
  StopWatch sw;
  unordered_map<int, CounterType> result_alis, result_bela;
  // elements whose fingerprints have been transmitted
  unordered_set<int> elements_alis, elements_bela;
  bool no_advance = false;
  int comm_rounds = 0, actual_comm_rounds = 0;
  bool success = false;

  while (!doro.empty() && comm_rounds < max_comm_rounds && status != Status::Finished) {
    actual_comm_rounds = std::max(actual_comm_rounds, comm_rounds + 1);
    // who is current decoder?
    party = (party == Party::Alis) ? Party::Bela : Party::Alis;
    const DecodeConfig& dconf = (party == Party::Alis) ? dconf_alis : dconf_bela;
    DoroDecoder<CounterType>& decoder = (party == Party::Alis) ? decoder_alis : decoder_bela,
      & other_decoder = (party == Party::Alis) ? decoder_bela : decoder_alis;
    unordered_set<int>& candidates = (party == Party::Alis) ? setA : setB;
    unordered_map<int, CounterType>* result;
    unordered_map<int, CounterType>& last_result = (party == Party::Alis) ? result_alis : result_bela;
    unordered_set<int>& elements = (party == Party::Alis) ? elements_alis : elements_bela;
    int num_peels = decoder.decode(&doro, candidates, &dconf, result);

    if (doro.empty() && decoder.unresolved_elements().empty()) {
      config["time"].push_back(sw.peek());
      config["num peels"].push_back(num_peels);
      success = success || get_sizes(decoder.result(), other_decoder.result(), setA_minus_B, setB_minus_A, config);
      break;
    }

    Skellam code_distribution = moment_fit_skellam(doro);
    auto code_frequency = code_distribution.pmf_map<CounterType>();
    double doro_cost = 0;
    RansWrapper rans_wrapper(code_frequency);
    RansCode doro_compressed_code = rans_wrapper.encode(doro.code());
    auto decompressed_code = rans_wrapper.decode(doro_compressed_code);

    double skellam_entropy = entropy(frequency_count(doro.code()));
    double theoretical_entropy = skellam_entropy * d;
    doro_cost = doro_compressed_code.size() * 8 + 64;
    // 8 bytes to transmit two Skellam parameters in float.
    config["theoretical entropy costs"].push_back(theoretical_entropy);
    config["doro costs"].push_back(doro_cost);
    DEBUG_VECTOR_EQUAL(decompressed_code, doro.code());

    int new_extra_size = 0;  // number of new elements decoded in this round
    for (auto [key, value] : *result) {
      if (value == 0) continue;
      if (!elements.contains(key)) {
        ++new_extra_size;
        elements.insert(key);
      }
    }

    int unresolved_size = decoder.unresolved_elements().size();
    double resolving_cost = unresolved_size * (s2);
    if (!decoder.unresolved_elements().empty()) {
      int num_unresolved = decoder.unresolved_elements().size();
      actual_comm_rounds = std::max(actual_comm_rounds, comm_rounds + 2); // need an additional resolving round
      int number_collisions = other_decoder.resolve_collision(decoder, setA, setB);
      if (number_collisions > 0) {
        resolving_cost += num_unresolved;
        actual_comm_rounds = std::max(actual_comm_rounds, comm_rounds + 3);   // resolving round needs feedback
      }
    }
    double finger_cost = new_extra_size * s;
    config["finger costs"].push_back(finger_cost);
    config["resolving costs"].push_back(resolving_cost);
    double comm_cost = doro_cost + finger_cost + resolving_cost;
    config["comm costs"].push_back(comm_cost);
    success = success || get_sizes(decoder.result(), other_decoder.result(), setA_minus_B, setB_minus_A, config);

    config["time"].push_back(sw.peek());
    config["num peels"].push_back(doro.num_peels());
    config["number of recenters"] = doro.num_recenters();

    if (decoder.result() == last_result || comm_rounds == resolving_round) {
      if (status == Status::CollisionAvoiding) {
        status = Status::CollisionResolving;
        decoder_alis.enter_resolving();
        decoder_bela.enter_resolving();
        config["round entering resolving"].push_back(actual_comm_rounds);
      }
      else if (!no_advance) no_advance = true;
      else if (decoder.result() == last_result) status = Status::Finished;
    }
    else no_advance = false;

    decoder.update_fingerprints(other_decoder);
    last_result = decoder.result();

    ++comm_rounds;
  }
  config["communication rounds"] = actual_comm_rounds;

  double total_comm_cost = 0;
  for (auto cost : config.at("comm costs")) {
    total_comm_cost += cost.get<double>();
  }
  config["total communication cost"] = total_comm_cost;
  config["success"] = success;

  std::ofstream fout(save_path);
  fout << std::setw(4) << config << endl;  // indent 4
  return 0;
}