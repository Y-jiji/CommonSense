#ifndef __PROBABILITY_HPP__
#define __PROBABILITY_HPP__

#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

#include "oniakDataStructure/ohist.h"

namespace Doro {
constexpr double DoubleEpsilon = 1e-18;

class Skellam {
public:
  Skellam(double mu1, double mu2) : mu1_(mu1), mu2_(mu2) {}

  // cyn_bessel_i is called "modified Bessel function of first type" in Wikipedia
  double pmf(int k) const {
    return std::exp(-mu1_ - mu2_) * std::pow(mu1_ / mu2_, k / 2) * std::cyl_bessel_i(std::abs(k), 2 * sqrt(mu1_ * mu2_));
  }

  // skellam may have negative values
  std::unordered_map<int, double> pmf_map() const {
    std::unordered_map<int, double> result;
    for (int k = 0; true; ++k) {
      double p = pmf(k);
      if (p < DoubleEpsilon)
        break;
      else
        result[k] = p;
    }
    for (int k = -1; true; --k) {
      double p = pmf(k);
      if (p < DoubleEpsilon)
        break;
      else
        result[k] = p;
    }
    return result;
  }

private:
  double mu1_, mu2_;
};

// correlation is because an item is only hashed to one bucket per hash function
Skellam create_corrected_skellam(double lambda, double x) {
  double mu = lambda * x / (2 * x - 1);
  return { mu, mu };
}

class Poisson {
public:
  Poisson(double lambda) : lambda_(lambda) {}

  double pmf(int k) const {
    return std::exp(-lambda_) * std::pow(lambda_, k) / std::tgamma(k + 1);
  }

  std::unordered_map<int, double> pmf_map() const {
    std::unordered_map<int, double> result;
    for (int k = 0; true; ++k) {
      double p = pmf(k);
      if (p < DoubleEpsilon)
        break;
      else
        result[k] = p;
    }
    return result;
  }

private:
  double lambda_;
};

std::unordered_map<int, double> get_pmf(double lambda, int x, bool counting) {
  if (counting) {
    Poisson poisson(lambda);
    return poisson.pmf_map();
  } else {
    Skellam skellam = create_corrected_skellam(lambda, x);
    return skellam.pmf_map();
  }

}

// all entropies are in bits
template <typename T, typename V, typename H>
double entropy(const std::unordered_map<T, V, H>& map) {
  double sum = 0, psum = 0;
  for (auto [key, p] : map) {
    sum += p * std::log2(p);
    psum += p;
  }
  if (std::abs(psum - 1) > 1e-6) {
    std::cerr << "Warning: Probability does not sum to 1." << std::endl;
  }
  return -sum;
}

// computes entropy of (X+Y), (X+Z)
double doro_entropy(const std::unordered_map<int, double>& rvx, const std::unordered_map<int, double>& rvy,
  const std::unordered_map<int, double>& rvz) {
  using int_pair = std::pair<int, int>;
  std::unordered_map<int_pair, double, ONIAK::TupleHash<int, int>> pmap;

  for (auto [vx, px] : rvx) {
    for (auto [vy, py] : rvy) {
      for (auto [vz, pz] : rvz) {
        int_pair key = { vx + vy, vx + vz };
        auto iter = pmap.find(key);
        double p = (iter != pmap.end()) ? iter->second : 0;
        p += px * py * pz;
        pmap[key] = p;
      }
    }
  }
  return entropy(pmap);
}

template<typename T>
int num_extra_elements(const std::unordered_map<int, T>& map1, const std::unordered_map<int, T>& map2) {
  int count = 0;
  for (auto [key, value] : map1) {
    if (map2.find(key) == map2.end() || map2.at(key) != value)
      count += 1;
  }
  return count;
}

}  // namespace Doro

#endif
