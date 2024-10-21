#ifndef __PROBABILITY_HPP__
#define __PROBABILITY_HPP__

#include <numeric>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

#include "oniakDataStructure/ohist.h"
#include "oniakMath/ostats.h"

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
  template <typename T = int>
  std::unordered_map<T, double> pmf_map() const {
    if (mu1_ < DoubleEpsilon && mu2_ < DoubleEpsilon) return {{ 0, 1.0 }};
    std::unordered_map<T, double> result;
    for (int k = 0; true; ++k) {
      double p = pmf(k);
      if (p < DoubleEpsilon || std::isnan(p))
        break;
      else
        result[k] = p;
    }
    for (int k = -1; true; --k) {
      double p = pmf(k);
      if (p < DoubleEpsilon || std::isnan(p))
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
    return std::exp(-lambda_ + k * std::log(lambda_) - std::lgamma(k + 1));
  }

  template <typename T = int>
  std::unordered_map<T, double> pmf_map() const {
    if (lambda_ < DoubleEpsilon) return {{0, 1.0}};
    std::unordered_map<T, double> result;
    bool started = false;
    for (int k = 0; true; ++k) {
      double p = pmf(k);
      if (started && (p < DoubleEpsilon || std::isnan(p)))
        break;
      if (p >= DoubleEpsilon)
        started = true;
      result[k] = p;
    }
    return result;
  }

private:
  double lambda_;
};

/* get pmf of counters*/
std::unordered_map<int, double> get_pmf(double lambda, int x, bool counting) {
  if (counting) {
    Poisson poisson(lambda);
    return poisson.pmf_map();
  } else {
    Skellam skellam = create_corrected_skellam(lambda, x);
    return skellam.pmf_map();
  }
}

// uniform distribution of [lb, ub)
template<typename T>
std::unordered_map<T, double> uniform_pmf(int lb, int ub) {
  std::unordered_map<T, double> result;
  for (int i = lb; i < ub; ++i) {
    result[i] = 1.0 / (ub - lb);
  }
  return result;
}

// all entropies are in bits
template <typename T, typename V, typename H>
double entropy(const std::unordered_map<T, V, H>& map) {
  double psum = std::accumulate(map.begin(), map.end(), 0.0, [](double sum, const auto& pair) {
    return sum + pair.second;
  });
  double sum = 0;
  for (auto [key, p] : map) {
    double p0 = p / psum;
    sum += p0 * std::log2(p0);
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

int log2_ceil(int x) {
  return static_cast<int>(std::ceil(std::log2(x)));
}

}  // namespace Doro

#endif
