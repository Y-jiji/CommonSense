#ifndef __PROBABILITY_HPP__
#define __PROBABILITY_HPP__

#include <cmath>
#include <unordered_map>
#include <vector>

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

Skellam create_corrected_skellam(double lambda, double x) {
  double mu = lambda * x / (2 * x - 1);
  return {mu, mu};
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

// all entropies are in bits
template <typename T>
double entropy(const std::vector<T>& vec) {
  double sum = 0;
  for (auto p : vec) {
    sum += p * std::log2(p);
  }
  return -sum;
}

template <typename T, typename V>
double entropy(const std::unordered_map<T, V>& map) {
  double sum = 0;
  for (auto [key, p] : map) {
    sum += p * std::log2(p);
  }
  return -sum;
}

}

#endif
