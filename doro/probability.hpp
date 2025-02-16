#ifndef __PROBABILITY_HPP__
#define __PROBABILITY_HPP__

#include <numeric>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>
#include <iostream>

#include "libONIAK/oniakDataStructure/ohist.h"
#include "libONIAK/oniakMath/ostats.h"

namespace Doro {
constexpr double DoubleEpsilon = 1e-18;

class Skellam {
public:
  Skellam(double mu1, double mu2) : mu1_(mu1), mu2_(mu2) {}

  // cyn_bessel_i is called "modified Bessel function of first type" in Wikipedia
  double pmf(int k) const {
    if (mu1_ < DoubleEpsilon) {
      if (k > 0) return std::numeric_limits<double>::quiet_NaN();
      else return std::exp(-mu2_ - k * std::log(mu2_) - std::lgamma(-k + 1));
    }
    else if (mu2_ < DoubleEpsilon) {
      if (k < 0) return std::numeric_limits<double>::quiet_NaN();
      else return std::exp(-mu1_ + k * std::log(mu1_) - std::lgamma(k + 1));
    }
    return std::exp(-mu1_ - mu2_) * std::pow(mu1_ / mu2_, static_cast<double>(k) / 2) * std::cyl_bessel_i(std::abs(k), 2 * sqrt(mu1_ * mu2_));
  }

  // skellam may have negative values
  template <typename T = int>
  std::unordered_map<T, double> pmf_map() const {
    if (mu1_ < DoubleEpsilon && mu2_ < DoubleEpsilon) return { { 0, 1.0 } };
    std::unordered_map<T, double> result;
    int mean = mu1_ - mu2_;
    for (k1_ = mean; true; ++k1_) {
      double p = pmf(k1_);
      if (p < DoubleEpsilon || std::isnan(p))
        break;
      else
        result[k1_] = p;
    }
    for (k2_ = mean - 1; true; --k2_) {
      double p = pmf(k2_);
      if (p < DoubleEpsilon || std::isnan(p))
        break;
      else
        result[k2_] = p;
    }
    return result;
  }

  // cdf, and both ends
  // cdf[k] = Pr(X < k), so Pr(X in [k1, k2)) = cdf[k2] - cdf[k1]
  std::tuple<std::unordered_map<int, double>, int, int>
    cdf_map() const {
    auto pmf = pmf_map<int>();
    std::unordered_map<int, double> cdf;
    double sum = 0;
    for (int k = k2_ + 1; k < k1_; ++k) {
      cdf[k] = sum;
      sum += pmf[k];
    }
    cdf[k1_] = sum;
    return { cdf, k2_, k1_ };
  }

private:
  double mu1_, mu2_;
  mutable int k1_, k2_; // high and low ends of distribution
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
    if (lambda_ < DoubleEpsilon) return { {0, 1.0} };
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
  }
  else {
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

inline double entropy(double p) {
  return -p * std::log2(p) - (1 - p) * std::log2(1 - p);
}

// computes entropy of (X+Y), (X+Z)
inline double doro_entropy(const std::unordered_map<int, double>& rvx, const std::unordered_map<int, double>& rvy,
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

inline int log2_ceil(int x) {
  return static_cast<int>(std::ceil(std::log2(x)));
}

// both input and output is double
template <typename F>
concept UnitaryFunction = std::is_invocable_r_v<double, F, double>;

template <UnitaryFunction F>
double convex_argmax(F f, double lb, double ub, double epsilon = 1e-6) {
  double mid = (lb + ub) / 2;
  double fmid = f(mid);
  double flb = f(lb);
  double fub = f(ub);

  while (fmid < flb && ub - lb > epsilon) {
    ub = mid;
    mid = (lb + ub) / 2;
    fmid = f(mid);
  }

  while (fmid < fub && ub - lb > epsilon) {
    lb = mid;
    mid = (lb + ub) / 2;
    fmid = f(mid);
  }

  while (ub - lb > epsilon) {
    double& farther = (ub - mid > mid - lb) ? ub : lb;
    double& closer = (ub - mid > mid - lb) ? lb : ub;
    double new_mid = (mid + farther) / 2;
    double fnew_mid = f(new_mid);
    if (fnew_mid > fmid) {
      closer = mid;
      mid = new_mid;
      fmid = fnew_mid;
    }
    else {
      farther = new_mid;
    }
  }
  return mid;
}

// average load per fingerprint bucket is |A/B| * |B/A| / |C|^2, which is alpha * beta^2
// however, the target epsilon per bucker is epsilon * beta. This cancels out one copy of beta.
inline int finger_l_size(double alpha, double beta, double epsilon) {
  if (alpha < 1e-30 || alpha > 1e30) return 0;
  return std::ceil(-std::log2(epsilon) + std::log2(alpha) + std::log2(beta));
}

// alpha is |A/B| / |B/A|, beta is |B/A| / |C|
inline double loss_func(double alpha, double beta, double b, double epsilon) {
  double first_nonzero_prob = 1 - std::exp(-beta);
  double entropy_first = entropy(first_nonzero_prob);
  double second_nonzero_prob = 1 - std::exp(-alpha * beta);
  double entropy_second = entropy(second_nonzero_prob * first_nonzero_prob);
  double recon_cost = alpha * beta * first_nonzero_prob * (finger_l_size(alpha, beta, epsilon) + log2_ceil(b / beta));
  // std::cout << "entropy_first: " << entropy_first/beta << ", entropy_second: " << entropy_second/beta << ", recon_cost: " << recon_cost/beta << std::endl;
  return (entropy_first + entropy_second + recon_cost) / beta;
}

}  // namespace Doro

#endif
