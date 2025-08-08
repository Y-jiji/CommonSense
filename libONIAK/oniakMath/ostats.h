#pragma once

#include <numeric>
#include <utility>
#include <vector>

namespace ONIAK {

template <typename T>
double sample_mean(const std::vector<T>& data) {
  return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

template <typename T>
std::pair<double, double> sample_mean_variance(const std::vector<T>& data) {
  double mean = sample_mean(data);
  double variance = std::accumulate(data.begin(), data.end(), 0.0, [mean](double acc, T x) {
    return acc + (x - mean) * (x - mean);
  }) / (data.size() - 1);
  return {mean, variance};
}

}
