#pragma once

#include <vector>
#include "../global.h"
#include <ranges>

namespace ONIAK
{

  template <typename DT>
  class Orange
  {
  public:
    Orange(DT start, DT end, int size, bool include_start = true, bool include_end = true) : size_(size), cur_(0)
    {
      assert(size >= 1);
      double num_steps = size - 1;
      num_steps = (include_start) ? num_steps : num_steps + 1;
      num_steps = (include_end) ? num_steps : num_steps + 1;
      // normalize to left-closed, right-open ranges
      double step = (end - start) / num_steps;
      start_ = (include_start) ? start : start + step;
      end_ = (include_end) ? end + step : end;
    }

    DT operator()()
    {
      if (cur_ >= size_)
      {
        // range exceeded.
        return end_;
      }

      // future: check for overflow at integer multiplication
      DT result = start_ + (end_ - start_) * cur_ / size_;
      cur_++;
      return result;
    }

  private:
    DT start_, end_;
    int size_, cur_;

    friend std::vector<DT> orange_to_vector(Orange orange) {
      std::vector<DT> result;
      result.reserve(orange.size_);
      for (int i = 0; i < orange.size_; i++)
      {
        result.push_back(orange());
      }
      return result;
    }
  };


  template <typename Container>
  inline auto nonzero_indices(const Container& c) {
    return std::views::enumerate(c) | std::views::filter([](const auto& p) { return std::get<1>(p) != 0; }) | std::views::keys;
  }
}

