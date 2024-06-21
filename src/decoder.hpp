#ifndef __DORO_DECODER_HPP__
#define __DORO_DECODER_HPP__

#include "doro.hpp"
#include "oniakDataStructure/oupq.h"

#include <algorithm>
#include <format>
#include <iostream>
#include <map>
#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace Doro
{
  template <typename ArrType = int32_t, ONIAK::UpdatePQBackend backend = ONIAK::UpdatePQBackend::PriorityQueue>
  class DoroDecoder
  {
  public:
    using DoroCodeT = DoroCode<ArrType>;
    using TwoDimVector = std::vector<std::vector<int>>;
    using UpdatePQ = ONIAK::UpdatePQ<backend>;

    enum class PursuitChoice
    {
      L1,
      L2
    };
    struct DecodeConfig
    {
      int tk, max_stage; // used for stage decode.
      // tk is stage size
      int ta; // terminate if any element thrashes for ta times.
      bool verbose;
      // lower and upper bound of value range
      double lb, ub;
      PursuitChoice pursuit_choice;
    };
    DoroDecoder() = default;

    int stage_decode(
        DoroCodeT *code,
        const std::vector<int> &setA, /*possible candidates of code*/
        DecodeConfig *config,
        std::unordered_map<int, ArrType> &result)
    {
      code_ = code;
      config_ = config;
      thrashing_.clear();
      result_.clear();
      bool is_l2 = config_->pursuit_choice == PursuitChoice::L2; // else is l1
      for (int element : setA)
      {
        thrashing_[element] = 0;
      }

      for (int rnd : std::views::iota(0, config_->max_stage))
      {
        std::vector<std::pair<ArrType, int>> signals;
        for (auto element_thrash : thrashing_)
        {
          int element = element_thrash.first;
          ArrType strength = is_l2 ? code_->sense(element) : code_->sense_l1(element);
          signals.emplace_back(strength, element);
        }
        std::nth_element(signals.begin(), signals.begin() + config_->tk, signals.end());
        signals.resize(config_->tk);
        bool finished = false;
        for (auto [strength, element] : signals)
        {
          auto cur_iter = result.find(element);
          ArrType cur_element_value = (cur_iter != result.end()) ? cur_iter->second : 0;
          // In L2, delta is computed from average instead of sum
          double normalized_strength = is_l2 ? static_cast<double>(strength) / static_cast<double>(code->k()) : strength;
          ArrType delta = get_delta(normalized_strength, cur_element_value);
          if (std::abs(delta) > 1e-6)
          {
            int thrash = thrashing_.find(element)->second;
            if (thrash > config_->ta)
            {
              finished = true;
              break;
            }
            thrashing_[element] = thrash + 1;
            code_->peel(element, delta);
            result[element] = cur_element_value + delta;

            if (config_->verbose)
            {
              std::cout << std::format("thrash: {}, power: {}, element: {}, cur_element_value: {}, delta: {}\n",
                                       thrash, strength, element, cur_element_value, delta);
            }
          }
        }
        if (config_->verbose)
        {
          std::cout << std::format("Round: {}\n", rnd);
          code_->show_result();
        }
        if (finished)
          return rnd;
      }
      return config_->max_stage;
    }

  private:
    ArrType
    get_delta(double delta, ArrType cur_value) const
    {
      double new_value = cur_value + delta;
      new_value = std::max(config_->lb, new_value);
      new_value = std::min(config_->ub, new_value);
      delta = new_value - cur_value;
      if constexpr (std::is_integral_v<ArrType>)
      {
        if (std::abs(delta) < 0.50001)
          delta = 0; // if delta is 0.5, default round to 1
        else         // but we let it be 0 to avoid oscillation
          delta = std::round(delta);
      }
      return delta;
    }

    DoroCodeT *code_;
    TwoDimVector neighbors_;
    std::unordered_map<int, ArrType> result_;
    std::unordered_map<int, int> thrashing_;
    UpdatePQ priority_queue_;
    DecodeConfig *config_;
  };
}
#endif