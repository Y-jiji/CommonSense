#pragma once

#include "../global.h"

namespace ONIAK {

template <typename DT> 
DT lpnorm(const RowVector<DT>& vector, double pow) {
  auto vec_pow = vector.array().abs().pow(pow);
  return std::pow(vec_pow.sum(), 1.0/pow);
}

}

