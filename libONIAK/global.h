#pragma once

#include <cassert>
#include <vector>
#include <eigen3/Eigen/Dense>

namespace ONIAK {

template <typename T>
using EigenMatrixRowMaj = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

template <typename T>
using RowVector = Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
using ColumnVector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;

// Here double means two layers instead of data type "double"
template <typename T>
using DoubleVector = std::vector<std::vector<T>>;

// Always-false dependent type: used to nullify default template specializations
template <typename T>
struct AlwaysFalse : std::false_type {};

}
