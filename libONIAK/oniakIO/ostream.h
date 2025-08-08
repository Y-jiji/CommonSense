#pragma once

#include <iostream>
#include <vector>
#include "../global.h"

namespace ONIAK {

template<typename T>
std::ostream& operator<<(std::ostream& os, const Eigen::Matrix<T, 1, Eigen::Dynamic>& row){
  for (auto val: row) {
    os << val << "\t";
  }
  return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& row){
  for (auto val: row) {
    os << val << "\t";
  }
  return os;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix){
  int rows = matrix.rows();
  for (int rr = 0; rr < rows; ++rr) {
    os << matrix.row(rr) << std::endl;
  }
  return os;
}

template<typename T>
void print_one_line(std::ostream& os, T value) {
    os << value << std::endl;
}

template<typename T, typename... Targs>
void print_one_line(std::ostream& os, T value, Targs... Fargs) {
    os << value << "\t";
    print_one_line(os, Fargs...);
}

// Prints table in csv format.
template<typename T>
void print_table(std::ostream& os, T value) {
    os << value << std::endl;
}

template<typename T, typename... Targs>
void print_table(std::ostream& os, T value, Targs... Fargs) {
    os << value << ",";
    print_table(os, Fargs...);
}

}

