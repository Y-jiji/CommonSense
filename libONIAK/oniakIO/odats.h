#pragma once

#include <fstream>
#include <iostream>
#include <type_traits>

#include "../global.h"
#include "../oniakPath/opath.h"
#include "odat_type.h"

namespace ONIAK
{

  // read a file in datatype FT and returns an eigen matrix of type T
  // returns an empty matrix if # columns in each row is not equal.
  template <typename T, typename FT = T>
  EigenMatrixRowMaj<T> read_file(std::ifstream &fin, size_t file_sz)
  {
    size_t type_sz = sizeof(FT);
    int sz;
    fin.read(reinterpret_cast<char *>(&sz), 4);
    int64_t cols = sz;
    int64_t bytes_per_line = cols * type_sz + 4;
    int64_t rows = file_sz / bytes_per_line;

    // test if division is complete
    assert(rows * bytes_per_line == static_cast<int64_t>(file_sz));

    EigenMatrixRowMaj<T> result(rows, cols);
    if constexpr (std::is_same_v<T, FT>)
    {
      // read without type cast
      for (int rr = 0; rr < rows - 1; ++rr)
      {
        fin.read(reinterpret_cast<char *>(&result(rr, 0)), type_sz * cols);

        fin.read(reinterpret_cast<char *>(&sz), 4);
        // abort if # cols not equal.
        if (sz != cols)
          return {};
      }
      fin.read(reinterpret_cast<char *>(&result(rows - 1, 0)), type_sz * cols);
    }
    else
    {
      // type cast is needed
      RowVector<FT> buffer(cols);
      for (int rr = 0; rr < rows - 1; ++rr)
      {
        fin.read(reinterpret_cast<char *>(buffer.data()), type_sz * cols);
        result.row(rr) = buffer.template cast<T>();

        fin.read(reinterpret_cast<char *>(&sz), 4);
        // abort if # cols not equal.
        if (sz != cols)
          return {};
      }
      fin.read(reinterpret_cast<char *>(buffer.data()), type_sz * cols);
      result.row(rows - 1) = buffer.template cast<T>();
    }

    return result;
  }

  template <typename T>
  EigenMatrixRowMaj<T> read_matrix(std::string filename)
  {
    auto extension = file_extension(filename);
    size_t file_sz = std::filesystem::file_size(filename);
    std::ifstream fin(filename);
    if (!fin)
    {
      std::cerr << "[Error] file " << filename << " is not existent!" << std::endl;
      return {};
    }

    if (extension == "fvecs")
    {
      return read_file<T, float>(fin, file_sz);
    }
    else if (extension == "ivecs")
    {
      return read_file<T, int32_t>(fin, file_sz);
    }
    else if (extension == "dvecs")
    {
      return read_file<T, double>(fin, file_sz);
    }
    else if (extension == "bvecs")
    {
      return read_file<T, uint8_t>(fin, file_sz);
    }
    else if (extension == "odat")
    {
      file_sz--;
      unsigned char dtype;
      fin.read(reinterpret_cast<char *>(&dtype), 1);

      switch (dtype)
      {
      case 0x51: // int32
        return read_file<T, int32_t>(fin, file_sz);
      case 0x52: // float
        return read_file<T, float>(fin, file_sz);
      case 0x61: // int64
        return read_file<T, int64_t>(fin, file_sz);
      case 0x62:
        return read_file<T, double>(fin, file_sz);
      case 0x33:
        return read_file<T, uint8_t>(fin, file_sz);
      default:
        return {};
      }
    }
    else
    { // unsupported format
      return {};
    }
  }

  template <typename FT, typename T>
  void write_file(std::string filename, const EigenMatrixRowMaj<T> &input)
  {
    std::ofstream fout(filename);
    unsigned char dtype = ODAT_code<FT>::value;
    fout.write(reinterpret_cast<char *>(&dtype), 1);

    EigenMatrixRowMaj<FT> converted_matrix;
    const EigenMatrixRowMaj<FT> *cview;
    if constexpr (std::is_same_v<T, FT>)
    {
      cview = &input;
    }
    else
    {
      converted_matrix = input.template cast<FT>();
      cview = &converted_matrix;
    }
    const EigenMatrixRowMaj<FT> &matrix = *cview;

    size_t type_sz = sizeof(FT);
    int sz = matrix.cols();
    for (int rr = 0; rr < matrix.rows(); ++rr)
    {
      fout.write(reinterpret_cast<char *>(&sz), 4);
      fout.write(reinterpret_cast<const char *>(&matrix(rr, 0)), type_sz * sz);
    }
  }

  // read until file end or given number of lines
  // when FT != T, assume implicit conversion is possible
  template <typename T, typename FT = T>
  DoubleVector<T> read_nonuniform_file(std::ifstream &fin, int lines)
  {
    size_t type_sz = sizeof(FT);
    DoubleVector<T> result;

    while (true)
    {
      int sz;
      fin.read(reinterpret_cast<char *>(&sz), 4);
      if (!fin || lines <= 0)
        break; // end of file or enough lines

      std::vector<FT> buffer;
      std::vector<T> vec;
      fin.read(reinterpret_cast<char *>(buffer.data()), type_sz * sz);
      if constexpr (std::is_same_v<T, FT>)
      {
        vec = std::move(buffer);
      }
      else
      {
        vec.assign(buffer.begin(), buffer.end());
      }
      result.emplace_back(std::move(vec));
      --lines;
    }
    return result;
  }

  // Read non uniform file into vector<vector<T>> 
  template <typename T>
  DoubleVector<T> read_double_vector(std::string filename, int lines = INT_MAX)
  {
    auto extension = file_extension(filename);
    std::ifstream fin(filename);
    if (!fin)
    {
      std::cerr << "[Error] file " << filename << " is not existent!" << std::endl;
      return {};
    }

    if (extension == "fvecs")
    {
      return read_nonuniform_file<T, float>(fin, lines);
    }
    else if (extension == "ivecs")
    {
      return read_nonuniform_file<T, int32_t>(fin, lines);
    }
    else if (extension == "dvecs")
    {
      return read_nonuniform_file<T, double>(fin, lines);
    }
    else if (extension == "bvecs")
    {
      return read_nonuniform_file<T, uint8_t>(fin, lines);
    }
    else if (extension == "odat")
    {
      unsigned char dtype;
      fin.read(reinterpret_cast<char *>(&dtype), 1);

      switch (dtype)
      {
      case 0x51: // int32
        return read_nonuniform_file<T, int32_t>(fin, lines);
      case 0x52: // float
        return read_nonuniform_file<T, float>(fin, lines);
      case 0x61: // int64
        return read_nonuniform_file<T, int64_t>(fin, lines);
      case 0x62:
        return read_nonuniform_file<T, double>(fin, lines);
      case 0x33:
        return read_nonuniform_file<T, uint8_t>(fin, lines);
      default:
        return {};
      }
    }
    else
    { // unsupported format
      return {};
    }
  }

  // Write vector<vector<T>> (no need to be uniform) into file
  template <typename FT, typename T>
  void write_file(std::string filename, const DoubleVector<T> &input)
  {
    std::ofstream fout(filename);
    unsigned char dtype = ODAT_code<FT>::value;
    fout.write(reinterpret_cast<char *>(&dtype), 1);
    size_t type_sz = sizeof(FT);

    for (const auto &vec : input)
    {
      int sz = vec.size();
      fout.write(reinterpret_cast<char *>(&sz), 4);
      std::vector<FT> buffer;
      const std::vector<FT> *view;

      if constexpr (std::is_same_v<T, FT>)
      {
        view = &vec;
      }
      else
      {
        view = &buffer;
        buffer.assign(vec.begin(), vec.end());
      }
      fout.write(reinterpret_cast<const char *>(view->data()), type_sz * sz);
    }
  }

  template <typename FT, typename T>
  void write_file(std::string filename, const std::vector<T> &input)
  {
    DoubleVector<T> double_vec = {input};
    write_file<FT>(filename, double_vec);
  }

} // namespace ONIAK

