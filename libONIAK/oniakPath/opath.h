#pragma once

#include <cassert>
#include <filesystem>
#include <string>

namespace ONIAK {

inline std::string file_extension(const std::string& path) {
  size_t dot_pos = path.find_last_of(".");
  assert(dot_pos != std::string::npos && dot_pos+1 < path.size());
  return path.substr(dot_pos+1);
}

// joins path
inline std::string join_path(std::filesystem::path p1, std::filesystem::path p2) {
    return p1 / p2;
}

template<typename... Targs>
std::string join_path(std::filesystem::path p1, std::filesystem::path p2, Targs... Fargs) {
  return join_path(p1/p2, Fargs...);
}

}

