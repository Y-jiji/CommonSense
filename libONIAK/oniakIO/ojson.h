#pragma once

#include "nlohmann/json.hpp"
#include "../oniakPath/opath.h"

#include <fstream>
#include <string>

namespace ONIAK {

inline nlohmann::json load_config(std::string filename) {
  using json = nlohmann::json;
  std::ifstream fin(filename);
  json config = json::parse(fin);

  std::string root_path = config["root path"];

  for (json::iterator it = config.begin(); it != config.end(); ++it) {
    std::string key = it.key();
    if (key.find("file") != std::string::npos) {
      std::string old_path = it.value();
      std::string new_path = join_path(root_path, old_path);
      config[key] = new_path;
    }
  }
  return config;
}

}

