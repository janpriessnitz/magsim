#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>
#include <tuple>
#include <vector>

#include "Util.h"

#include "json.hpp"
using json = nlohmann::json;

typedef std::tuple<real, uint64_t, uint64_t> annealing_step_t;
typedef std::vector<annealing_step_t> annealing_sched_t;

class Config {
public:
  Config(int argc, char **argv);

  template<typename T>
  T Get(const std::string & key) const {
    try {
      return data_.at(key).get<T>();
    } catch (...) {
      fprintf(stderr, "failed to get config key \"%s\"\n", key.c_str());
      throw;
    }
  }

  template<typename T>
  T Get(const std::string & key, T default_value) const {
    try {
      return data_.at(key).get<T>();
    } catch (...) {
      return default_value;
    }
  }

  annealing_sched_t GetAnnealingSched() const;
  std::vector<mat3d> GetSymmetries() const;
  std::vector<std::tuple<vec3d, real>> GetExchange() const;

  json data_;
  std::string out_dir_;
  std::string restart_file_;
};

#endif // UTIL_H