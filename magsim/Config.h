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
  Config(const std::string &config_fname);

  json Get(const std::string & key) const;
  annealing_sched_t GetAnnealingSched() const;
  std::vector<mat3d> GetSymmetries() const;
  std::vector<std::tuple<vec3d, real>> GetExchange() const;

  json data_;
};

#endif // UTIL_H