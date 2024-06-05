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

  annealing_sched_t Config::GetAnnealingSched() const;
  std::vector<mat3d> Config::GetSymmetries() const;
  std::vector<std::tuple<vec3d, real>> Config::GetExchange() const;

  json data_;
};

#endif // UTIL_H