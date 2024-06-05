#include "Config.h"

#include "Constants.h"
#include "TupleReader.h"

#include <fstream>

Config::Config(const std::string &config_fname) {
  std::ifstream f(config_fname);
  data_ = json::parse(f);
}

annealing_sched_t Config::GetAnnealingSched() const {
  TupleReader reader(data_["annealing_sched_fname"]);
  annealing_sched_t sched;
  for (int i = 0; i < reader.NumRows(); ++i) {
    annealing_step_t step = {reader.GetDouble(i, 0), reader.GetInt(i, 1), reader.GetInt(i, 2)};
    sched.push_back(step);
  }
  return sched;
}

std::vector<mat3d> Config::GetSymmetries() const {
  TupleReader reader = TupleReader(data_["symmetries_fname"]);
  std::vector<mat3d> res;
  int n_syms = reader.GetInt(0, 0);
  for (int i_sym = 0; i_sym < n_syms; ++i_sym) {
    int sr = i_sym*3 + 1;
    res.push_back({
      {reader.GetDouble(sr, 0), reader.GetDouble(sr, 1), reader.GetDouble(sr, 2)},
      {reader.GetDouble(sr+1, 0), reader.GetDouble(sr+1, 1), reader.GetDouble(sr+1, 2)},
      {reader.GetDouble(sr+2, 0), reader.GetDouble(sr+2, 1), reader.GetDouble(sr+2, 2)},
    });
  }
  return res;
}

std::vector<std::tuple<vec3d, real>> Config::GetExchange() const {
  TupleReader reader = TupleReader(data_["exchange_fname"]);
  std::vector<std::tuple<vec3d, real>> res;
  for (int i = 0; i < reader.NumRows(); ++i) {
    vec3d vec = {reader.GetDouble(i, 0), reader.GetDouble(i, 1), reader.GetDouble(i, 2)};
    // maptype = 2 in UppASD
    // vec = vec*base_mat_;
    real energy = reader.GetDouble(i, 3)*constants::Ry;
    energy *= 2;  // consistent with UppASD, TODO: Remove
    res.push_back({vec, energy});
  }
  return res;
}
