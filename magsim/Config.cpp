#include "Config.h"

#include "Constants.h"
#include "TupleReader.h"

#include <fstream>

Config::Config(int argc, char **argv) {
  std::string config_file = "config.json";
  if (argc > 1) {
    config_file = argv[1];
  }
  std::ifstream f(config_file);
  data_ = json::parse(f);

  out_dir_ = "output/";
  if (argc > 2) {
    out_dir_ = argv[2];
  }
  if (!std::filesystem::create_directory(out_dir_)) {
    fprintf(stderr, "failed to create output directory %s\n", out_dir_.c_str());
  }

  if (argc > 3) {
    restart_file_ = argv[3];
  }
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
  TupleReader reader = TupleReader(Get<std::string>("symmetry_file"));
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
  TupleReader reader = TupleReader(data_["exchange_file"]);
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
