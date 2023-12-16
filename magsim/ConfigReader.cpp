#include "ConfigReader.h"

#include "Constants.h"
#include "MapReader.h"


ConfigReader::ConfigReader(const std::string &filename)
  : filename_(filename)
{}

Config ConfigReader::ReadConfig() {
  MapReader reader(filename_);
  Config conf;
  conf.method = reader.GetChar("method");
  conf.lattice_w = reader.GetInt("lattice_w");
  conf.lattice_h = reader.GetInt("lattice_h");
  conf.J = reader.GetDouble("J");
  conf.D = reader.GetDouble("D");
  conf.B = {0, 0, reader.GetDouble("Bz")};
  conf.deltaSpin = reader.GetDouble("deltaSpin");
  conf.annealing_sched = ReadAnnealingSched(reader.GetString("annealing_sched_file"));
  conf.lattice_dump_file = reader.GetString("lattice_dump_file");

  return conf;
}

annealing_sched_t ConfigReader::ReadAnnealingSched(std::string filename) {
  TupleReader reader(filename);
  annealing_sched_t sched;
  for (int i = 0; i < reader.NumRows(); ++i) {
    annealing_step_t step = {reader.GetDouble(i, 0), reader.GetInt(i, 1), reader.GetInt(i, 2)};
    sched.push_back(step);
  }
  return sched;
}

std::vector<mat3d> ConfigReader::ReadSymmetries(const std::string & fname) {
  TupleReader reader = TupleReader(fname);
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

std::vector<std::tuple<vec3d, real>> ConfigReader::ReadExchange(const std::string & fname) {
  TupleReader reader = TupleReader(fname);
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
