#include "ConfigReader.h"

#include "MapReader.h"

ConfigReader::ConfigReader(const std::string &filename)
  : filename_(filename)
{}

Config ConfigReader::ReadConfig() {
  MapReader reader(filename_);
  Config conf;
  conf.lattice_w = reader.GetInt("lattice_w");
  conf.lattice_h = reader.GetInt("lattice_h");
  conf.J = reader.GetFloat("J");
  conf.D = reader.GetFloat("D");
  conf.B = {0, 0, reader.GetFloat("Bz")};
  conf.T_init = reader.GetFloat("T_init");
  conf.T_step_ratio = reader.GetFloat("T_step_ratio");
  conf.T_steps = reader.GetInt("T_steps");
  conf.deltaSpin = reader.GetFloat("deltaSpin");
  conf.metropolis_reporting_macrostep = reader.GetInt("metropolis_reporting_macrostep");
  conf.metropolis_equilibrium_macrosteps = reader.GetInt("metropolis_equilibrium_macrosteps");
  return conf;
}