#include "ConfigReader.h"

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
  conf.J = reader.GetFloat("J");
  conf.D = reader.GetFloat("D");
  conf.B = {0, 0, reader.GetFloat("Bz")};
  conf.deltaSpin = reader.GetFloat("deltaSpin");
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