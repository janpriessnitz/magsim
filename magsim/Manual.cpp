
#include "Config.h"
#include "SimulationFactory.h"
#include "SpinDynamics.h"
#include "Metropolis.h"
#include "LatticeGenerator.h"
#include "CoRuCoGenerator.h"
#include "XyzReader.h"

#include <cstdio>
#include <string>
#include <chrono>
#include <filesystem>
#include <unordered_map>

#include <iostream>

SpinLattice LoadLattice(const Config &config) {

  std::vector<std::string> species_list;
  species_list = config.Get<std::vector<std::string>>("species");

  std::vector<real> anisotropies = config.Get<std::vector<real>>("anisotropy");
  std::unordered_map<std::string, real> species_anisotropy;

  for (size_t i = 0; i < species_list.size(); ++i) {
    species_anisotropy[species_list[i]] = anisotropies[i]*constants::eV;
  }

  SpinLattice res;
  XyzReader xyz_r("pos.xyz");

  auto xyz_data = xyz_r.GetData();

  // printf("loaded xyz\n");

  res.positions_.resize(xyz_data.size());
  res.anisotropy_.resize(xyz_data.size());
  for (size_t i = 0; i < res.positions_.size(); ++i) {
    res.positions_[i] = std::get<1>(xyz_data[i]);
    res.anisotropy_[i] = species_anisotropy[std::get<0>(xyz_data[i])];
  }

  // printf("loading jpairs\n");

  TupleReader ex_r("Jpairs.txt");
  res.exchange_.resize(res.positions_.size());
  for (size_t i = 0; i < ex_r.NumRows(); ++i) {
    long long ia = ex_r.GetInt(i, 0);
    long long ib = ex_r.GetInt(i, 1);
    double en = ex_r.GetDouble(i, 2);
    res.exchange_[ia].push_back({ib, en*constants::eV});
  }

  double z_boundary = std::get<2>(res.positions_[0]);

  res.spins_.resize(res.positions_.size());
  res.avg_spins_.resize(res.positions_.size());

  bool domain_wall = config.Get<int>("domain_wall") != 0;
  for (size_t i = 0; i < res.spins_.size(); ++i) {
    if (std::get<2>(res.positions_[i]) <= z_boundary) {
      if (domain_wall) {
        res.spins_[i] = {0,0,-1};
      } else {
        res.spins_[i] = {0,0,1};
      }
    } else {
      res.spins_[i] = {0,0,1};
    }
  }

  return res;
}


int main(int argc, char **argv) {
  Config c(argc, argv);

  printf("generating spin lattice\n");
  SpinLattice lat = LoadLattice(c);

  auto sim = ConstructSimulation(c, &lat);

  // printf("dumping positions xyz\n");
  // lat.DumpPositions(out_dir + "/positions.out");
  // lat.DumpXYZ(out_dir + "/positions.xyz");

  // if (reader.GetInt("dump_exchange")) {
  //   lat.DumpExchange(out_dir + "/exchange.out");
  // }

  if (!c.restart_file_.empty()) {
    printf("loading restart file\n");
    lat.LoadLattice(c.restart_file_);
  }

  printf("starting sim\n");

  sim->Run();
  return 0;
}
