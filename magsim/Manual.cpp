
#include "ConfigReader.h"
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


SpinLattice LoadLattice(const MapReader &config_reader) {
  
  std::string species_str = config_reader.GetString("species");
  std::vector<std::string> species_list;
  size_t pos = 0;
  while ((pos = species_str.find(",")) != std::string::npos) {
      std::string specie = species_str.substr(0, pos);
      species_list.push_back(specie);
      species_str.erase(0, pos + 1);
  }
  species_list.push_back(species_str);

  std::unordered_map<std::string, real> species_anisotropy;
  std::string anisotropy_str = config_reader.GetString("anisotropy");
  size_t i = 0;
  for (; (pos = anisotropy_str.find(",")) != std::string::npos; ++i) {
      std::string tok = anisotropy_str.substr(0, pos);
      real anis = strtod(tok.c_str(), nullptr);
      species_anisotropy[species_list[i]] = anis*constants::eV;
      anisotropy_str.erase(0, pos + 1);
  }
  real anis = strtod(anisotropy_str.c_str(), nullptr);
  species_anisotropy[species_list[i]] = anis*constants::eV;
      
  // // species_anisotropy["Fe"] = -2.4e-6*constants::eV;  // TODO: find better source for bcc Fe anisotropy
  // species_anisotropy["Fe"] = -2.4e-6*constants::eV;  // TODO: find better source for bcc Fe anisotropy
  // // species_anisotropy["Fe"] = 0;  // TODO: find better source for bcc Fe anisotropy
  // species_anisotropy["Pt"] = 0;

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

  bool domain_wall = config_reader.GetInt("domain_wall") != 0;
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
  std::string out_dir = "output/";
  if (argc > 2) {
    out_dir = argv[2];
  }
  std::string restart_fname = "";
  if (argc > 3) {
    restart_fname = argv[3];
  }

  if (!std::filesystem::create_directory(out_dir)) {
    fprintf(stderr, "failed to create output directory %s\n", out_dir.c_str());
  }

  MapReader reader(argv[1]);

  printf("generating spin lattice\n");
  SpinLattice lat = LoadLattice(reader);

  auto sim = ConstructSimulation(reader, &lat);

  // printf("dumping positions xyz\n");
  // lat.DumpPositions(out_dir + "/positions.out");
  // lat.DumpXYZ(out_dir + "/positions.xyz");

  // if (reader.GetInt("dump_exchange")) {
  //   lat.DumpExchange(out_dir + "/exchange.out");
  // }

  // if (restart_fname.length()) {
  //   printf("loading restart file\n");
  //   lat.LoadLattice(restart_fname);
  // }

  printf("starting sim\n");

  int64_t num_step = reader.GetInt("num_step");
  int64_t num_substep = reader.GetInt("num_substep");

  auto start = std::chrono::high_resolution_clock::now();
  for (int j = 0; j < num_step; ++j) {
    for (int i = 0; i < num_substep; ++i) {
      sim->DoStep();
    }
    sim->lattice_->PrintEnergy();
    auto avgm = sim->lattice_->AvgM();
    printf("%s %lf\n", to_string(avgm).c_str(), mag(avgm));
    bool dump_avgs = true;
    sim->lattice_->DumpLattice(out_dir + "/lattice.out" + std::to_string(j), dump_avgs);
    sim->lattice_->DumpProfile(out_dir + "/profile.out" + std::to_string(j), reader.GetChar("domain_wall_direction"), dump_avgs);
    sim->lattice_->ResetAverages();
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  printf("took %lu ms\n", duration.count());
  global_timer.PrintStatistics();
  return 0;
}
