#ifndef SPIN_LATTICE_H
#define SPIN_LATTICE_H

#include <vector>
#include <string>
#include <cstdio>

#include "Config.h"
#include "Util.h"
#include "Timer.h"


class SpinLattice {
public:
  static SpinLattice GenerateFefcc();

  void LoadLattice(const std::string &fname);

  size_t NumSpins() const;

  void DumpLattice(const std::string &fname, bool dump_average = false) const;
  void DumpPositions(const std::string &fname) const;
  void DumpXYZ(const std::string &fname) const;
  void DumpProfile(const std::string &fname, char direction, bool dump_average = false) const;
  void DumpExchange(const std::string &fname) const;

  void ComputeHeffs(const std::vector<vec3d> spins, std::vector<vec3d> & Heffs) const;
  void ComputeHeffs(std::vector<vec3d> & Heffs) const;

  void SampleAverages();
  void ResetAverages();

  void PrintEnergy() const;

  std::vector<vec3d> spins_;
  std::vector<vec3d> avg_spins_;
  size_t n_avgs_ = 0;

  std::vector<real> anisotropy_;
  std::vector<std::vector<std::tuple<size_t, real>>> exchange_;

  vec3d AvgM() const;

  void ComputeAnis(std::vector<vec3d> & Heffs, const std::vector<vec3d> spins) const;
  void ComputeExch(std::vector<vec3d> & Heffs, const std::vector<vec3d> spins) const;

  std::vector<vec3d> positions_;
};

#endif // SIMULATION_H
