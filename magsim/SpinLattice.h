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
  SpinLattice(Timer & timer);
  ~SpinLattice();

  static SpinLattice GenerateFefcc();

  void DumpLattice(const std::string &fname) const;
  void DumpPositions(const std::string &fname) const;
  void DumpXYZ(const std::string &fname) const;
  void DumpProfile(const std::string &fname, char direction) const;
  void DumpExchange(const std::string &fname) const;

  void ComputeHeffs(const std::vector<vec3d> spins, std::vector<vec3d> & Heffs) const;
  void ComputeHeffs(std::vector<vec3d> & Heffs) const;

  void PrintEnergy() const;

  std::vector<vec3d> spins_;

  real anisotropy_ = 0;
  std::vector<std::vector<std::tuple<size_t, real>>> exchange_;

  vec3d AvgM() const;

  void ComputeAnis(std::vector<vec3d> & Heffs, const std::vector<vec3d> spins) const;
  void ComputeExch(std::vector<vec3d> & Heffs, const std::vector<vec3d> spins) const;

  std::vector<vec3d> positions_;

  Timer & timer_;
};

#endif // SIMULATION_H
