#ifndef SPIN_LATTICE_H
#define SPIN_LATTICE_H

#include <vector>
#include <string>
#include <cstdio>

#include "Config.h"
#include "Util.h"


class SpinLattice {
public:
  SpinLattice();
  ~SpinLattice();

  static SpinLattice GenerateFefcc();

  void DumpLattice(const std::string &fname) const;
  void DumpHeffs(const std::string &fname) const;
  void DumpPositions(const std::string &fname) const;
  void DumpExchange(const std::string &fname) const;

  const std::vector<vec3d> & ComputeHeffs();

  std::vector<vec3d> spins_;

  real anisotropy_ = 0;
  std::vector<std::vector<std::tuple<size_t, real>>> exchange_;

  vec3d AvgM() const;

  void ComputeAnis();
  void ComputeExch();

  std::vector<vec3d> Heffs_;

  std::vector<vec3d> positions_;
};

#endif // SIMULATION_H