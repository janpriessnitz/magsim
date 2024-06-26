#ifndef CORUCOGENERATOR_H
#define CORUCOGENERATOR_H

#include <vector>
#include <unordered_map>
#include <optional>

#include "Constants.h"
#include "LatticeGenerator.h"
#include "SpinLattice.h"
#include "Util.h"
#include "MapReader.h"
#include "Config.h"

class CoRuCoGenerator : LatticeGenerator {
public:
  CoRuCoGenerator(const Config & config);

  SpinLattice Generate() const;

  std::vector<vec3d> GeneratePositions() const;
  std::vector<std::vector<std::tuple<size_t, real>>> GenerateExchange(
    const PointLookup & point_lookup,
    const std::vector<std::tuple<vec3d, real>> & ints,
    const std::vector<mat3d> & syms) const;
  std::vector<vec3d> GenerateSpins(const std::vector<vec3d> & positions) const;

  vec3d area_dims_;

  real Co_anis_;  // J

  std::vector<mat3d> syms_;
  std::vector<std::tuple<vec3d, real>> exchange_ints_;

  double interface_exchange_energy_ = 0.3/4/constants::J_eV;

  char spin_direction_;
};

#endif //CORUCOGENERATOR_H