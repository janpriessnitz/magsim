#ifndef BULKCOGENERATOR_H
#define BULKCOGENERATOR_H

#include <vector>
#include <unordered_map>
#include <optional>

#include "Constants.h"
#include "LatticeGenerator.h"
#include "SpinLattice.h"
#include "Util.h"
#include "Config.h"


class BulkCoGenerator : LatticeGenerator {
public:
  BulkCoGenerator(const Config & config);

  SpinLattice Generate() const;

  std::vector<vec3d> GeneratePositions() const;
  std::vector<std::vector<std::tuple<size_t, real>>> GenerateExchange(
    const PointLookup & point_lookup,
    const std::vector<std::tuple<vec3d, real>> & ints,
    const std::vector<mat3d> & syms) const;
  std::vector<vec3d> GenerateSpins(const std::vector<vec3d> & positions) const;

  std::pair<std::optional<size_t>, int> GetPoint(const PointLookup & lookup, const vec3d & pos) const;

  std::vector<vec3d> ApplySymmetry(const vec3d & vec, const std::vector<mat3d> & syms) const;

  vec3d base1_ = {1, 0, 0};
  vec3d base2_ = {0.5, 0.8660254037844386, 0};
  vec3d base3_ = {0, 0, 1.632993161855452};
  std::vector<vec3d> spin_pos_list = {
    {0, 0, 0},
    (1/3.0)*(base1_ + base2_) + (1/2.0)*base3_
  };

  vec3d area_dims_;
  real Co_anis_;  // J

  std::vector<mat3d> syms_;
  std::vector<std::tuple<vec3d, real>> exchange_ints_;

  double interface_exchange_energy_ = 0.3/4/constants::J_eV;

  char spin_direction_;
};

#endif //CORUCOGENERATOR_H