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


class CoRuCoGenerator : LatticeGenerator {
public:
  CoRuCoGenerator(const MapReader & config, Timer & timer);

  SpinLattice Generate() const;

  std::vector<vec3d> GeneratePositions() const;
  std::vector<std::vector<std::tuple<size_t, real>>> GenerateExchange(
    const PointLookup & point_lookup,
    const std::vector<std::tuple<vec3d, real>> & ints,
    const std::vector<mat3d> & syms) const;
  std::vector<vec3d> GenerateSpins(const std::vector<vec3d> & positions) const;

  std::optional<size_t> GetPoint(const PointLookup & lookup, const vec3d & pos) const;

  std::vector<vec3d> ApplySymmetry(const vec3d & vec, const std::vector<mat3d> & syms) const;

  int nx_, ny_, nz_;
  vec3d base1_ = {1, 0, 0};
  vec3d base2_ = {0.5, 0.8660254037844386, 0};
  vec3d base3_ = {0, 0, 1.632993161855452};
  mat3d base_mat_ = {base1_, base2_, base3_};
  std::vector<vec3d> spin_pos_list = {
    {0, 0, 0},
    (1/3.0)*(base1_ + base2_) + (1/2.0)*base3_
  };
  vec3d area_dims_;

  mat3d base_mat_inv_ = inverse(base_mat_);

  real Co_anis_;  // J
  std::string exchange_fname_;
  std::string symmetry_fname_;

  bool periodic_x_;
  bool periodic_y_;
  bool periodic_z_;

  double interface_exchange_energy_ = 0.3/4/constants::J_eV;

  char spin_direction_;

  Timer & timer_;
};

#endif //CORUCOGENERATOR_H