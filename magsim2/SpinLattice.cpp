#include "SpinLattice.h"

#include <cstdio>


SpinLattice::SpinLattice() {
}

SpinLattice::~SpinLattice() {
}

size_t N = 2;
static size_t index(int x, int y, int z, int t) {
  x = (x + N)%N;
  y = (y + N)%N;
  z = (z + N)%N;
  return 2*(N*(N*x + y) + z) + t;
}

real J1 = 1;

SpinLattice SpinLattice::GenerateFefcc() {


  SpinLattice lat;
  lat.anisotropy_ = 0.020;
  lat.exchange_.resize(2*N*N*N);
  lat.spins_.resize(2*N*N*N);
  lat.Heffs_.resize(2*N*N*N);


  for (int x = 0; x < 10; ++x) {
    for (int y = 0; y < 10; ++y) {
      for(int z = 0; z < 10; ++z) {
        // type 0
        size_t ind0 = index(x, y, z, 0);
        lat.exchange_[ind0].push_back({index(x, y, z+1, 1), J1});
        lat.exchange_[ind0].push_back({index(x, y, z, 1), J1});
        lat.exchange_[ind0].push_back({index(x, y+1, z+1, 1), J1});
        lat.exchange_[ind0].push_back({index(x, y+1, z, 1), J1});
        lat.exchange_[ind0].push_back({index(x+1, y, z+1, 1), J1});
        lat.exchange_[ind0].push_back({index(x+1, y, z, 1), J1});
        lat.exchange_[ind0].push_back({index(x+1, y+1, z+1, 1), J1});
        lat.exchange_[ind0].push_back({index(x+1, y+1, z, 1), J1});

        // type 1
        size_t ind1 = index(x, y, z, 1);
        lat.exchange_[ind1].push_back({index(x, y, z-1, 0), J1});
        lat.exchange_[ind1].push_back({index(x, y, z, 0), J1});
        lat.exchange_[ind1].push_back({index(x, y-1, z-1, 0), J1});
        lat.exchange_[ind1].push_back({index(x, y-1, z, 0), J1});
        lat.exchange_[ind1].push_back({index(x-1, y, z-1, 0), J1});
        lat.exchange_[ind1].push_back({index(x-1, y, z, 0), J1});
        lat.exchange_[ind1].push_back({index(x-1, y-1, z-1, 0), J1});
        lat.exchange_[ind1].push_back({index(x-1, y-1, z, 0), J1});

        // initial config
        lat.spins_[ind0] = {1, 0, 1};
        lat.spins_[ind1] = {1, 0, 1};

        lat.spins_[ind0] = (1/mag(lat.spins_[ind0]))*lat.spins_[ind0];
        lat.spins_[ind1] = (1/mag(lat.spins_[ind1]))*lat.spins_[ind1];
      }
    }
  }

  return lat;
}


void SpinLattice::DumpLattice(const std::string &fname) const {
  FILE *fp = fopen(fname.c_str(), "w");
  for (size_t i = 0; i < spins_.size(); ++i) {
    auto [x, y, z] = spins_[i];
    fprintf(fp, "%lf %lf %lf %lf\n", x, y, z, mag(spins_[i]));
  }
  fclose(fp);
}

void SpinLattice::DumpHeffs(const std::string &fname) const {
  FILE *fp = fopen(fname.c_str(), "w");
  for (size_t i = 0; i < Heffs_.size(); ++i) {
    auto [x, y, z] = Heffs_[i];
    fprintf(fp, "%lf %lf %lf %lg\n", x, y, z, mag(Heffs_[i]));
  }
  fclose(fp);
}

void SpinLattice::DumpPositions(const std::string &fname) const {
  FILE *fp = fopen(fname.c_str(), "w");
  for (size_t i = 0; i < positions_.size(); ++i) {
    auto [x, y, z] = positions_[i];
    fprintf(fp, "%lf %lf %lf\n", x, y, z);
  }
  fclose(fp);
}

void SpinLattice::DumpExchange(const std::string &fname) const {
  FILE *fp = fopen(fname.c_str(), "w");
  for (size_t i = 0; i < exchange_.size(); ++i) {
    for (size_t j = 0; j < exchange_[i].size(); ++j) {
      auto [partner_ind, int_energy] = exchange_[i][j];
      fprintf(fp, "%lu %lu %lg\n", i, partner_ind, int_energy);
    }
  }
  fclose(fp);
}


const std::vector<vec3d> & SpinLattice::ComputeHeffs() {
  this->ComputeAnis();
  this->ComputeExch();

  return this->Heffs_;
}

// TODO: anisotropy in general direction
// Ham: E = K*(SxA)^2
// Heff = 2*K*(SxA)xA
void SpinLattice::ComputeAnis() {
  #pragma omp parallel for simd
  for (size_t i = 0; i < spins_.size(); ++i) {
    real zmag = std::get<2>(spins_[i]);
    Heffs_[i] = {0, 0, 2*anisotropy_*zmag};
  }
}

void SpinLattice::ComputeExch() {
  #pragma omp parallel for simd
  for (size_t i = 0; i < spins_.size(); ++i) {
    vec3d J_field = {0, 0, 0};
    for (size_t j = 0; j < exchange_[i].size(); ++j) {
      auto [spin_ind, J] = exchange_[i][j];
      J_field = J_field - J*spins_[spin_ind];
    }
    Heffs_[i] = Heffs_[i] + J_field;
  }
}

void SpinLattice::PrintEnergy() const {
  real anis_en = 0;
  for (size_t i = 0; i < spins_.size(); ++i) {
    real zmag = std::get<2>(spins_[i]);
    anis_en += anisotropy_*zmag*zmag;
  }
  anis_en /= spins_.size();

  real exch_en = 0;
  for (size_t i = 0; i < spins_.size(); ++i) {
    vec3d J_field = {0, 0, 0};
    for (size_t j = 0; j < exchange_[i].size(); ++j) {
      auto [spin_ind, J] = exchange_[i][j];
      exch_en -= J*scal_prod(spins_[i], spins_[spin_ind]);
    }
  }
  exch_en /= spins_.size();
  printf("anis: %lg, exch: %lg\n", anis_en, exch_en);
}

vec3d SpinLattice::AvgM() const {
  vec3d res = {};
  for (size_t i = 0; i < spins_.size(); ++i) {
    res = res + spins_[i];
  }
  res = (1.0/spins_.size())*res;
  return res;
}
