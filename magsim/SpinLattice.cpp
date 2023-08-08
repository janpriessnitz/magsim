#include "SpinLattice.h"
#include "Constants.h"

#include <cstdio>
#include <unordered_map>
#include <chrono>

SpinLattice::SpinLattice(Timer & timer)
  : timer_(timer)
{
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

  Timer t;
  SpinLattice lat(t);
  lat.anisotropy_ = 0.020;
  lat.exchange_.resize(2*N*N*N);
  lat.spins_.resize(2*N*N*N);


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
  auto start = std::chrono::high_resolution_clock::now();
  FILE *fp = fopen(fname.c_str(), "w");
  for (size_t i = 0; i < spins_.size(); ++i) {
    auto [x, y, z] = spins_[i];
    fprintf(fp, "%lf %lf %lf\n", x, y, z);
  }
  fclose(fp);
  auto stop = std::chrono::high_resolution_clock::now();
  timer_.AddTime(stop - start, Timer::Section::Dump);
}

void SpinLattice::DumpPositions(const std::string &fname) const {
  FILE *fp = fopen(fname.c_str(), "w");
  for (size_t i = 0; i < positions_.size(); ++i) {
    auto [x, y, z] = positions_[i];
    fprintf(fp, "%lf %lf %lf\n", x, y, z);
  }
  fclose(fp);
}

void SpinLattice::DumpProfile(const std::string &fname, char direction) const {
  auto start = std::chrono::high_resolution_clock::now();

  std::unordered_map<real, real> sums;
  std::unordered_map<real, size_t> counts;
  for (size_t ind = 0; ind < spins_.size(); ++ind) {
    auto [sx, sy, sz] = spins_[ind];
    auto [x, y, z] = positions_[ind];
    if (direction == 'x') {
      sums[x] += sz;
      counts[x]++;
    } else if (direction == 'z') {
      sums[z] += sz;
      counts[z]++;
    } else {
      fprintf(stderr, "DumpProfile: unknown direction %c", direction);
      exit(1);
    }
  }

  FILE *fp = fopen(fname.c_str(), "w");
  for(const auto & s : sums) {
    fprintf(fp, "%lf %lf\n", s.first, s.second/counts[s.first]);
  }
  fclose(fp);
  auto stop = std::chrono::high_resolution_clock::now();
  timer_.AddTime(stop - start, Timer::Section::Dump);
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

void SpinLattice::ComputeHeffs(const std::vector<vec3d> spins, std::vector<vec3d> & Heffs) const {
  auto start = std::chrono::high_resolution_clock::now();

  Heffs.resize(spins.size());

  this->ComputeAnis(Heffs, spins);
  this->ComputeExch(Heffs, spins);
  auto end = std::chrono::high_resolution_clock::now();
  timer_.AddTime(end - start, Timer::Section::Heff);
}

void SpinLattice::ComputeHeffs(std::vector<vec3d> & Heffs) const {
  return ComputeHeffs(spins_, Heffs);
}

// TODO: anisotropy in general direction
// Ham: E = K*(S.A)^2
// Heff = 2*K*(S.A).A
void SpinLattice::ComputeAnis(std::vector<vec3d> & Heffs, const std::vector<vec3d> spins) const {
  #pragma omp parallel for simd
  for (size_t i = 0; i < spins_.size(); ++i) {
    real zmag = std::get<2>(spins_[i]);
    Heffs[i] = {0, 0, 2*anisotropy_*zmag};
  }
}

void SpinLattice::ComputeExch(std::vector<vec3d> & Heffs, const std::vector<vec3d> spins) const {
  #pragma omp parallel for simd
  for (size_t i = 0; i < spins_.size(); ++i) {
    vec3d J_field = {0, 0, 0};
    for (size_t j = 0; j < exchange_[i].size(); ++j) {
      auto [spin_ind, J] = exchange_[i][j];
      J_field = J_field - J*spins_[spin_ind];
    }
    Heffs[i] = Heffs[i] + J_field;
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
  exch_en /= spins_.size()*2;  // double counting
  printf("anis: %lg mRy, exch: %lg mRy\n", anis_en/constants::Ry*1000, exch_en/constants::Ry*1000);
}

vec3d SpinLattice::AvgM() const {
  vec3d res = {};
  for (size_t i = 0; i < spins_.size(); ++i) {
    res = res + spins_[i];
  }
  res = (1.0/spins_.size())*res;
  return res;
}
