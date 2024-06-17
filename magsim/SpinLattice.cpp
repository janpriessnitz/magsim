#include "SpinLattice.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <unordered_map>

#include "Constants.h"

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
  lat.anisotropy_.resize(2*N*N*N, 0.020);
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

void SpinLattice::LoadLattice(const std::string &fname) {
  FILE *fp = fopen(fname.c_str(), "r");
  for (size_t i = 0; i < spins_.size(); ++i) {
    int n = fscanf(fp, "%lg %lg %lg", &std::get<0>(spins_[i]), &std::get<1>(spins_[i]), &std::get<2>(spins_[i]));
    if (n != 3) {
      fprintf(stderr, "could not load lattice\n");
      exit(1);
    }
  }
  fclose(fp);
}


size_t SpinLattice::NumSpins() const {
  return spins_.size();
}

void SpinLattice::DumpLattice(const std::string &fname, bool dump_average) const {
  auto start = std::chrono::high_resolution_clock::now();
  FILE *fp = fopen(fname.c_str(), "w");
  if (!fp) {
    fprintf(stderr, "could not open dump file %s\n", fname.c_str());
    exit(1);
  }
  for (size_t i = 0; i < spins_.size(); ++i) {
    vec3d spin = dump_average ? avg_spins_[i]/n_avgs_ : spins_[i];
    fprintf(fp, "%s\n", to_string(spin).c_str());
  }
  fclose(fp);
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Dump);
}

void SpinLattice::DumpPositions(const std::string &fname) const {
  FILE *fp = fopen(fname.c_str(), "w");
  if (fp == nullptr) {
    fprintf(stderr, "could not open file for dumping positions: %s\n", fname.c_str());
    return;
  }
  for (size_t i = 0; i < positions_.size(); ++i) {
    fprintf(fp, "%s\n", to_string(positions_[i]).c_str());
    // printf("%s\n", to_string(positions_[i]).c_str());
  }
  fclose(fp);
}

void SpinLattice::DumpXYZ(const std::string &fname) const {
  FILE *fp = fopen(fname.c_str(), "w");
  fprintf(fp, "%lu\ncomment\n", positions_.size());

  for (size_t i = 0; i < positions_.size(); ++i) {
    fprintf(fp, "%s\n", to_string(positions_[i]).c_str());
  }
  fclose(fp);
}


void SpinLattice::DumpProfile(const std::string &fname, char direction, bool dump_average) const {
  auto start = std::chrono::high_resolution_clock::now();

  std::unordered_map<real, std::pair<vec3d, size_t>> sums;
  for (size_t ind = 0; ind < spins_.size(); ++ind) {
    vec3d spin;
    if (dump_average) {
      spin = avg_spins_[ind]/n_avgs_;
    } else {
      spin = spins_[ind];
    }
    auto [x, y, z] = positions_[ind];
    if (direction == 'x') {
      sums[x].first += spin;
      sums[x].second++;
    } else if (direction == 'y') {
      sums[y].first += spin;
      sums[y].second++;
    } else if (direction == 'z') {
      sums[z].first += spin;
      sums[z].second++;
    } else {
      fprintf(stderr, "DumpProfile: unknown direction %c", direction);
      exit(1);
    }
  }
  // order by position
  std::vector<std::pair<real, std::pair<vec3d, size_t>>> sum_list;
  sum_list.insert(sum_list.begin(), sums.begin(), sums.end());
  std::sort(sum_list.begin(), sum_list.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  FILE *fp = fopen(fname.c_str(), "w");
  for(const auto & s : sum_list) {
    fprintf(fp, "%lg %s\n", s.first, to_string(s.second.first/s.second.second).c_str());
  }
  fclose(fp);
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Dump);
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
  global_timer.AddTime(end - start, Timer::Section::Heff);
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
    Heffs[i] = {0, 0, 2*anisotropy_[i]*zmag};
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

void SpinLattice::SampleAverages() {
  auto start = std::chrono::high_resolution_clock::now();
  #pragma omp parallel for simd
  for (size_t i = 0; i < spins_.size(); ++i) {
    avg_spins_[i] += spins_[i];
  }
  ++n_avgs_;
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Averages);

}

void SpinLattice::ResetAverages() {
  auto start = std::chrono::high_resolution_clock::now();
  #pragma omp parallel for simd
  for (size_t i = 0; i < spins_.size(); ++i) {
    avg_spins_[i] = {0, 0, 0};
  }
  n_avgs_ = 0;
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Averages);
}

void SpinLattice::PrintEnergy() const {
  real anis_en = 0;
  for (size_t i = 0; i < spins_.size(); ++i) {
    real zmag = std::get<2>(spins_[i]);
    anis_en += anisotropy_[i]*zmag*zmag;
  }

  real exch_en = 0;
  for (size_t i = 0; i < spins_.size(); ++i) {
    vec3d J_field = {0, 0, 0};
    for (size_t j = 0; j < exchange_[i].size(); ++j) {
      auto [spin_ind, J] = exchange_[i][j];
      exch_en -= J*scal_prod(spins_[i], spins_[spin_ind]);
    }
  }
  exch_en /= 2;  // double counting
  printf("tot  - anis: %lg mRy, exch: %lg mRy\n", anis_en/constants::Ry*1000, exch_en/constants::Ry*1000);
  anis_en /= spins_.size();
  exch_en /= spins_.size();
  printf("p/at - anis: %lg mRy, exch: %lg mRy\n", anis_en/constants::Ry*1000, exch_en/constants::Ry*1000);
}

vec3d SpinLattice::AvgM() const {
  vec3d res = {};
  for (size_t i = 0; i < spins_.size(); ++i) {
    res = res + spins_[i];
  }
  res = (1.0/spins_.size())*res;
  return res;
}
