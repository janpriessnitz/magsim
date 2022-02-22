#include "SpinLattice.h"

#include <cstdio>


SpinLattice::SpinLattice(const Config &conf)
  : w_(conf.lattice_w)
  , h_(conf.lattice_h)
  , J_(conf.J)
  , D_(conf.D)
  , B_(conf.B)
{
  lattice_.resize(w_*h_);
  init_random();
}

SpinLattice::~SpinLattice() {
}

vec3d SpinLattice::get(int64_t x, int64_t y) {
  return lattice_[index(x, y)];
}

void SpinLattice::set(int64_t x, int64_t y, vec3d spin) {
  lattice_[index(x, y)] = spin;
}

uint64_t SpinLattice::index(int64_t x, int64_t y) {
  // periodic coordinates
  x = (x + w_) % w_;
  y = (y + h_) % h_;
  return x + y*w_;
}

real SpinLattice::getEnergy() {
  real E = 0;
  for (int64_t x = 0; x < w_; ++x) {
    for (int64_t y = 0; y < h_; ++y) {
      E += -J_*scal_prod(get(x, y), get(x+1, y));
      E += -J_*scal_prod(get(x, y), get(x, y+1));
      E += -J_*scal_prod(get(x, y), get(x-1, y));
      E += -J_*scal_prod(get(x, y), get(x, y-1));
      E += -scal_prod({0, D_, 0}, vec_prod(get(x, y), get(x+1, y)));
      E += -scal_prod({0, -D_, 0}, vec_prod(get(x, y), get(x-1, y)));
      E += -scal_prod({-D_, 0, 0}, vec_prod(get(x, y), get(x, y+1)));
      E += -scal_prod({D_, 0, 0}, vec_prod(get(x, y), get(x, y-1)));
      E += -scal_prod(B_, get(x, y));
    }
  }
  return E;
}

real SpinLattice::getEnergyDelta(int64_t x, int64_t y, vec3d newspin) {
  real deltaE = 0;
  deltaE += -J_*scal_prod(newspin, get(x+1, y));
  deltaE += -J_*scal_prod(newspin, get(x, y+1));
  deltaE += -J_*scal_prod(newspin, get(x-1, y));
  deltaE += -J_*scal_prod(newspin, get(x, y-1));
  deltaE += -scal_prod({0, D_, 0}, vec_prod(newspin, get(x+1, y)));
  deltaE += -scal_prod({0, -D_, 0}, vec_prod(newspin, get(x-1, y)));
  deltaE += -scal_prod({-D_, 0, 0}, vec_prod(newspin, get(x, y+1)));
  deltaE += -scal_prod({D_, 0, 0}, vec_prod(newspin, get(x, y-1)));
  deltaE += -scal_prod(B_, newspin);

  deltaE -= -J_*scal_prod(get(x, y), get(x+1, y));
  deltaE -= -J_*scal_prod(get(x, y), get(x, y+1));
  deltaE -= -J_*scal_prod(get(x, y), get(x-1, y));
  deltaE -= -J_*scal_prod(get(x, y), get(x, y-1));
  deltaE -= -scal_prod({0, D_, 0}, vec_prod(get(x, y), get(x+1, y)));
  deltaE -= -scal_prod({0, -D_, 0}, vec_prod(get(x, y), get(x-1, y)));
  deltaE -= -scal_prod({-D_, 0, 0}, vec_prod(get(x, y), get(x, y+1)));
  deltaE -= -scal_prod({D_, 0, 0}, vec_prod(get(x, y), get(x, y-1)));
  deltaE -= -scal_prod(B_, get(x, y));
  return deltaE;
}

void SpinLattice::dump(std::string filename) {
  FILE *fp = fopen(filename.c_str(), "w");
  fprintf(fp, "%ld %ld %le\n", w_, h_, getEnergy());
  for (uint64_t i; i < w_*h_; ++i) {
    fprintf(fp, "%lf %lf %lf\n", std::get<0>(lattice_[i]), std::get<1>(lattice_[i]), std::get<2>(lattice_[i]));
  }
  fclose(fp);
}


void SpinLattice::init_random() {
  for (uint64_t i = 0; i < w_*h_; ++i) {
    lattice_[i] = rnd_vec();
  }
}
