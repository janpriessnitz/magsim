#include "HeatBath.h"

#include <cmath>

HeatBath::HeatBath(const Config& conf)
  : lattice_(conf)
  , conf_(conf)
{
  // E_history_fp_ = fopen("Ehistory.out", "w");
}

HeatBath::~HeatBath() {
  // fclose(E_history_fp_);
}


void HeatBath::setT(real T) {
  T_ = T;
}

real HeatBath::getT() {
  return T_;
}

SpinLattice* HeatBath::getLattice() {
  return &lattice_;
}

real HeatBath::doStep() {
  int64_t rx = rnd_int(0, lattice_.w_);
  int64_t ry = rnd_int(0, lattice_.h_);
  vec3d oldSpin = lattice_.get(rx, ry);

  vec3d field = lattice_.getLocalField(rx, ry);
  real magfield = mag(field);  // equals E, since mag(spin) = 1
  real beta = 1/T_;
  real R = rnd_uni(0, 1);
  real Rprime = rnd_uni(0, 1);

  real phi = Rprime*2*M_PI;
  real costheta;
  if (magfield*beta > 500) {
    // prevents overflow in exp function
    costheta = -1/(magfield*beta)*(magfield*beta + log(1-R));
  } else {
    costheta = -1/(magfield*beta)*log(exp(magfield*beta)*(1 - R) + R*exp(-magfield*beta));
  }
  // printf("costhetaprime %lf\n", costheta);
  real sintheta = sqrt(1 - costheta*costheta);
  real cosphi = cos(phi);
  real sinphi = sin(phi);

  vec3d spinPrime = {cosphi*sintheta, sinphi*sintheta, costheta};
  // vec3d spinPrime = {0, 0, 1};

  real fieldcostheta = std::get<2>(field)/magfield;
  real fieldsintheta = sqrt(1 - fieldcostheta*fieldcostheta);
  real fieldcosphi = std::get<0>(field)/fieldsintheta/magfield;
  real fieldsinphi = std::get<1>(field)/fieldsintheta/magfield;

  // vec3d newSpin = {
  //   fieldcostheta*std::get<0>(spinPrime) + 0*std::get<1>(spinPrime) + fieldsintheta*fieldcosphi*std::get<2>(spinPrime),
  //   0*std::get<0>(spinPrime) + fieldcostheta*std::get<1>(spinPrime) + fieldsintheta*fieldsinphi*std::get<2>(spinPrime),
  //   -fieldsintheta*fieldcosphi*std::get<0>(spinPrime) - fieldsintheta*fieldsinphi*std::get<1>(spinPrime) + fieldcostheta*std::get<2>(spinPrime),
  // };

  vec3d newSpin = {
    fieldcostheta*fieldcosphi*std::get<0>(spinPrime) - fieldsinphi*std::get<1>(spinPrime) + fieldsintheta*fieldcosphi*std::get<2>(spinPrime),
    fieldsinphi*fieldcostheta*std::get<0>(spinPrime) + fieldcosphi*std::get<1>(spinPrime) + fieldsintheta*fieldsinphi*std::get<2>(spinPrime),
    -fieldsintheta*std::get<0>(spinPrime) - 0*std::get<1>(spinPrime) + fieldcostheta*std::get<2>(spinPrime),
  };


  lattice_.set(rx, ry, newSpin);
  // printf("costheta %lf %lf\n", costheta, scal_prod(newSpin, field)/mag(newSpin)/mag(field));
  // printf("field %lf %lf %lf, spinPrime %lf %lf %lf, spin %lf %lf %lf\n", std::get<0>(field)/magfield, std::get<1>(field)/magfield, std::get<2>(field)/magfield,
  // std::get<0>(spinPrime), std::get<1>(spinPrime), std::get<2>(spinPrime),
  // std::get<0>(newSpin), std::get<1>(newSpin), std::get<2>(newSpin)
  // );

  return scal_prod(newSpin, field) - scal_prod(oldSpin, field);
}
