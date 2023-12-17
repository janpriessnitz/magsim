#include "Heatbath.h"

#include "Constants.h"

Heatbath::Heatbath(SpinLattice *lattice)
  : Simulation(lattice)
{
  Heffs_.resize(lattice->NumSpins());
}

// From this paper:
// Miyatake, Y., Yamamoto, M., Kim, J. J., Toyonaga, M., & Nagai, O. (1986).
// On the implementation of the “heat bath” algorithms for Monte Carlo simulations of classical Heisenberg spin systems.
// Journal of Physics C: Solid State Physics, 19(14), 2539–2546. doi:10.1088/0022-3719/19/14/020
void Heatbath::DoStep() {
  size_t n_spins = lattice_->NumSpins();
  lattice_->ComputeHeffs(Heffs_);

  auto start = std::chrono::high_resolution_clock::now();

  #pragma omp parallel for
  for (size_t i = 0; i < n_spins; ++i) {
    vec3d old_spin = lattice_->spins_[i];
    vec3d Heff = Heffs_[i];
    auto [Heff_x, Heff_y, Heff_z] = Heff;
    real Heff_mag = mag(Heff);

    real beta = 1/(constants::boltzmann*temperature_);
    real R = rnd_uni(0, 1, rng_engs_[omp_get_thread_num()]);
    real Rprime = rnd_uni(0, 1, rng_engs_[omp_get_thread_num()]);

    real phi = Rprime*2*M_PI;
    real HK = Heff_mag*beta;
    // original eq. for costheta from Miyatake paper
    // real costheta = 1/(HK)*log(exp(HK)*(1 - R) + R*exp(-HK));
    // modified eq. for costheta taken from Uppasd which avoids exp(HK) which might be too big
    real costheta = 1 +(1/HK)*log((1 - exp(-2*HK))*R + exp(-2*HK));
    real sintheta = sqrt(1 - costheta*costheta);
    real cosphi = cos(phi);
    real sinphi = sin(phi);
    vec3d spin_prime = {cosphi*sintheta, sinphi*sintheta, costheta};
    auto [spin_prime_x, spin_prime_y, spin_prime_z] = spin_prime;
    real Heffcostheta = Heff_z/Heff_mag;
    real Heffsintheta = sqrt(1 - Heffcostheta*Heffcostheta);
    real Heffcosphi, Heffsinphi;

    if (Heffsintheta > 1e-5) {
      // normal behaviour
      Heffcosphi = Heff_x/Heffsintheta/Heff_mag;
      Heffsinphi = Heff_y/Heffsintheta/Heff_mag;
    } else {
      // avoid dividing by zero if vec{Heff} || z-axis;
      Heffcosphi = 1;
      Heffsinphi = 0;
    }

    // if(Heff_mag < 1e-14) {
    //   printf("!!!!! big bad Heff %lg\n", Heff_mag);
    //   printf("Heff cos/sin phi %lg %lg\n", Heffcosphi, Heffsinphi);
    //   printf("Heff cos/sin theta %lg %lg\n", Heffcostheta, Heffsintheta);
    //   printf("cos/sin phi %lg %lg\n", cosphi, sinphi);
    //   printf("cos/sin theta %lg %lg\n", costheta, sintheta);
    //   printf("Heff_mag*beta %lg\n", Heff_mag*beta);
    // }

    // proper Euler angle (ZXZ) rotation with alpha=phi, beta=theta, gamma=0
    vec3d new_spin = {
      Heffcosphi*Heffcostheta*spin_prime_x - Heffsinphi*spin_prime_y + Heffcosphi*Heffsintheta*spin_prime_z,
      Heffsinphi*Heffcostheta*spin_prime_x + Heffcosphi*spin_prime_y + Heffsinphi*Heffsintheta*spin_prime_z,
                -Heffsintheta*spin_prime_x - 0*spin_prime_y          +            Heffcostheta*spin_prime_z,
    };
    // printf("new spin %s\n", to_string(new_spin).c_str());
    lattice_->spins_[i] = new_spin;
  }
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Integrator);

  lattice_->SampleAverages();
}


