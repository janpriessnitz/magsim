#include "SimulationFactory.h"

#include "Heatbath.h"
#include "Metropolis.h"
#include "SpinDynamics.h"

std::unique_ptr<Simulation> ConstructSimulation(const Config & config, SpinLattice *lattice) {
  char mode = config.Get<std::string>("mode")[0];
  if (mode == 'S') {
    auto sim = std::make_unique<SpinDynamics>(config, lattice);
    return sim;
  } else if (mode == 'M') {
    auto sim = std::make_unique<Metropolis>(config, lattice);
    return sim;
  } else if (mode == 'H') {
    auto sim = std::make_unique<Heatbath>(config, lattice);
    return sim;
  } else {
    fprintf(stderr, "Unknown simulation mode: %c\n", mode);
    exit(1);
  }
}
