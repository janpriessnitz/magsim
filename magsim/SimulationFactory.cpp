#include "SimulationFactory.h"

#include "Heatbath.h"
#include "Metropolis.h"
#include "SpinDynamics.h"

std::unique_ptr<Simulation> ConstructSimulation(const Config & config, SpinLattice *lattice) {
  char mode = config.Get("mode").get<std::string>()[0];
  if (mode == 'S') {
    auto sim = std::make_unique<SpinDynamics>(lattice);
    sim->alpha_ = config.Get("damping");
    sim->timestep_ = config.Get("timestep");
    sim->temperature_ = config.Get("temperature");
    return sim;
  } else if (mode == 'M') {
    auto sim = std::make_unique<Metropolis>(lattice);
    sim->temperature_ = config.Get("temperature");
    return sim;
  } else if (mode == 'H') {
    auto sim = std::make_unique<Heatbath>(lattice);
    sim->temperature_ = config.Get("temperature");
    return sim;
  } else {
    fprintf(stderr, "Unknown simulation mode: %c\n", mode);
    exit(1);
  }
}
