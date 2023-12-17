#include "SimulationFactory.h"

#include "Heatbath.h"
#include "Metropolis.h"
#include "SpinDynamics.h"

std::unique_ptr<Simulation> ConstructSimulation(const MapReader & config, SpinLattice *lattice) {
  char mode = config.GetChar("mode");
  if (mode == 'S') {
    auto sim = std::make_unique<SpinDynamics>(lattice);
    sim->alpha_ = config.GetDouble("damping");
    sim->timestep_ = config.GetDouble("timestep");
    sim->temperature_ = config.GetDouble("temperature");
    return sim;
  } else if (mode == 'M') {
    auto sim = std::make_unique<Metropolis>(lattice);
    sim->temperature_ = config.GetDouble("temperature");
    return sim;
  } else if (mode == 'H') {
    auto sim = std::make_unique<Heatbath>(lattice);
    sim->temperature_ = config.GetDouble("temperature");
    return sim;
  } else {
    fprintf(stderr, "Unknown simulation mode: %c\n", mode);
    exit(1);
  }
}
