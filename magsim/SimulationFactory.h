#ifndef SIMULATION_FACTORY_H
#define SIMULATION_FACTORY_H

#include <memory>

#include "MapReader.h"
#include "Simulation.h"
#include "SpinLattice.h"

std::unique_ptr<Simulation> ConstructSimulation(const Config & config, SpinLattice *lattice);

#endif