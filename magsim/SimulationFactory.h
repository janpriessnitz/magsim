#ifndef SIMULATION_FACTORY_H
#define SIMULATION_FACTORY_H

#include <memory>

#include "MapReader.h"
#include "Metropolis.h"
#include "Simulation.h"
#include "SpinDynamics.h"
#include "SpinLattice.h"

std::unique_ptr<Simulation> ConstructSimulation(const MapReader & config, SpinLattice *lattice);

#endif