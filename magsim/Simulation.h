#ifndef SIMULATION_H
#define SIMULATION_H

#include "SpinLattice.h"

class Simulation {
public:
    virtual ~Simulation() {}
    virtual void setT(real T) = 0;
    virtual real getT() = 0;
    virtual SpinLattice* getLattice() = 0;
    virtual real doStep() = 0;

};

#endif // SIMULATION_H