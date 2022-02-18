#ifndef HEATBATH_H
#define HEATBATH_H

class HeatBathAlgo {
    Metropolis(const Config &conf);
  ~Metropolis();

  real do_step();
  void equilibrize();  
};

#endif