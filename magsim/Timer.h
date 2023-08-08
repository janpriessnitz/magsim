#ifndef TIMER_H
#define TIMER_H

#include <vector>
#include <chrono>

#include "Util.h"

class Timer {
public:
  enum class Section {
    Heff = 0,
    Integrator,
    Dump,
    Temperature,
    GenExchange
  };
  const static size_t n_section = 5;

  Timer();
  void AddTime(std::chrono::nanoseconds duration_ns, Section section);
  void PrintStatistics() const;


  // indices correspond to the Section enum elements
  uint64_t durations_ns_[n_section];
};

#endif
