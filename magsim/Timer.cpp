#include "Timer.h"

Timer::Timer() {
  for (size_t i = 0; i < n_section; ++i) {
    durations_ns_[i] = 0;
  }
}

void Timer::AddTime(std::chrono::nanoseconds duration_ns, Section section) {
  durations_ns_[static_cast<size_t>(section)] += duration_ns.count();
}

void Timer::PrintStatistics() const {
  printf("Timer statistics:\n");
  printf("Heff: %lf s\n", durations_ns_[static_cast<size_t>(Section::Heff)]/1e9);
  printf("Integrator: %lf s\n", durations_ns_[static_cast<size_t>(Section::Integrator)]/1e9);
  printf("Dump: %lf s\n", durations_ns_[static_cast<size_t>(Section::Dump)]/1e9);
  printf("Temperature: %lf s\n", durations_ns_[static_cast<size_t>(Section::Temperature)]/1e9);
  printf("GenExchange: %lf s\n", durations_ns_[static_cast<size_t>(Section::GenExchange)]/1e9);
  printf("Averages: %lf s\n", durations_ns_[static_cast<size_t>(Section::Averages)]/1e9);
}

Timer global_timer;