#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <string>
#include "Config.h"

class ConfigReader {
public:
  ConfigReader(const std::string &filename);
  Config ReadConfig();

private:
  annealing_sched_t ReadAnnealingSched(std::string filename);

  std::string filename_;
};

#endif