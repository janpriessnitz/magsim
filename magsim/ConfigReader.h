#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <string>
#include "Config.h"

class ConfigReader {
public:
  ConfigReader(const std::string &filename);
  Config ReadConfig();

  static std::vector<mat3d> ReadSymmetries(const std::string & fname);
  static std::vector<std::tuple<vec3d, real>> ReadExchange(const std::string & fname);
private:
  annealing_sched_t ReadAnnealingSched(std::string filename);

  std::string filename_;
};

#endif