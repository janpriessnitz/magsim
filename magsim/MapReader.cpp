
#include "MapReader.h"

#include <stdexcept>

MapReader::MapReader(const std::string &filename) : TupleReader(filename) {
  data_map_ = GetPairMap();
}
std::string MapReader::GetString(const std::string &key) const {
  try {
    return data_map_.at(key);
  } catch (std::out_of_range &e) {
    fprintf(stderr, "Missing \"%s\" in file \"%s\"\n", key.c_str(), filename_.c_str());
    exit(1);
  }
}
long long int MapReader::GetInt(const std::string &key) const {
  return std::strtoll(GetString(key).c_str(), nullptr, 10);
}
double MapReader::GetDouble(const std::string &key) const {
  return std::strtod(GetString(key).c_str(), nullptr);
}

char MapReader::GetChar(const std::string &key) const {
  return GetString(key)[0];
}

std::map<std::string, std::string> MapReader::GetPairMap() const {
  std::map<std::string, std::string> result;
  for ( auto &r : data_) {
    if (r.size() < 2) continue;
    result[r[0]] = r[1];
  }
  return result;
}