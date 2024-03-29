#ifndef SIMBL_MAPREADER_H
#define SIMBL_MAPREADER_H


#include "TupleReader.h"

class MapReader : TupleReader {
public:
  explicit MapReader(const std::string &filename);

  std::string GetString(const std::string &key) const;
  long long int GetInt(const std::string &key) const;
  double GetDouble(const std::string &key) const;
  char GetChar(const std::string &key) const;

  std::map<std::string, std::string> GetPairMap() const;

  std::map<std::string, std::string> data_map_;
};


#endif//SIMBL_MAPREADER_H
