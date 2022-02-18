//
// Created by jp on 09.12.20.
//

#ifndef SIMBL_TUPLEREADER_H
#define SIMBL_TUPLEREADER_H


#include <map>
#include <string>
#include <vector>

class TupleReader {
public:
  TupleReader(const std::string &filename);
  std::map<std::string, std::string> GetPairMap();

  long long GetInt(size_t r, size_t c);
  double GetDouble(size_t r, size_t c);

  std::vector<std::vector<std::string>> data;

protected:
  std::string filename_;
};


#endif//SIMBL_TUPLEREADER_H
