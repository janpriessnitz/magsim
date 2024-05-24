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

  long long GetInt(size_t r, size_t c);
  double GetDouble(size_t r, size_t c);
  std::string GetString(size_t r, size_t c);
  size_t NumRows();

  std::vector<std::vector<std::string>> data_;

protected:
  std::string filename_;
};


#endif//SIMBL_TUPLEREADER_H
