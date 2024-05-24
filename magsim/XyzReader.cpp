

#include "XyzReader.h"

#include <cstdio>

XyzReader::XyzReader(const std::string &filename) : TupleReader(filename) {
  ReadXyzData();
}

XyzData XyzReader::GetData() const {
  return xyz_data_; 
}

void XyzReader::ReadXyzData() {
  size_t num_atoms = GetInt(0, 0);
  if (NumRows() - 2 < num_atoms) {
    fprintf(stderr, "xyz file invalid: expecting %lu atoms, but the file has only %lu lines", num_atoms, NumRows());
  }

  xyz_data_.resize(num_atoms);
  for (size_t i = 0; i < num_atoms; ++i) {
    xyz_data_[i] = {GetString(i+2, 0), {GetDouble(i+2, 1), GetDouble(i+2, 2), GetDouble(i+2, 3)}};
  }
}