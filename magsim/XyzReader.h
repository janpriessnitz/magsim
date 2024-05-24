
#ifndef SIMBL_XYZREADER_H
#define SIMBL_XYZREADER_H

#include <string>

#include "TupleReader.h"
#include "Util.h"

typedef std::vector<std::tuple<std::string, vec3d>> XyzData;

class XyzReader : TupleReader {
public:
  explicit XyzReader(const std::string &filename);

  XyzData GetData() const;

private:
  void ReadXyzData();

  XyzData xyz_data_;
};


#endif//SIMBL_XYZREADER_H
