#include "SiconosMatrix.h"

// Default (private) constructor
SiconosMatrix::SiconosMatrix(bool isblock): isBlockMatrix(isblock)
{
  dim.resize(2);
}

SiconosMatrix::~SiconosMatrix() {}

