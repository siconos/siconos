#include "SiconosMatrixException.h"

SiconosMatrixException::SiconosMatrixException() :
  SiconosException("Siconos Matrix Exception") {}

SiconosMatrixException::SiconosMatrixException(string report) :
  SiconosException("Siconos Matrix Exception : " + report) {}

SiconosMatrixException::~SiconosMatrixException() {}

void SiconosMatrixException::selfThrow()
{
  throw SiconosMatrixException();
}

void SiconosMatrixException::selfThrow(string report)
{
  throw SiconosMatrixException(report);
}
