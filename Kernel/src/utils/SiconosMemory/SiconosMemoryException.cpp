#include "SiconosMemoryException.h"

SiconosMemoryException::SiconosMemoryException() :
  SiconosException("Siconos Memory Exception  (saved values of previous states of simulation)") {}

SiconosMemoryException::SiconosMemoryException(string report) :
  SiconosException("Siconos Memory Exception (saved values of previous states of simulation): " + report) {}

SiconosMemoryException::~SiconosMemoryException() {}

void SiconosMemoryException::selfThrow()
{
  throw SiconosMemoryException();
}


void SiconosMemoryException::selfThrow(string report)
{
  throw SiconosMemoryException(report);
}
