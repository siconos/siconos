#include "SiconosMemoryException.h"
using namespace std;

SiconosMemoryException::SiconosMemoryException() :
  SiconosException("Siconos Memory Exception  (saved values of previous states of simulation)") {}

SiconosMemoryException::SiconosMemoryException(const string& report) :
  SiconosException("Siconos Memory Exception (saved values of previous states of simulation): " + report) {}

SiconosMemoryException::~SiconosMemoryException() {}

void SiconosMemoryException::selfThrow()
{
  throw SiconosMemoryException();
}


void SiconosMemoryException::selfThrow(const string& report)
{
  throw SiconosMemoryException(report);
}
