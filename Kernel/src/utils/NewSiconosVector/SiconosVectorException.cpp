#include "SiconosVectorException.h"


SiconosVectorException::SiconosVectorException() :
  SiconosException("Siconos Vector Exception") {}


SiconosVectorException::SiconosVectorException(string report) :
  SiconosException("Siconos Vector Exception : " + report) {}


SiconosVectorException::~SiconosVectorException() {}


void SiconosVectorException::selfThrow()
{
  throw SiconosVectorException();
}


void SiconosVectorException::selfThrow(string report)
{
  throw SiconosVectorException(report);
}
