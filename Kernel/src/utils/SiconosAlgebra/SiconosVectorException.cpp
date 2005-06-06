#include "SiconosVectorException.h"
using namespace std;

SiconosVectorException::SiconosVectorException() :
  SiconosException("Siconos Vector Exception") {}


SiconosVectorException::SiconosVectorException(const string& report) :
  SiconosException("Siconos Vector Exception : " + report) {}


SiconosVectorException::~SiconosVectorException() {}


void SiconosVectorException::selfThrow()
{
  throw SiconosVectorException();
}


void SiconosVectorException::selfThrow(const string& report)
{
  throw SiconosVectorException(report);
}
