#include "SiconosException.h"

SiconosException::SiconosException()
{
  this->reportMsg = "Siconos Exception";
}

SiconosException::SiconosException(string report)
{
  this->reportMsg = report;
}

SiconosException::~SiconosException() {}
