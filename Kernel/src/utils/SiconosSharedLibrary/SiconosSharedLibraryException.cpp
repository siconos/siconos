#include "SiconosSharedLibraryException.h"

SiconosSharedLibraryException::SiconosSharedLibraryException() :
  SiconosException("Shared Library Exception") {}

SiconosSharedLibraryException::SiconosSharedLibraryException(string report) :
  SiconosException("Shared Library Exception : " + report) {}

SiconosSharedLibraryException::~SiconosSharedLibraryException() {}

void SiconosSharedLibraryException::selfThrow()
{
  throw SiconosSharedLibraryException();
}


void SiconosSharedLibraryException::selfThrow(string report)
{
  throw SiconosSharedLibraryException(report);
}
