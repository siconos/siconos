#include "SiconosSharedLibraryException.h"
using namespace std;

SiconosSharedLibraryException::SiconosSharedLibraryException() :
  SiconosException("Shared Library Exception") {}

SiconosSharedLibraryException::SiconosSharedLibraryException(const string & report) :
  SiconosException("Shared Library Exception : " + report) {}

SiconosSharedLibraryException::~SiconosSharedLibraryException() {}

void SiconosSharedLibraryException::selfThrow()
{
  throw SiconosSharedLibraryException();
}


void SiconosSharedLibraryException::selfThrow(const string& report)
{
  throw SiconosSharedLibraryException(report);
}
