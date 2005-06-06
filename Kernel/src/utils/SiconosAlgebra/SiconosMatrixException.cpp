#include "SiconosMatrixException.h"
using namespace std;

SiconosMatrixException::SiconosMatrixException() :
  SiconosException("Siconos Matrix Exception") {}

SiconosMatrixException::SiconosMatrixException(const string& report) :
  SiconosException("Siconos Matrix Exception : " + report) {}

SiconosMatrixException::~SiconosMatrixException() {}

void SiconosMatrixException::selfThrow()
{
  throw SiconosMatrixException();
}

void SiconosMatrixException::selfThrow(const string& report)
{
  throw SiconosMatrixException(report);
}
