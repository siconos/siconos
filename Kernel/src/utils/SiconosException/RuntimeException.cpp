#include "RuntimeException.h"

RuntimeException::RuntimeException():
  SiconosException("Runtime Exception") {}

RuntimeException::RuntimeException(string report):
  SiconosException("Runtime Exception : " + report) {}

RuntimeException::~RuntimeException() {}

void RuntimeException::selfThrow()
{
  throw RuntimeException();
}

void RuntimeException::selfThrow(string report)
{
  throw RuntimeException(report);
}
