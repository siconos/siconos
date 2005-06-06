#include "RuntimeException.h"
using namespace std;

RuntimeException::RuntimeException():
  SiconosException("Runtime Exception") {}

RuntimeException::RuntimeException(const string& report):
  SiconosException("Runtime Exception : " + report) {}

RuntimeException::~RuntimeException() {}

void RuntimeException::selfThrow()
{
  throw RuntimeException();
}

void RuntimeException::selfThrow(const string& report)
{
  throw RuntimeException(report);
}
