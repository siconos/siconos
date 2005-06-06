#include "XMLException.h"
using namespace std;

XMLException::XMLException():
  SiconosException("XML Exception") {}

XMLException::XMLException(const string& report):
  SiconosException("XML Exception : " + report) {}

XMLException::~XMLException() {}

void XMLException::selfThrow()
{
  throw XMLException();
}

void XMLException::selfThrow(const string& report)
{
  throw XMLException(report);
}
