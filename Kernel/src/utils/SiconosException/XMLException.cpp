#include "XMLException.h"
#include <iostream>
XMLException::XMLException():
  SiconosException("XML Exception") {}

XMLException::XMLException(string report):
  SiconosException("XML Exception : " + report) {}

XMLException::~XMLException() {}

void XMLException::selfThrow()
{
  throw XMLException();
}

void XMLException::selfThrow(string report)
{
  //std::cout << report ; exit(0);
  throw XMLException(report);
}
