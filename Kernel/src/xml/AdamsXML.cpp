#include "AdamsXML.h"
using namespace std;

AdamsXML::AdamsXML() : OneStepIntegratorXML()
{}

AdamsXML::AdamsXML(xmlNode * AdamsNode,  map<int, bool> definedDSNumbersMap)
  : OneStepIntegratorXML(AdamsNode, definedDSNumbersMap)
{}

