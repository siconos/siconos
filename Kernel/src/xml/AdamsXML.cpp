#include "AdamsXML.h"

AdamsXML::AdamsXML() : OneStepIntegratorXML()
{}

AdamsXML::AdamsXML(xmlNode * AdamsNode,  map<int, bool> definedDSNumbersMap)
  : OneStepIntegratorXML(AdamsNode, definedDSNumbersMap)
{}

