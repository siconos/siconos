
#include "LsodarXML.h"

LsodarXML::LsodarXML() : OneStepIntegratorXML()
{}

LsodarXML::LsodarXML(xmlNode * LsodarNode,  map<int, bool> definedDSNumbers)
  : OneStepIntegratorXML(LsodarNode, definedDSNumbers)
{}

