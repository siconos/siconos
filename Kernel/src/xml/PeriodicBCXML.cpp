
#include "PeriodicBCXML.h"

PeriodicBCXML::PeriodicBCXML()
{
}

PeriodicBCXML::PeriodicBCXML(xmlNode * PeriodicBCNode)
  : BoundaryConditionXML(PeriodicBCNode)
{
}

PeriodicBCXML::~PeriodicBCXML()
{
}

void PeriodicBCXML::updateBoundaryConditionXML(xmlNode* node) //, BoundaryCondition* bc)
{
  IN("PeriodicBCXML::updateBoundaryConditionXML\n");
  this->rootBCNode = node;
  OUT("PeriodicBCXML::updateBoundaryConditionXML\n");
}

