
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

//$Log: PeriodicBCXML.cpp,v $
//Revision 1.6  2004/09/10 08:05:24  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.5  2004/07/29 14:25:44  jbarbier
