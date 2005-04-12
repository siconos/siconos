
#include "NLinearBCXML.h"

NLinearBCXML::NLinearBCXML()
{}

NLinearBCXML::NLinearBCXML(xmlNode * NLinearBCNode)
  : BoundaryConditionXML(NLinearBCNode)
{}

NLinearBCXML::~NLinearBCXML()
{}

void NLinearBCXML::updateBoundaryConditionXML(xmlNode* node) //, BoundaryCondition* bc)
{
  IN("NLinearBCXML::updateBoundaryConditionXML\n");
  this->rootBCNode = node;
  OUT("NLinearBCXML::updateBoundaryConditionXML\n");
}

