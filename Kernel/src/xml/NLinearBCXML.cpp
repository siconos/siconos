
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

//$Log: NLinearBCXML.cpp,v $
//Revision 1.7  2004/09/10 08:04:51  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.6  2004/07/29 14:25:43  jbarbier
