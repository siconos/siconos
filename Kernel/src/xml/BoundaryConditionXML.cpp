#include "BoundaryConditionXML.h"

BoundaryConditionXML::BoundaryConditionXML()
{}

BoundaryConditionXML::BoundaryConditionXML(xmlNode * rootBCNode)
{
  this->rootBCNode = rootBCNode;
}

BoundaryConditionXML::~BoundaryConditionXML()
{
}
