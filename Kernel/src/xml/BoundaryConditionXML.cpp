//$Id: BoundaryConditionXML.cpp,v 1.5 2004/07/29 14:25:41 jbarbier Exp $
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
//$Log: BoundaryConditionXML.cpp,v $
//Revision 1.5  2004/07/29 14:25:41  jbarbier
//- $Log$ and $Id$ added
//
