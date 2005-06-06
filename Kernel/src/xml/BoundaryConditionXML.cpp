#include "BoundaryConditionXML.h"
using namespace std;

BoundaryConditionXML::BoundaryConditionXML()
{}

BoundaryConditionXML::BoundaryConditionXML(xmlNode * rootBCNode)
{
  rootBCNode = rootBCNode;
}

BoundaryConditionXML::~BoundaryConditionXML()
{
}
