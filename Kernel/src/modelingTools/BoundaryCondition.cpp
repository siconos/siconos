
#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition()
{
  this->bcXML = NULL;
}

BoundaryCondition::BoundaryCondition(BoundaryConditionXML* bcxml)
{
  this->bcXML = bcxml;
}

BoundaryCondition::~BoundaryCondition()
{}


void BoundaryCondition::fillBCWithBCXML()
{
  OUT("BoundaryCondition::fillBCWithBCXML\n");
  if (this->bcXML != NULL)
  {}
  else RuntimeException::selfThrow("BoundaryCondition::fillBCWithBCXML - BoundaryConditionXML object not exists");
}

