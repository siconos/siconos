//$Id: BoundaryCondition.cpp,v 1.5 2004/07/29 14:25:34 jbarbier Exp $

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

//$Log: BoundaryCondition.cpp,v $
//Revision 1.5  2004/07/29 14:25:34  jbarbier
//- $Log$ and $Id$ added
//
