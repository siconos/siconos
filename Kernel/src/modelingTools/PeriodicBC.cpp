
#include "PeriodicBC.h"
#include "check.h"

PeriodicBC::PeriodicBC(): BoundaryCondition()
{
  this->boundaryType = PERIODICBC;
}

PeriodicBC::PeriodicBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = PERIODICBC;
}

PeriodicBC::~PeriodicBC()
{}

void PeriodicBC::fillBCWithBCXML()
{
  if (this->bcXML != NULL)
  {
    OUT("PeriodicBC::fillBCWithBCXML\n");
  }
  else RuntimeException::selfThrow("PeriodicBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void PeriodicBC::saveBCToXML()
{
  if (this->bcXML != NULL)
  {
    OUT("PeriodicBC::saveBCToXML\n");
  }
  else RuntimeException::selfThrow("PeriodicBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void PeriodicBC::createBoundaryCondition(BoundaryConditionXML * bcXML)//, DynamicalSystem* ds)
{
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = PERIODICBC;
    this->fillBCWithBCXML();
  }
  else
  {}
}


PeriodicBC* PeriodicBC::convert(BoundaryCondition* bc)
{
  cout << "PeriodicBC::convert (BoundaryCondition* bc)" << endl;
  PeriodicBC* pbc = dynamic_cast<PeriodicBC*>(bc);
  return pbc;
}

