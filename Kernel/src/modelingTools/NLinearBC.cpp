
#include "NLinearBC.h"
using namespace std;

NLinearBC::NLinearBC(): BoundaryCondition()
{
  this->boundaryType = NLINEARBC;
}

NLinearBC::NLinearBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = NLINEARBC;
}

NLinearBC::~NLinearBC()
{}

void NLinearBC::fillBCWithBCXML()
{
  OUT("NLinearBC::fillBCWithBCXML\n");
  if (this->bcXML != NULL)
  {

  }
  else RuntimeException::selfThrow("NLinearBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void NLinearBC::saveBCToXML()
{
  OUT("NLinearBC::saveBCToXML\n");
  if (this->bcXML != NULL)
  {

  }
  else RuntimeException::selfThrow("NLinearBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void NLinearBC::createBoundaryCondition(BoundaryConditionXML * bcXML)
{
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = NLINEARBC;
    this->fillBCWithBCXML();
  }
  else
  {}
}

NLinearBC* NLinearBC::convert(BoundaryCondition* bc)
{
  cout << "NLinearBC::convert (BoundaryCondition* bc)" << endl;
  NLinearBC* nlbc = dynamic_cast<NLinearBC*>(bc);
  return nlbc;
}

