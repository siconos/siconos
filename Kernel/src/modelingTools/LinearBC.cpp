
#include "LinearBC.h"
#include "check.h"

LinearBC::LinearBC(): BoundaryCondition()
{
  IN("LinearBC::LinearBC()\n");
  this->boundaryType = LINEARBC;
  this->omega = /*SiconosVector*/SimpleVector::/*SiconosVector*/SimpleVector();
  this->omega0 = SiconosMatrix::SiconosMatrix();
  this->omegaT = SiconosMatrix::SiconosMatrix();
  OUT("LinearBC::LinearBC()\n");
}

LinearBC::LinearBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = LINEARBC;
}

LinearBC::~LinearBC()
{}


void LinearBC::fillBCWithBCXML()
{
  OUT("LinearBC::fillBCWithBCXML\n");
  if (this->bcXML != NULL)
  {
    this->omega = static_cast<LinearBCXML*>(this->bcXML)->getOmega();
    this->omegaT = static_cast<LinearBCXML*>(this->bcXML)->getOmegaT();
    this->omega0 = static_cast<LinearBCXML*>(this->bcXML)->getOmega0();
  }
  else RuntimeException::selfThrow("LinearBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void LinearBC::saveBCToXML()
{
  OUT("LinearBC::saveBCToXML\n");
  if (this->bcXML != NULL)
  {
    static_cast<LinearBCXML*>(this->bcXML)->setOmega(&(this->omega));
    static_cast<LinearBCXML*>(this->bcXML)->setOmegaT(&(this->omegaT));
    static_cast<LinearBCXML*>(this->bcXML)->setOmega0(&(this->omega0));
  }
  else RuntimeException::selfThrow("LinearBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void LinearBC::createBoundaryCondition(BoundaryConditionXML * bcXML,
                                       SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
{
  IN("LinearBC::createBoundaryCondition\n");
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = LINEARBC;
    this->fillBCWithBCXML();
  }
  else if (omega != NULL && omega0 != NULL && omegaT != NULL)
  {
    this->omega = *omega;
    this->omega0 = *omega0;
    this->omegaT = *omegaT;
  }
  else RuntimeException::selfThrow("LinearBC::createBoundaryCondition - The omega, omega0 and/or omegaT matrices is/are missing");
  OUT("LinearBC::createBoundaryCondition\n");
}


LinearBC* LinearBC::convert(BoundaryCondition* bc)
{
  cout << "LinearBC::convert (BoundaryCondition* bc)" << endl;
  LinearBC* lbc = dynamic_cast<LinearBC*>(bc);
  return lbc;
}

