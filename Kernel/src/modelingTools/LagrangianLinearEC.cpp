
#include "LagrangianLinearEC.h"
#include "check.h"

LagrangianLinearEC::LagrangianLinearEC(): EqualityConstraint()
{
  this->type = LAGRANGIANEC;
}

LagrangianLinearEC::LagrangianLinearEC(EqualityConstraintXML *ecxml): EqualityConstraint(ecxml)
{
  this->type = LAGRANGIANEC;
}

LagrangianLinearEC::~LagrangianLinearEC()
{}

void LagrangianLinearEC::createEqualityConstraint(EqualityConstraintXML *ecXML,
    int number,  SiconosMatrix *G,
    vector<DSInputOutput*> *dsioVector)
{
  if (ecXML != NULL)
  {
    this->ecXML = ecXML;
    this->type = LAGRANGIANEC;
    this->fillEqualityConstraintWithEqualityConstraintXML();
  }
  else
  {
    this->ecXML = NULL;
    this->type = LAGRANGIANEC;
    this->number = number;
    this->G = *G;
    this->dsioVector = *dsioVector;
  }
}
