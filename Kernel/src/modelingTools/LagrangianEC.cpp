
#include "LagrangianEC.h"
using namespace std;

LagrangianEC::LagrangianEC(): EqualityConstraint()
{
  this->type = LAGRANGIANEC;
}

LagrangianEC::LagrangianEC(EqualityConstraintXML *ecxml): EqualityConstraint(ecxml)
{
  this->type = LAGRANGIANEC;
}

LagrangianEC::~LagrangianEC()
{}

void LagrangianEC::createEqualityConstraint(EqualityConstraintXML *ecXML,
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
