#include "LinearTIEC.h"

LinearTIEC::LinearTIEC(): LinearEC()
{
  this->type = LINEARTIEC;
}

LinearTIEC::LinearTIEC(EqualityConstraintXML *ecxml): LinearEC(ecxml)
{
  this->type = LINEARTIEC;
}

LinearTIEC::~LinearTIEC()
{}

void LinearTIEC::createEqualityConstraint(EqualityConstraintXML *ecXML ,
    int number,  SiconosMatrix *G,
    vector<DSInputOutput*> *dsioVector)
{
  if (ecXML != NULL)
  {
    this->ecXML = ecXML;
    this->type = LINEARTIEC;
    this->fillEqualityConstraintWithEqualityConstraintXML();
  }
  else
  {
    this->ecXML = NULL;
    this->type = LINEARTIEC;
    this->number = number;
    this->G = *G;
    this->dsioVector = *dsioVector;
  }
}

