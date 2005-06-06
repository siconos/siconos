#include "LinearEC.h"
using namespace std;

LinearEC::LinearEC(): EqualityConstraint()
{
  this->type = LINEAREC;
}

LinearEC::LinearEC(EqualityConstraintXML *ecxml): EqualityConstraint(ecxml)
{
  this->type = LINEAREC;
}

LinearEC::~LinearEC()
{}

void LinearEC::createEqualityConstraint(EqualityConstraintXML *ecXML ,
                                        int number,  SiconosMatrix *G,
                                        vector<DSInputOutput*> *dsioVector)
{
  if (ecXML != NULL)
  {
    this->ecXML = ecXML;
    this->type = NLINEAREC;
    this->fillEqualityConstraintWithEqualityConstraintXML();
  }
  else
  {
    this->ecXML = NULL;
    this->type = NLINEAREC;
    this->number = number;
    this->G = *G;
    this->dsioVector = *dsioVector;
  }
}

