#include "ComplementarityConditionNSL.h"

#include "check.h"


ComplementarityConditionNSL::ComplementarityConditionNSL(): NonSmoothLaw()
{
  this->nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
}

ComplementarityConditionNSL::ComplementarityConditionNSL(NonSmoothLawXML* nslawxml): NonSmoothLaw(nslawxml)
{
  this->nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
}

ComplementarityConditionNSL::~ComplementarityConditionNSL()
{}

bool ComplementarityConditionNSL::isVerified(void) const
{
  bool res = false;

  // to do

  return res;
}

void ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML()
{
  IN("ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
  NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  if (this->nslawxml != NULL)
  {}
  else RuntimeException::selfThrow("ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML - ComplementarityConditionNSLXML object not exists");
  OUT("ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
}

void ComplementarityConditionNSL::saveNonSmoothLawToXML()
{
  IN("ComplementarityConditionNSL::saveNonSmoothLawToXML\n");
  NonSmoothLaw::saveNonSmoothLawToXML();
  if (this->nslawxml != NULL)
  {}
  else RuntimeException::selfThrow("ComplementarityConditionNSL::saveNonSmoothLawToXML - ComplementarityConditionNSLXML object not exists");
  OUT("ComplementarityConditionNSL::saveNonSmoothLawToXML\n");
}

void ComplementarityConditionNSL::createNonSmoothLaw(ComplementarityConditionNSLXML * nslawXML)//, Interaction * interaction)
{
  if (nslawXML != NULL)
  {
    this->nslawxml = nslawXML;
    this->nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
    this->fillNonSmoothLawWithNonSmoothLawXML();
  }
  else
  {}
}

ComplementarityConditionNSL* ComplementarityConditionNSL::convert(NonSmoothLaw* nsl)
{
  cout << "ComplementarityConditionNSL::convert (NonSmoothLaw* nsl)" << endl;
  ComplementarityConditionNSL* ccnsl = dynamic_cast<ComplementarityConditionNSL*>(nsl);
  return ccnsl;
}


