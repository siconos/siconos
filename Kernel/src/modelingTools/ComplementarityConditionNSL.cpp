#include "ComplementarityConditionNSL.h"
using namespace std;

ComplementarityConditionNSL::ComplementarityConditionNSL(): NonSmoothLaw()
{
  nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
}

ComplementarityConditionNSL::ComplementarityConditionNSL(NonSmoothLawXML* nslawxml):
  NonSmoothLaw(nslawxml)
{
  nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
}

ComplementarityConditionNSL::~ComplementarityConditionNSL()
{}

bool ComplementarityConditionNSL::isVerified() const
{
  bool res = false;
  // to do
  return res;
}

ComplementarityConditionNSL* ComplementarityConditionNSL::convert(NonSmoothLaw* nsl)
{
  cout << "ComplementarityConditionNSL::convert (NonSmoothLaw* nsl)" << endl;
  ComplementarityConditionNSL* ccnsl = dynamic_cast<ComplementarityConditionNSL*>(nsl);
  return ccnsl;
}


