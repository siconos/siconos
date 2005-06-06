#include "NewtonImpactLawNSL.h"
using namespace std;

NewtonImpactLawNSL::NewtonImpactLawNSL(): NonSmoothLaw(), e(0.0)
{
  nsLawType = NEWTONIMPACTLAWNSLAW;
}

NewtonImpactLawNSL::NewtonImpactLawNSL(NonSmoothLawXML* nslawxml):
  NonSmoothLaw(nslawxml), e(0.0)
{
  nsLawType = NEWTONIMPACTLAWNSLAW;
  if (nslawxml != NULL)
    e = (static_cast<NewtonImpactLawNSLXML*>(nslawxml))->getE();
  else RuntimeException::selfThrow("NewtonImpactLawNSL:: xml constructor, xml file=NULL");
}

NewtonImpactLawNSL::NewtonImpactLawNSL(const double& newE):
  NonSmoothLaw(), e(newE)
{
  nsLawType = NEWTONIMPACTLAWNSLAW;
}

NewtonImpactLawNSL::~NewtonImpactLawNSL()
{}

bool NewtonImpactLawNSL::isVerified() const
{
  bool res = false;
  // to do
  return res;
}

void NewtonImpactLawNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the NewtonImpactLawNSL" << endl;
  cout << "| The Newton coefficient of restitution e : " << e << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void NewtonImpactLawNSL::saveNonSmoothLawToXML()
{
  IN("NewtonImpactLawNSL::saveNonSmoothLawToXML\n");
  static_cast<NewtonImpactLawNSLXML*>(this->nslawxml)->setE(e);
  OUT("NewtonImpactLawNSL::saveNonSmoothLawToXML\n");
}

NewtonImpactLawNSL* NewtonImpactLawNSL::convert(NonSmoothLaw* nsl)
{
  cout << "NewtonImpactLawNSL::convert (NonSmoothLaw* nsl)" << endl;
  NewtonImpactLawNSL* nilnsl = dynamic_cast<NewtonImpactLawNSL*>(nsl);
  return nilnsl;
}


