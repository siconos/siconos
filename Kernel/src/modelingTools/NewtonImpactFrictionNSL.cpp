#include "NewtonImpactFrictionNSL.h"
using namespace std;

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL():
  NonSmoothLaw(), en(0.0), et(0.0), mu(0.0)
{
  nsLawType = NEWTONIMPACTFRICTIONNSLAW;
}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(NonSmoothLawXML* nslawxml):
  NonSmoothLaw(nslawxml), en(0.0), et(0.0), mu(0.0)
{
  nsLawType = NEWTONIMPACTFRICTIONNSLAW;
  if (nslawxml != NULL)
  {
    en = (static_cast<NewtonImpactFrictionNSLXML*>(nslawxml))->getEn();
    et = (static_cast<NewtonImpactFrictionNSLXML*>(nslawxml))->getEt();
    mu = (static_cast<NewtonImpactFrictionNSLXML*>(nslawxml))->getMu();
  }
  else RuntimeException::selfThrow("NewtonImpactFrictionNSL:: xml constructor, xml file=NULL");
}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(const double& newEn, const double& newEt, const double& newMu):
  NonSmoothLaw(), en(newEn), et(newEt), mu(newMu)
{
  nsLawType = NEWTONIMPACTFRICTIONNSLAW;
}

NewtonImpactFrictionNSL::~NewtonImpactFrictionNSL()
{}

bool NewtonImpactFrictionNSL::isVerified(void) const
{
  bool res = false;
  // to do
  return res;
}

void NewtonImpactFrictionNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the NewtonImpactFrictionNSL" << endl;
  cout << "| The normal Newton coefficient of restitution en : " << en << endl;
  cout << "| The tangential Newton coefficient of restitution et : " << et << endl;
  cout << "| The friction coefficient mu : " << mu << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void NewtonImpactFrictionNSL::saveNonSmoothLawToXML()
{
  IN("NewtonImpactFrictionNSL::saveNonSmoothLawToXML\n");
  static_cast<NewtonImpactFrictionNSLXML*>(nslawxml)->setEn(en);
  static_cast<NewtonImpactFrictionNSLXML*>(nslawxml)->setEt(et);
  static_cast<NewtonImpactFrictionNSLXML*>(nslawxml)->setMu(mu);
  OUT("NewtonImpactFrictionNSL::saveNonSmoothLawToXML\n");
}

NewtonImpactFrictionNSL* NewtonImpactFrictionNSL::convert(NonSmoothLaw* nsl)
{
  cout << "NewtonImpactFrictionNSL::convert (NonSmoothLaw* nsl)" << endl;
  NewtonImpactFrictionNSL* nilnsl = dynamic_cast<NewtonImpactFrictionNSL*>(nsl);
  return nilnsl;
}


