#include "NewtonImpactFrictionNSL.h"

#include "check.h"

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(): NonSmoothLaw()
{
  this->en = 0.0;
  this->et = 0.0;
  this->mu = 0.0;
  this->nsLawType = NEWTONIMPACTFRICTIONNSLAW;
}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(NonSmoothLawXML* nslawxml): NonSmoothLaw(nslawxml)
{
  this->en = 0.0;
  this->et = 0.0;
  this->mu = 0.0;
  this->nsLawType = NEWTONIMPACTFRICTIONNSLAW;
}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(double en, double et, double mu)
{
  this->en = en;
  this->et = et;
  this->mu = mu;
  this->nsLawType = NEWTONIMPACTFRICTIONNSLAW;
}

NewtonImpactFrictionNSL::~NewtonImpactFrictionNSL()
{}

bool NewtonImpactFrictionNSL::isVerified(void) const
{
  bool res = false;

  // to do

  return res;
}


void NewtonImpactFrictionNSL::fillNonSmoothLawWithNonSmoothLawXML()
{
  IN("NewtonImpactFrictionNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
  NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  if (this->nslawxml != NULL)
  {
    this->en = (static_cast<NewtonImpactFrictionNSLXML*>(this->nslawxml))->getEn();
    this->et = (static_cast<NewtonImpactFrictionNSLXML*>(this->nslawxml))->getEt();
    this->mu = (static_cast<NewtonImpactFrictionNSLXML*>(this->nslawxml))->getMu();
    //    this->display();
  }
  else RuntimeException::selfThrow("NewtonImpactFrictionNSL::fillNonSmoothLawWithNonSmoothLawXML - object NonSmoothLawXML does not exist");
  OUT("NewtonImpactFrictionNSL::fillNonSmoothLawWithNonSmoothLawXML\n");

}

void NewtonImpactFrictionNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the NewtonImpactFrictionNSL" << endl;
  cout << "| The normal Newton coefficient of restitution en : " << this->en << endl;
  cout << "| The tangential Newton coefficient of restitution et : " << this->et << endl;
  cout << "| The friction coefficient mu : " << this->mu << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void NewtonImpactFrictionNSL::saveNonSmoothLawToXML()
{
  IN("NewtonImpactFrictionNSL::saveNonSmoothLawToXML\n");
  static_cast<NewtonImpactFrictionNSLXML*>(this->nslawxml)->setEn(this->en);
  static_cast<NewtonImpactFrictionNSLXML*>(this->nslawxml)->setEt(this->et);
  static_cast<NewtonImpactFrictionNSLXML*>(this->nslawxml)->setMu(this->mu);
  OUT("NewtonImpactFrictionNSL::saveNonSmoothLawToXML\n");
}

void NewtonImpactFrictionNSL::createNonSmoothLaw(NewtonImpactFrictionNSLXML * nslawXML, double en, double et, double mu)
{
  if (nslawXML != NULL)
  {
    this->nslawxml = nslawXML;
    this->nsLawType = NEWTONIMPACTFRICTIONNSLAW;
    this->en = 0.0;
    this->et = 0.0;
    this->mu = 0.0;
    this->fillNonSmoothLawWithNonSmoothLawXML();
  }
  else
  {
    this->en = en;
    this->et = et;
    this->mu = mu;
  }
}


NewtonImpactFrictionNSL* NewtonImpactFrictionNSL::convert(NonSmoothLaw* nsl)
{
  cout << "NewtonImpactFrictionNSL::convert (NonSmoothLaw* nsl)" << endl;
  NewtonImpactFrictionNSL* nilnsl = dynamic_cast<NewtonImpactFrictionNSL*>(nsl);
  return nilnsl;
}


