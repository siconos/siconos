#include "NewtonImpactLawNSL.h"

#include "check.h"

NewtonImpactLawNSL::NewtonImpactLawNSL(): NonSmoothLaw()
{
  this->e = 0.0;
  this->nsLawType = NEWTONIMPACTLAWNSLAW;
}

NewtonImpactLawNSL::NewtonImpactLawNSL(NonSmoothLawXML* nslawxml): NonSmoothLaw(nslawxml)
{
  this->e = 0.0;
  this->nsLawType = NEWTONIMPACTLAWNSLAW;
}

NewtonImpactLawNSL::NewtonImpactLawNSL(double e)
{
  this->e = e;
  this->nsLawType = NEWTONIMPACTLAWNSLAW;
}

NewtonImpactLawNSL::~NewtonImpactLawNSL()
{}

bool NewtonImpactLawNSL::isVerified(void) const
{
  bool res = false;

  // to do

  return res;
}


void NewtonImpactLawNSL::fillNonSmoothLawWithNonSmoothLawXML()
{
  IN("NewtonImpactLawNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
  NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  if (this->nslawxml != NULL)
  {
    this->e = (static_cast<NewtonImpactLawNSLXML*>(this->nslawxml))->getE();
    //    this->display();
  }
  else RuntimeException::selfThrow("NewtonImpactLawNSL::fillNonSmoothLawWithNonSmoothLawXML - object NonSmoothLawXML does not exist");
  OUT("NewtonImpactLawNSL::fillNonSmoothLawWithNonSmoothLawXML\n");

}

void NewtonImpactLawNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the NewtonImpactLawNSL" << endl;
  cout << "| The Newton coefficient of restitution e : " << this->e << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void NewtonImpactLawNSL::saveNonSmoothLawToXML()
{
  IN("NewtonImpactLawNSL::saveNonSmoothLawToXML\n");
  static_cast<NewtonImpactLawNSLXML*>(this->nslawxml)->setE(this->e);
  OUT("NewtonImpactLawNSL::saveNonSmoothLawToXML\n");
}

void NewtonImpactLawNSL::createNonSmoothLaw(NewtonImpactLawNSLXML * nslawXML, double e)//, Interaction * interaction)
{
  if (nslawXML != NULL)
  {
    this->nslawxml = nslawXML;
    this->nsLawType = NEWTONIMPACTLAWNSLAW;
    this->e = 0.0;
    this->fillNonSmoothLawWithNonSmoothLawXML();
  }
  else
  {
    this->e = e;
  }
}


NewtonImpactLawNSL* NewtonImpactLawNSL::convert(NonSmoothLaw* nsl)
{
  cout << "NewtonImpactLawNSL::convert (NonSmoothLaw* nsl)" << endl;
  NewtonImpactLawNSL* nilnsl = dynamic_cast<NewtonImpactLawNSL*>(nsl);
  return nilnsl;
}


