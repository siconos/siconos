#include "RelayNSL.h"

#include "check.h"

RelayNSL::RelayNSL(): NonSmoothLaw()
{
  c = 0.0;
  d = 0.0;
  this->nsLawType = RELAYNSLAW;
}

RelayNSL::RelayNSL(NonSmoothLawXML* nslawxml): NonSmoothLaw(nslawxml)
{
  c = 0.0;
  d = 0.0;
  this->nsLawType = RELAYNSLAW;
}

RelayNSL::RelayNSL(double c, double d)
{
  this->c = c;
  this->d = d;
  this->nsLawType = RELAYNSLAW;
}

RelayNSL::~RelayNSL()
{}

bool RelayNSL::isVerified(void) const
{
  bool res = false;

  // to do

  return res;
}

void RelayNSL::fillNonSmoothLawWithNonSmoothLawXML()
{
  IN("RelayNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
  NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  if (this->nslawxml != NULL)
  {
    this->c = (static_cast<RelayNSLXML*>(this->nslawxml))->getC();
    this->d = (static_cast<RelayNSLXML*>(this->nslawxml))->getD();

    //    this->display();
  }
  else RuntimeException::selfThrow("RelayNSL::fillNonSmoothLawWithNonSmoothLawXML - object NonSmoothLawXML does not exist");
  OUT("RelayNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
}

void RelayNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the RelayNSL" << endl;
  cout << "| c : " << this->c << endl;
  cout << "| d : " << this->d << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void RelayNSL::saveNonSmoothLawToXML()
{
  IN("RelayNSL::saveNonSmoothLawToXML\n");
  static_cast<RelayNSLXML*>(this->nslawxml)->setC(this->c);
  static_cast<RelayNSLXML*>(this->nslawxml)->setD(this->d);
  //NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  OUT("RelayNSL::saveNonSmoothLawToXML\n");
}

void RelayNSL::createNonSmoothLaw(RelayNSLXML * nslawXML, double c, double d)//, Interaction * interaction)
{
  if (nslawXML != NULL)
  {
    this->nslawxml = nslawXML;
    this->nsLawType = RELAYNSLAW;
    this->fillNonSmoothLawWithNonSmoothLawXML();
  }
  else
  {
    this->c = c;
    this->d = d;
  }
}


RelayNSL* RelayNSL::convert(NonSmoothLaw* nsl)
{
  cout << "RelayNSL::convert (NonSmoothLaw* nsl)" << endl;
  RelayNSL* rnsl = dynamic_cast<RelayNSL*>(nsl);
  return rnsl;
}

