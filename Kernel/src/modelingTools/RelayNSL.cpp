#include "RelayNSL.h"
using namespace std;

RelayNSL::RelayNSL():
  NonSmoothLaw(), c(0.0), d(0.0)
{
  nsLawType = RELAYNSLAW;
}

RelayNSL::RelayNSL(NonSmoothLawXML* nslawxml):
  NonSmoothLaw(nslawxml), c(0.0), d(0.0)
{
  nsLawType = RELAYNSLAW;
  if (nslawxml != NULL)
  {
    c = (static_cast<RelayNSLXML*>(nslawxml))->getC();
    d = (static_cast<RelayNSLXML*>(nslawxml))->getD();
  }
  else RuntimeException::selfThrow("RelayNSL::xml constructor, xml file=NULL");
}

RelayNSL::RelayNSL(const double& newC, const double& newD):
  NonSmoothLaw(), c(newC), d(newD)
{
  nsLawType = RELAYNSLAW;
}

RelayNSL::~RelayNSL()
{}

bool RelayNSL::isVerified(void) const
{
  bool res = false;
  // to do
  return res;
}

void RelayNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the RelayNSL" << endl;
  cout << "| c : " << c << endl;
  cout << "| d : " << d << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void RelayNSL::saveNonSmoothLawToXML()
{
  IN("RelayNSL::saveNonSmoothLawToXML\n");
  static_cast<RelayNSLXML*>(nslawxml)->setC(c);
  static_cast<RelayNSLXML*>(nslawxml)->setD(d);
  OUT("RelayNSL::saveNonSmoothLawToXML\n");
}

RelayNSL* RelayNSL::convert(NonSmoothLaw* nsl)
{
  cout << "RelayNSL::convert (NonSmoothLaw* nsl)" << endl;
  RelayNSL* rnsl = dynamic_cast<RelayNSL*>(nsl);
  return rnsl;
}

