#include "Adams.h"
using namespace std;

// --- xml constructor ---
Adams::Adams(OneStepIntegratorXML* osixml): OneStepIntegrator(osixml), r(-1)
{
  integratorType = ADAMS_INTEGRATOR;
  if (osixml != NULL)
  {
    if ((static_cast<AdamsXML*>(integratorXml))->hasR() == true)
    {
      r = (static_cast<AdamsXML*>(integratorXml))->getR();
    }
  }
  else RuntimeException::selfThrow("Adams::Adams() - xml constructor - IntegratorXML object not exists");
}

// --- Minimum data constructor ---
Adams::Adams(TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(td, ds), r(-1)
{
  integratorType = ADAMS_INTEGRATOR;
}

// --- Destructor ---
Adams::~Adams()
{}

// --- Casting for Python ---
Adams* Adams::convert(OneStepIntegrator* osi)
{
  cout << "Adams::convert (OneStepIntegrator* osi)" << endl;
  Adams* adams = dynamic_cast<Adams*>(osi);
  return adams;
}

// --- Default constructor ---
Adams::Adams(): OneStepIntegrator(), r(-1)
{
  integratorType = ADAMS_INTEGRATOR;
}
