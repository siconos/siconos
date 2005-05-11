#include "Adams.h"
#include "check.h"


Adams::Adams(OneStepIntegratorXML* osixml): OneStepIntegrator(osixml), r(-1)
{
  this->integratorType = ADAMS_INTEGRATOR;
  if (osixml != 0)
  {
    if ((static_cast<AdamsXML*>(this->integratorxml))->hasR() == true)
    {
      this->r = (static_cast<AdamsXML*>(this->integratorxml))->getR();
    }
  }
  else RuntimeException::selfThrow("Adams::Adams() - xml constructor - IntegratorXML object not exists");
}

Adams::Adams(TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(td, ds), r(-1)
{
  this->integratorType = ADAMS_INTEGRATOR;
}

Adams::~Adams()
{}

void Adams::saveIntegratorToXML()
{
  IN("Adams::saveIntegratorToXML\n");
  OneStepIntegrator::saveIntegratorToXML();
  if (this->integratorxml != 0)
  {
    //(static_cast<AdamsXML*>(this->integratorxml))->setR( this->r );
  }
  else RuntimeException::selfThrow("Adams::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Adams::saveIntegratorToXML\n");
}

void Adams::initialize()
{
  IN("Adams::initialize\n");
  OneStepIntegrator::initialize();
  OUT("Adams::initialize\n");
}


Adams* Adams::convert(OneStepIntegrator* osi)
{
  cout << "Adams::convert (OneStepIntegrator* osi)" << endl;
  Adams* adams = dynamic_cast<Adams*>(osi);
  return adams;
}

// Default constructor
Adams::Adams(): OneStepIntegrator(), r(-1)
{
  this->integratorType = ADAMS_INTEGRATOR;
}
