
#include "Adams.h"
#include "check.h"

Adams::Adams(): OneStepIntegrator()
{
  this->integratorType = ADAMS_INTEGRATOR;
  this->r = -1;
}

Adams::Adams(OneStepIntegratorXML* osixml, TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(osixml, td, ds)
{
  this->integratorType = ADAMS_INTEGRATOR;
}

Adams::~Adams()
{}

void Adams::saveIntegratorToXML()
{
  IN("Adams::saveIntegratorToXML\n");
  OneStepIntegrator::saveIntegratorToXML();
  if (this->integratorxml != NULL)
  {
    //(static_cast<AdamsXML*>(this->integratorxml))->setR( this->r );
  }
  else RuntimeException::selfThrow("Adams::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Adams::saveIntegratorToXML\n");
}

void Adams::createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation*td, DynamicalSystem* ds)//, Strategy * strategy)
{
  if (osiXML != NULL)
  {
    this->integratorxml = osiXML;
    this->timeDiscretisation = td;
    this->ds = ds;

    this->fillIntegratorWithIntegratorXML();
  }
  else
  {
    this->integratorxml = NULL;
    this->integratorType = ADAMS_INTEGRATOR;
    this->timeDiscretisation = td;
    this->ds = ds;
  }
}

void Adams::fillIntegratorWithIntegratorXML()
{
  IN("Adams::fillIntegratorWithIntegratorXML\n");
  OneStepIntegrator::fillIntegratorWithIntegratorXML();
  if (this->integratorxml != NULL)
  {
    if ((static_cast<AdamsXML*>(this->integratorxml))->hasR() == true)
    {
      this->r = (static_cast<AdamsXML*>(this->integratorxml))->getR();
    }
    else
    {
      cout << "Warning :  r is not defined in the XML file" << endl;
    }
  }
  else RuntimeException::selfThrow("Adams::fillIntegratorWithIntegratorXML - IntegratorXML object not exists");
  OUT("Adams::fillIntegratorWithIntegratorXML\n");
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

