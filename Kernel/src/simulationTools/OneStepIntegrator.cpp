#include "OneStepIntegrator.h"
//#include "Strategy.h"

#include "check.h"

OneStepIntegrator::OneStepIntegrator()
{
  this->integratorxml = NULL;
  this->timeDiscretisation = NULL;
  this->ds = NULL;
  this->r = 1;
}

OneStepIntegrator::OneStepIntegrator(OneStepIntegratorXML* osixml, TimeDiscretisation* td, DynamicalSystem* ds)
{
  this->integratorxml = osixml;
  this->timeDiscretisation = td;
  this->ds = ds;
  cout << "  - the DynamicalSystem Linked to this OneStepIntegrator has the number " << this->ds->getNumber() << " and his id is " << this->ds->getId() << endl;
}

OneStepIntegrator::~OneStepIntegrator()
{}


void OneStepIntegrator::initialize()
{
  IN("OneStepIntegrator::initialize \n");
  this->ds->initMemory(this->r);
  OUT("OneStepIntegrator::initialize\n");
}



void OneStepIntegrator::integrate()
{
  IN("OneStepIntegrator::integrate\n");
  OUT("OneStepIntegrator::integrate\n");
}

void OneStepIntegrator::nextStep(void)
{
  IN("OneStepIntegrator::nextStep\n");
  this->ds->swapInMemory();
  OUT("OneStepIntegrator::nextStep\n");

}


void OneStepIntegrator::computeFreeState()
{
  // to do
}


void OneStepIntegrator::updateState()
{
  // to do
}


void OneStepIntegrator::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the OneStepIntegrator " << endl;
  cout << "| integratorType : " << this->integratorType << endl;
  cout << "| ds : " << this->ds->getId() << endl;
  if (this->integratorType != MOREAU_INTEGRATOR)
    cout << "| r : " << this->r << endl;
  this->timeDiscretisation->display();
  cout << "-----------------------------------------------------" << endl << endl;
}



void OneStepIntegrator::fillIntegratorWithIntegratorXML()
{
  IN("OneStepIntegrator::fillIntegratorWithIntegratorXML\n");
  if (this->integratorxml != NULL)
  {
    if (this->integratorxml->hasR()) this->r = this->integratorxml->getR();
    else cout << "Warning : the r value of the OneStepIntegrator is not defined, optional attribute." << endl;
  }
  else RuntimeException::selfThrow("OneStepIntegrator::fillIntegratorWithIntegratorXML - OneStepIntegratorXML object not exists");
  OUT("OneStepIntegrator::fillIntegratorWithIntegratorXML\n");

}

void OneStepIntegrator::saveIntegratorToXML()
{
  IN("OneStepIntegrator::saveIntegratorToXML\n");
  if (this->integratorxml != NULL)
  {
    vector<int> dsConcerned;
    dsConcerned.push_back(this->ds->getNumber());
    this->integratorxml->setDSConcerned(&dsConcerned);

    // r is saved only if the integrator is not a Moreau integrator !
    if (this->integratorType != MOREAU_INTEGRATOR) this->integratorxml->setR(this->r);
  }
  else RuntimeException::selfThrow("OneStepIntegrator::saveIntegratorToXML - OneStepIntegratorXML object not exists");
  OUT("OneStepIntegrator::saveIntegratorToXML\n");

}
