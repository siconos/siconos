#include "OneStepIntegrator.h"
//#include "Strategy.h"

#include "check.h"

// --- Xml constructor ---
OneStepIntegrator::OneStepIntegrator(OneStepIntegratorXML* osixml): integratorType("none"), ds(0), r(1), timeDiscretisation(0), integratorxml(osixml)
{
  if (this->integratorxml != 0)
  {
    if (this->integratorxml->hasR()) this->r = this->integratorxml->getR();
  }
  else RuntimeException::selfThrow("OneStepIntegrator::fillIntegratorWithIntegratorXML - OneStepIntegratorXML object not exists");
}

// --- Constructor from a minimum set of data ---
OneStepIntegrator::OneStepIntegrator(TimeDiscretisation* td, DynamicalSystem* newDs): integratorType("none"), ds(newDs), r(1), timeDiscretisation(td), integratorxml(0)
{
  //to complete ...
}

/*OneStepIntegrator::OneStepIntegrator(OneStepIntegratorXML* osixml,TimeDiscretisation* td, DynamicalSystem* newDs): integratorType("none"), ds(newDs), r(1), timeDiscretisation(td), integratorxml(osixml)
{
  if(integratorxml != 0) if( this->integratorxml->hasR() ) this->r = this->integratorxml->getR();
  else RuntimeException::selfThrow("OneStepIntegrator::OneStepIntegrator() - xml constructor + data - OneStepIntegratorXML object not exists");
}
*/
// --- Destructor ---
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
  if (timeDiscretisation != 0) this->timeDiscretisation->display();
  cout << "-----------------------------------------------------" << endl << endl;
}

void OneStepIntegrator::saveIntegratorToXML()
{
  IN("OneStepIntegrator::saveIntegratorToXML\n");
  if (this->integratorxml != 0)
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

//-- Default constructor --
OneStepIntegrator::OneStepIntegrator(): integratorType("none"), ds(0), r(1), timeDiscretisation(0), integratorxml(0)
{}
