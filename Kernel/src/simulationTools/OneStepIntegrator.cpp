#include "OneStepIntegrator.h"
using namespace std;

// --- Xml constructor ---
OneStepIntegrator::OneStepIntegrator(OneStepIntegratorXML* osixml): integratorType("none"), ds(NULL), sizeMem(1), timeDiscretisation(NULL), integratorXml(osixml)
{
  if (integratorXml != NULL)
  {
    if (integratorXml->hasR()) sizeMem = integratorXml->getR();
  }
  else RuntimeException::selfThrow("OneStepIntegrator: xml constructor, xml object = NULL");
}

// --- Constructor from a minimum set of data ---
OneStepIntegrator::OneStepIntegrator(TimeDiscretisation* td, DynamicalSystem* newDs): integratorType("none"), ds(newDs), sizeMem(1), timeDiscretisation(td), integratorXml(NULL)
{}

// --- Destructor ---
OneStepIntegrator::~OneStepIntegrator()
{}


void OneStepIntegrator::initialize()
{
  ds->initMemory(sizeMem);
}

void OneStepIntegrator::nextStep()
{
  ds->swapInMemory();
}

void OneStepIntegrator::computeFreeState()
{
  RuntimeException::selfThrow("OneStepIntegrator:computeFreeState, not yet implemented for this type of integrator" + getType());
}

void OneStepIntegrator::integrate()
{
  RuntimeException::selfThrow("OneStepIntegrator:integrate, not yet implemented for this type of integrator" + getType());
}

void OneStepIntegrator::updateState()
{
  RuntimeException::selfThrow("OneStepIntegrator:updateState, not yet implemented for this type of integrator" + getType());
}

void OneStepIntegrator::display() const
{
  cout << "==== OneStepIntegrator display =====" << endl;
  cout << "| integratorType : " << integratorType << endl;
  cout << "| DS id is: " << endl;
  if (ds != NULL) cout << ds->getId() << endl;
  else cout << "-> NULL" << endl;
  cout << "| sizeMem: " << sizeMem << endl;
  if (timeDiscretisation != NULL) timeDiscretisation->display();
  else cout << "-> NULL" << endl;
  cout << "====================================" << endl;
}

void OneStepIntegrator::saveIntegratorToXML()
{
  IN("OneStepIntegrator::saveIntegratorToXML\n");
  if (integratorXml != 0)
  {
    vector<int> dsConcerned;
    dsConcerned.push_back(ds->getNumber());
    integratorXml->setDSConcerned(&dsConcerned);

    // r is saved only if the integrator is not a Moreau integrator !
    if (integratorType != MOREAU_INTEGRATOR) integratorXml->setR(sizeMem);
  }
  else RuntimeException::selfThrow("OneStepIntegrator::saveIntegratorToXML - OneStepIntegratorXML object = NULL");
  OUT("OneStepIntegrator::saveIntegratorToXML\n");

}

//-- Default constructor --
OneStepIntegrator::OneStepIntegrator(): integratorType("none"), ds(NULL), sizeMem(1), timeDiscretisation(NULL), integratorXml(NULL)
{}
