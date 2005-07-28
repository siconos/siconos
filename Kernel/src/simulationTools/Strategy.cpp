#include "Strategy.h"
// includes to be deleted thanks to factories:
#include "Moreau.h"
#include "Lsodar.h"
#include "Adams.h"
#include "LCP.h"
#include "CFD.h"
#include "QP.h"
#include "Relay.h"



using namespace std;

// --- Default constructor ---
Strategy::Strategy(): strategyType("none"), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(NULL), model(NULL),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

// --- Constructors from a given set of data ---

Strategy::Strategy(TimeDiscretisation* newTd, vector<OneStepIntegrator*> newOsiVector, OneStepNSProblem* newNspb, Model* newModel):
  strategyType("none"), timeDiscretisation(newTd), integratorVector(newOsiVector), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

Strategy::Strategy(TimeDiscretisation* newTd, vector<OneStepIntegrator*> newOsiVector, Model* newModel):
  strategyType("none"), timeDiscretisation(newTd), integratorVector(newOsiVector), nsProblem(NULL), strategyxml(NULL), model(newModel),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

Strategy::Strategy(TimeDiscretisation* newTd, OneStepNSProblem* newNspb, Model* newModel):
  strategyType("none"), timeDiscretisation(newTd), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

Strategy::Strategy(TimeDiscretisation* newTd, Model* newModel):
  strategyType("none"), timeDiscretisation(newTd), nsProblem(NULL), strategyxml(NULL), model(newModel),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

Strategy::Strategy(vector<OneStepIntegrator*> newOsiVector, OneStepNSProblem* newNspb, Model* newModel):
  strategyType("none"), timeDiscretisation(NULL), integratorVector(newOsiVector), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

Strategy::Strategy(vector<OneStepIntegrator*> newOsiVector, Model* newModel):
  strategyType("none"), timeDiscretisation(NULL), integratorVector(newOsiVector),
  nsProblem(NULL), strategyxml(NULL), model(newModel),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

Strategy::Strategy(OneStepNSProblem* newNspb, Model* newModel):
  strategyType("none"), timeDiscretisation(NULL), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscrAllocatedIn(false), isNsPbAllocatedIn(false)
{
  isStrategyComplete();
}

// --- xml constructor ---
Strategy::Strategy(StrategyXML* strxml, Model *newModel): strategyType("none"), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(strxml), model(newModel),
  isTimeDiscrAllocatedIn(true), isNsPbAllocatedIn(true)
{
  IN("Strategy::xml constructor\n");
  if (strategyxml != 0)
  {
    // memory allocation/construction for time discretisation
    timeDiscretisation = new TimeDiscretisation(strategyxml->getTimeDiscretisationXMLPtr(), this);

    if (model != NULL)
    {
      int dsNb;
      DynamicalSystem *dsPtr;

      // --- OneStepIntegrators ---
      // Get the OSI vector from xml
      vector<OneStepIntegratorXML*> osiXMLVector = strategyxml->getOneStepIntegratorXML();

      // For each OSI ...
      for (unsigned int i = 0; i < osiXMLVector.size(); i++)
      {
        // Get the number of the DynamicalSystem concerned
        dsNb = (osiXMLVector[i]->getDSConcerned())[0];
        // Get this DS
        dsPtr = model->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtrNumber(dsNb);

        if (dsPtr == NULL) RuntimeException::selfThrow("Strategy::xml constructor - DS = NULL");
        // memory allocation/construction of the osi
        // Moreau
        if (osiXMLVector[i]->getType() == MOREAU_TAG)
        {
          integratorVector.push_back(new Moreau(osiXMLVector[i], timeDiscretisation, dsPtr));
          isIntegratorVectorAllocatedIn.push_back(true);
        }
        // Lsodar
        else if (osiXMLVector[i]->getType() == LSODAR_TAG)
        {
          integratorVector.push_back(new Lsodar(osiXMLVector[i]));
          isIntegratorVectorAllocatedIn.push_back(true);
        }
        // Adams
        else if (osiXMLVector[i]->getType() == ADAMS_TAG)
        {
          integratorVector.push_back(new Adams(osiXMLVector[i]));
          isIntegratorVectorAllocatedIn.push_back(true);
        }
        else RuntimeException::selfThrow("Strategy::xml constructor - wrong type of Integrator");
      }

      // OneStepNSProblem
      if (strategyxml->hasOneStepNSProblemXML())
      {
        // we get all the numbers of the Interactions to link
        vector<int> interactionNumbers = strategyxml->getOneStepNSProblemXMLPtr()->getInteractionConcerned();

        // OneStepNSProblem - LCP memory allocation/construction
        if (strategyxml->getOneStepNSProblemXMLPtr()->getType() == LCP_TAG)
          nsProblem = new LCP(strategyxml->getOneStepNSProblemXMLPtr(), this);
        // OneStepNSProblem - CFD
        else if (strategyxml->getOneStepNSProblemXMLPtr()->getType() == CFD_TAG)
        {
          // CFD memory allocation/construction
          nsProblem = new CFD(strategyxml->getOneStepNSProblemXMLPtr());
          for (unsigned int i = 0; i < interactionNumbers.size(); i++)
            nsProblem->addInteraction(model->getNonSmoothDynamicalSystemPtr()->getInteractionPtrNumber(interactionNumbers[i]));
          nsProblem->setStrategy(this);
        }
        // OneStepNSProblem - QP
        else if (strategyxml->getOneStepNSProblemXMLPtr()->getType() == QP_TAG)
        {
          // QP memory allocation/construction
          nsProblem = new QP(strategyxml->getOneStepNSProblemXMLPtr());
          for (unsigned int i = 0; i < interactionNumbers.size(); i++)
            nsProblem->addInteraction(model->getNonSmoothDynamicalSystemPtr()->getInteractionPtrNumber(interactionNumbers[i]));
          nsProblem->setStrategy(this);
        }
        // OneStepNSProblem - Relay
        else if (strategyxml->getOneStepNSProblemXMLPtr()->getType() == RELAY_TAG)
        {
          // relay memory allocation/construction
          nsProblem = new Relay(strategyxml->getOneStepNSProblemXMLPtr());
          for (unsigned int i = 0; i < interactionNumbers.size(); i++)
            nsProblem->addInteraction(model->getNonSmoothDynamicalSystemPtr()->getInteractionPtrNumber(interactionNumbers[i]));
          nsProblem->setStrategy(this);
        }
        else RuntimeException::selfThrow("Strategy::xml constructor - wrong type of NSProblem");
      }
    }
    isStrategyComplete();
  }
  else  RuntimeException::selfThrow("Strategy:: xml constructor - xml file = NULL");
  OUT("Strategy::xml constructor\n");
}

// --- Destructor ---
Strategy::~Strategy()
{
  if (isTimeDiscrAllocatedIn)
  {
    delete timeDiscretisation;
    timeDiscretisation = NULL;
  }
  if (isNsPbAllocatedIn)
  {
    delete nsProblem;
    nsProblem = NULL;
  }
  if (integratorVector.size() > 0)
  {
    for (unsigned int i = 0; i < integratorVector.size(); i++)
    {
      if (isIntegratorVectorAllocatedIn[i])
      {
        delete integratorVector[i];
        integratorVector[i] = NULL;
      }
    }
    integratorVector.clear();
  }
}

// Check whether strategy is complete or not

bool Strategy::isStrategyComplete() const
{
  bool isComplete = 1;
  if (timeDiscretisation == NULL) cout << "Warning: strategy may be incomplete: no time discretisation" << endl;
  isComplete = 0;
  if (integratorVector.size() == 1 && integratorVector[0] == NULL)
    cout << "Warning: strategy may be incomplete: no integrator" << endl;
  isComplete = 0;
  if (nsProblem == NULL) cout << "Warning: strategy may be incomplete: no NS problem" << endl;
  isComplete = 0;
  if (model == NULL) cout << "Warning: strategy may be incomplete: not linked with any model" << endl;
  isComplete = 0;
  return(isComplete);
}

// Getters/setters

void Strategy::setTimeDiscretisationPtr(TimeDiscretisation* td)
{
  if (isTimeDiscrAllocatedIn) delete timeDiscretisation;
  timeDiscretisation = td;
  isTimeDiscrAllocatedIn = false;
}

void Strategy::setOneStepNSProblemPtr(OneStepNSProblem* nspb)
{
  if (isNsPbAllocatedIn) delete nsProblem;
  nsProblem = nspb;
  isNsPbAllocatedIn = false;
}

void Strategy::setOneStepIntegrators(const vector<OneStepIntegrator*> vOSI)
{
  for (unsigned int i = 0; i < integratorVector.size(); i++)
  {
    if (isIntegratorVectorAllocatedIn[i])
    {
      delete integratorVector[i];
      integratorVector[i] = NULL;
    }
  }
  integratorVector.clear();
  integratorVector = vOSI;
  isIntegratorVectorAllocatedIn.clear();
  isIntegratorVectorAllocatedIn.resize(vOSI.size(), false);
};

OneStepIntegrator* Strategy::getOneStepIntegrator(const int& nb) const
{
  if ((unsigned int)nb >= integratorVector.size())
    RuntimeException::selfThrow("Strategy - getIntegrator : \'nb\' is out of range");
  return integratorVector[nb];
}

void Strategy::computeFreeState()
{
  IN("Strategy::computeFreeState\n");
  for (unsigned int i = 0; i < integratorVector.size(); i++)
  {
    integratorVector[i]->computeFreeState();
  }
  OUT("Strategy::computeFreeState\n");
}

void Strategy::nextStep()
{
  // increment time step
  timeDiscretisation->increment();
  model->setCurrentT(model->getCurrentT() + timeDiscretisation->getH());

  for (unsigned int i = 0; i < integratorVector.size(); i++)
    integratorVector[i]->nextStep();
  if (nsProblem != NULL)
    nsProblem->nextStep();

}

void Strategy::computeOneStepNSProblem()
{
  if (nsProblem != NULL) nsProblem->compute(model->getCurrentT());
}

void Strategy::update()
{
  // compute input (lambda -> r)
  if (nsProblem != NULL) nsProblem->updateInput();

  // compute state for each dynamical system
  for (unsigned int i = 0; i < integratorVector.size(); i++)
  {
    integratorVector[i]->updateState();
  }

  // compute output (y, ydot)
  if (nsProblem != NULL) nsProblem->updateOutput();

}

void Strategy::initialize()
{
  // initialization of the OneStepIntegrators
  for (unsigned int i = 0; i < integratorVector.size(); i++)
    integratorVector[i]->initialize();

  // initialization of  OneStepNonSmoothProblem
  if (nsProblem != NULL)
    nsProblem->initialize();
}

OneStepIntegrator* Strategy::getIntegratorOfDSPtr(const int& numberDS) const
{
  vector<OneStepIntegrator*>::const_iterator it  = integratorVector.begin();
  while ((*it)->getDynamicalSystemPtr()->getNumber() != numberDS && it != integratorVector.end())
    it++;

  if (it == integratorVector.end())
    RuntimeException::selfThrow("Strategy::getIntegratorOfDSPtr(numberDS), no integrator corresponds to this dynamical sytem");

  return (*it);
}

OneStepIntegrator* Strategy::getIntegratorOfDSPtr(DynamicalSystem * ds) const
{
  vector<OneStepIntegrator*>::const_iterator it = integratorVector.begin();
  while ((*it)->getDynamicalSystemPtr() != ds && it != integratorVector.end())
    it++;

  if (it == integratorVector.end())
    RuntimeException::selfThrow("Strategy::getIntegratorOfDSPtr(ds), no integrator corresponds to this dynamical sytem");

  return (*it);
}

void Strategy::newtonSolve(const double& criterion, const int& maxStep)
{

  bool isNewtonConverge = false;
  long nbNewtonStep = 0; // number of Newton iterations
  while ((!isNewtonConverge) && (nbNewtonStep <= maxStep))
  {
    nbNewtonStep++;
    computeFreeState();
    computeOneStepNSProblem();
    update();
    isNewtonConverge = newtonCheckConvergence(criterion);
  }
  if (isNewtonConverge)
    cout << "Newton process: convergence reached after " << nbNewtonStep << " Newtons steps" << endl;
  else
    cout << "Newton process stopped: reach max step number" << endl ;

  // time step increment
  model->setCurrentT(model->getCurrentT() + timeDiscretisation->getH());
}

bool Strategy::newtonCheckConvergence(const double& criterion)
{
  bool checkConvergence = false;
  // get the non smooth dynamical system
  NonSmoothDynamicalSystem* nsds = model-> getNonSmoothDynamicalSystemPtr();
  // get the nsds indicator of convergence
  double nsdsConverge = nsds -> nsdsConvergenceIndicator();
  if (nsdsConverge < criterion) checkConvergence = true ;

  return(checkConvergence);
}

void Strategy::saveStrategyToXML()
{
  IN("Strategy::saveStrategyToXML\n");
  if (strategyxml != NULL)
  {
    int size, i;
    size = integratorVector.size();
    for (i = 0; i < size; i++)
    {
      if (integratorVector[i]->getType() == MOREAU_INTEGRATOR)
        (static_cast<Moreau*>(integratorVector[i]))->saveIntegratorToXML();
      else if (integratorVector[i]->getType() == ADAMS_INTEGRATOR)
        (static_cast<Adams*>(integratorVector[i]))->saveIntegratorToXML();
      else if (integratorVector[i]->getType() == LSODAR_INTEGRATOR)
        (static_cast<Lsodar*>(integratorVector[i]))->saveIntegratorToXML();
      else RuntimeException::selfThrow("Strategy::saveStrategyToXML - wrong type of OneStepIntegrator");
    }

    if (getStrategyXMLPtr()->hasOneStepNSProblemXML())
    {
      if (nsProblem->getType() == LCP_OSNSP)
        (static_cast<LCP*>(nsProblem))->saveNSProblemToXML();
      else if (nsProblem->getType() == CFD_OSNSP)
        (static_cast<CFD*>(nsProblem))->saveNSProblemToXML();
      else if (nsProblem->getType() == QP_OSNSP)
        (static_cast<QP*>(nsProblem))->saveNSProblemToXML();
      else if (nsProblem->getType() == RELAY_OSNSP)
        (static_cast<Relay*>(nsProblem))->saveNSProblemToXML();
      else
        RuntimeException::selfThrow("Strategy::saveStrategyToXML - wrong type of OneStepNSProblem");
    }
  }
  else RuntimeException::selfThrow("Strategy::saveStrategyToXML - StrategyXML = NULL");
  OUT("Strategy::saveStrategyToXML\n");
}

OneStepIntegrator* Strategy::addAdams(TimeDiscretisation* td, DynamicalSystem* ds)
{
  if (hasDynamicalSystemIntegrator(ds))
    RuntimeException::selfThrow("Strategy::addAdams : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
  OneStepIntegrator* osi;
  osi = new Adams(td, ds);
  integratorVector.push_back(osi);
  isIntegratorVectorAllocatedIn.push_back(true);
  return osi;
}

OneStepIntegrator* Strategy::addMoreau(TimeDiscretisation* td, DynamicalSystem* ds, const double& theta)
{
  if (hasDynamicalSystemIntegrator(ds))
    RuntimeException::selfThrow("Strategy::addMoreau : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
  OneStepIntegrator* osi;
  osi = new Moreau(td, ds, theta);
  integratorVector.push_back(osi);
  isIntegratorVectorAllocatedIn.push_back(true);
  return osi;
}

OneStepIntegrator* Strategy::addLsodar(TimeDiscretisation* td, DynamicalSystem* ds)
{
  if (hasDynamicalSystemIntegrator(ds))
    RuntimeException::selfThrow("Strategy::addLsodar : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
  OneStepIntegrator* osi;
  osi = new Lsodar(td, ds);
  integratorVector.push_back(osi);
  isIntegratorVectorAllocatedIn.push_back(true);
  return osi;
}

bool Strategy::hasDynamicalSystemIntegrator(DynamicalSystem* ds) const
{
  for (unsigned int i = 0; i < integratorVector.size(); i++)
  {
    if (ds == integratorVector[i]->getDynamicalSystemPtr()) return true;
  }
  return false;
}

