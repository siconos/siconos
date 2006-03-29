/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "Strategy.h"
// includes to be deleted thanks to factories:
#include "Moreau.h"
#include "Lsodar.h"
#include "LCP.h"
#include "FrictionContact2D.h"
#include "FrictionContact3D.h"
#include "QP.h"
#include "Relay.h"



using namespace std;

// --- Default constructor ---
Strategy::Strategy(Model * mainModel):
  name("unnamed"), strategyType("undefined"), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(NULL), model(mainModel), isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{

  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model = NULL.");
  model->setStrategyPtr(this);
  //  isStrategyComplete();
}

// --- Default constructor ---
Strategy::Strategy(Model& newModel):
  name("unnamed"), strategyType("undefined"), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(NULL), model(NULL), isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  model = &newModel;
  model->setStrategyPtr(this);
  //  isStrategyComplete();
}

// --- Constructors from a given set of data ---

// Vector of OSI + OSNS + Model
Strategy::Strategy(vector<OneStepIntegrator*> newOsiVector, OneStepNSProblem* newNspb, Model* newModel):
  name("unnamed"), strategyType("undefined"), timeDiscretisation(NULL), integratorVector(newOsiVector), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model = NULL.");
  model->setStrategyPtr(this);
  //isStrategyComplete();
}

// Vector of OSI + Model
Strategy::Strategy(vector<OneStepIntegrator*> newOsiVector, Model* newModel):
  name("unnamed"), strategyType("undefined"), timeDiscretisation(NULL), integratorVector(newOsiVector), nsProblem(NULL), strategyxml(NULL), model(newModel),
  isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model = NULL.");
  model->setStrategyPtr(this);
  //  isStrategyComplete();
}

// OSNS + Model
Strategy::Strategy(OneStepNSProblem* newNspb, Model* newModel):
  name("unnamed"), strategyType("undefined"), timeDiscretisation(NULL), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model = NULL.");
  model->setStrategyPtr(this);
  //isStrategyComplete();
}

// --- xml constructor ---
Strategy::Strategy(StrategyXML* strxml, Model *newModel): strategyType("undefined"), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(strxml), model(newModel),
  isTimeDiscretisationAllocatedIn(true), isNSProblemAllocatedIn(false)
{
  IN("Strategy::xml constructor\n");
  if (strategyxml != 0)
  {
    // memory allocation/construction for time discretisation
    timeDiscretisation = new TimeDiscretisation(strategyxml->getTimeDiscretisationXMLPtr(), this);

    if (model == NULL)
      RuntimeException::selfThrow("Strategy constructor - model = NULL.");
    model->setStrategyPtr(this);
    int dsNb;
    DynamicalSystem *dsPtr;

    // --- OneStepIntegrators ---
    // Get the OSI vector from xml
    vector<OneStepIntegratorXML*> osiXMLVector = strategyxml->getOneStepIntegratorXML();
    unsigned int sizeOsi = osiXMLVector.size();
    integratorVector.reserve(sizeOsi);

    // For each OSI ...
    for (unsigned int i = 0; i < sizeOsi; i++)
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
      else RuntimeException::selfThrow("Strategy::xml constructor - wrong type of Integrator");
    }

    // OneStepNSProblem
    if (strategyxml->hasOneStepNSProblemXML())
    {
      // we get all the numbers of the Interactions to link
      //        vector<int> interactionNumbers = strategyxml->getOneStepNSProblemXMLPtr()->getInteractionConcerned();

      // OneStepNSProblem - LCP memory allocation/construction
      string type = strategyxml->getOneStepNSProblemXMLPtr()->getNSProblemType();
      if (type == LCP_TAG)
        nsProblem = new LCP(strategyxml->getOneStepNSProblemXMLPtr(), this);
      // OneStepNSProblem - FrictionContact2D
      else if (type == FrictionContact2D_TAG)
        nsProblem = new FrictionContact2D(strategyxml->getOneStepNSProblemXMLPtr(), this);
      // OneStepNSProblem - FrictionContact3D
      else if (type == FrictionContact3D_TAG)
        nsProblem = new FrictionContact3D(strategyxml->getOneStepNSProblemXMLPtr(), this);
      // Uncomment following lines when QP and Relay will be fully implemented
      //        // OneStepNSProblem - QP
      //        else if( strategyxml->getOneStepNSProblemXMLPtr()->getType() == QP_TAG )
      //    {
      //      // QP memory allocation/construction
      //      nsProblem = new QP(strategyxml->getOneStepNSProblemXMLPtr());
      //      for(unsigned int i = 0; i < interactionNumbers.size(); i++ )
      //        nsProblem->addInteraction( model->getNonSmoothDynamicalSystemPtr()->getInteractionPtrNumber(interactionNumbers[i]) );
      //      nsProblem->setStrategy(this);
      //    }
      //        // OneStepNSProblem - Relay
      //        else if( strategyxml->getOneStepNSProblemXMLPtr()->getType() == RELAY_TAG )
      //    {
      //      // relay memory allocation/construction
      //      nsProblem = new Relay(strategyxml->getOneStepNSProblemXMLPtr());
      //      for( unsigned int i = 0; i < interactionNumbers.size(); i++ )
      //        nsProblem->addInteraction( model->getNonSmoothDynamicalSystemPtr()->getInteractionPtrNumber(interactionNumbers[i]) );
      //      nsProblem->setStrategy(this);
      //    }
      else RuntimeException::selfThrow("Strategy::xml constructor - wrong type of NSProblem: inexistant or not yet implemented");
      isNSProblemAllocatedIn = true;
    }

    //isStrategyComplete();
  }
  else  RuntimeException::selfThrow("Strategy:: xml constructor - xml file = NULL");
  OUT("Strategy::xml constructor\n");
}

// --- Destructor ---
Strategy::~Strategy()
{
  if (isTimeDiscretisationAllocatedIn) delete timeDiscretisation;
  timeDiscretisation = NULL;
  if (isNSProblemAllocatedIn) delete nsProblem;
  nsProblem = NULL;
  if (integratorVector.size() > 0)
  {
    for (unsigned int i = 0; i < integratorVector.size(); i++)
    {
      if (isIntegratorVectorAllocatedIn[i]) delete integratorVector[i];
      integratorVector[i] = NULL;
    }
    integratorVector.clear();
  }
}

// Check whether strategy is complete or not
bool Strategy::isStrategyComplete() const
{
  bool isComplete = 1;
  //if (integratorVector.size()==1 && integratorVector[0] == NULL)
  //  cout << "Warning: strategy may be incomplete: no integrator" << endl;isComplete =0;
  //if (nsProblem == NULL) cout << "Warning: strategy may be incomplete: no NS problem" << endl;isComplete =0;
  return(isComplete);
}

// Getters/setters

void Strategy::setTimeDiscretisationPtr(TimeDiscretisation* td)
{
  // Warning: this function may be used carefully because of the links between Model and TimeDiscretisation
  // The strategy of the td input MUST be the current strategy.
  //
  if (isTimeDiscretisationAllocatedIn) delete timeDiscretisation;
  timeDiscretisation = td;
  isTimeDiscretisationAllocatedIn = false;
}

void Strategy::setOneStepNSProblemPtr(OneStepNSProblem* nspb)
{
  if (isNSProblemAllocatedIn) delete nsProblem;
  nsProblem = nspb;
  isNSProblemAllocatedIn = false;
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
    integratorVector[i]->computeFreeState();
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
  else RuntimeException::selfThrow("Strategy - computeOneStepNSProblem, OneStepNSProblem == NULL ");
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

void Strategy::newtonSolve(const double& criterion, const unsigned int& maxStep)
{

  bool isNewtonConverge = false;
  unsigned int nbNewtonStep = 0; // number of Newton iterations
  while ((!isNewtonConverge) && (nbNewtonStep <= maxStep))
  {
    nbNewtonStep++;
    computeFreeState();
    computeOneStepNSProblem();
    update();
    isNewtonConverge = newtonCheckConvergence(criterion);
  }
  if (!isNewtonConverge)
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
      else if (integratorVector[i]->getType() == LSODAR_INTEGRATOR)
        (static_cast<Lsodar*>(integratorVector[i]))->saveIntegratorToXML();
      else RuntimeException::selfThrow("Strategy::saveStrategyToXML - wrong type of OneStepIntegrator");
    }

    if (getStrategyXMLPtr()->hasOneStepNSProblemXML())
    {
      string NSPType = nsProblem->getType();
      if (NSPType == "LCP")
        (static_cast<LCP*>(nsProblem))->saveNSProblemToXML();
      else if (NSPType == "FrictionContact2D" || NSPType == "FrictionContact3D")
        (static_cast<FrictionContact*>(nsProblem))->saveNSProblemToXML();
      else if (NSPType == "QP")
        (static_cast<QP*>(nsProblem))->saveNSProblemToXML();
      else if (NSPType == "Relay")
        (static_cast<Relay*>(nsProblem))->saveNSProblemToXML();
      else
        RuntimeException::selfThrow("Strategy::saveStrategyToXML - wrong type of OneStepNSProblem");
    }
  }
  else RuntimeException::selfThrow("Strategy::saveStrategyToXML - StrategyXML = NULL");
  OUT("Strategy::saveStrategyToXML\n");
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
