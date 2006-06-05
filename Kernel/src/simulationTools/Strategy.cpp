/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
Strategy::Strategy(Model * mainModel, const string& id):
  name("unnamed"), strategyType(id), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(NULL), model(mainModel), isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{

  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model == NULL.");
  model->setStrategyPtr(this);
  //  isStrategyComplete();
}

// --- Default constructor ---
Strategy::Strategy(Model& newModel, const string& id):
  name("unnamed"), strategyType(id), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(NULL), model(NULL), isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  model = &newModel;
  model->setStrategyPtr(this);
  //  isStrategyComplete();
}

// --- Constructors from a given set of data ---

// set of OSI + OSNS + Model
Strategy::Strategy(const OSISet& newOSIs, OneStepNSProblem* newNspb, Model* newModel, const string& id):
  name("unnamed"), strategyType(id), timeDiscretisation(NULL), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model == NULL.");
  model->setStrategyPtr(this);
  allOSI = newOSIs;
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
    isOSIAllocatedIn[*it] = false;
  //isStrategyComplete();
}

// Vector of OSI + Model
// Warning: copy of newOsiVector into integratorVector. That may be a bad construction way. Remove this constructor??
Strategy::Strategy(const OSISet& newOSIs, Model* newModel, const string& id):
  name("unnamed"), strategyType(id), timeDiscretisation(NULL), nsProblem(NULL), strategyxml(NULL), model(newModel),
  isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model == NULL.");
  model->setStrategyPtr(this);
  allOSI = newOSIs;
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
    isOSIAllocatedIn[*it] = false;
  //  isStrategyComplete();
}

// OSNS + Model
Strategy::Strategy(OneStepNSProblem* newNspb, Model* newModel, const string& id):
  name("unnamed"), strategyType(id), timeDiscretisation(NULL), nsProblem(newNspb), strategyxml(NULL), model(newModel),
  isTimeDiscretisationAllocatedIn(false), isNSProblemAllocatedIn(false)
{
  if (model == NULL)
    RuntimeException::selfThrow("Strategy constructor - model == NULL.");
  model->setStrategyPtr(this);
  //isStrategyComplete();
}

// --- xml constructor ---
Strategy::Strategy(StrategyXML* strxml, Model *newModel, const string& id): strategyType(id), timeDiscretisation(NULL), nsProblem(NULL),
  strategyxml(strxml), model(newModel),
  isTimeDiscretisationAllocatedIn(true), isNSProblemAllocatedIn(false)
{
  if (strategyxml != 0)
  {
    // memory allocation/construction for time discretisation
    timeDiscretisation = new TimeDiscretisation(strategyxml->getTimeDiscretisationXMLPtr(), this);

    if (model == NULL)
      RuntimeException::selfThrow("Strategy constructor - model = NULL.");
    model->setStrategyPtr(this);

    // --- OneStepIntegrators ---
    SetOfOSIXML OSIXMLList = strategyxml->getOneStepIntegratorsXML();
    SetOfOSIXMLIt it;
    CheckInsertOSI checkOSI;

    string typeOfOSI;
    for (it = OSIXMLList.begin(); it != OSIXMLList.end(); ++it)
    {
      typeOfOSI = (*it)->getType();
      // if OSI is a Moreau
      if (typeOfOSI == MOREAU_TAG)
        checkOSI = allOSI.insert(new Moreau(*it, this));

      else if (typeOfOSI == LSODAR_TAG) // if OSI is a Lsodar-type
        checkOSI = allOSI.insert(new Lsodar(*it, this));

      else RuntimeException::selfThrow("Strategy::xml constructor - unknown one-step integrator type: " + typeOfOSI);

      // checkOSI.first is an iterator that points to the OSI inserted into the set.
      isOSIAllocatedIn[*(checkOSI.first)] = true ;
    }

    // --- OneStepNSProblem ---
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
}

// --- Destructor ---
Strategy::~Strategy()
{
  if (isTimeDiscretisationAllocatedIn) delete timeDiscretisation;
  timeDiscretisation = NULL;
  if (isNSProblemAllocatedIn) delete nsProblem;
  nsProblem = NULL;

  // == delete OSI ==
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
  {
    if (isOSIAllocatedIn[*it]) delete *it;
  }

  allOSI.clear();
  isOSIAllocatedIn.clear();
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

void Strategy::setOneStepIntegrators(const OSISet& newVect)
{
  // clear old set
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
  {
    if (isOSIAllocatedIn[*it]) delete *it;
  }

  allOSI.clear();
  isOSIAllocatedIn.clear();

  // copy the new one
  allOSI = newVect;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
    isOSIAllocatedIn[*it] = false;
}

void Strategy::addOneStepIntegratorPtr(OneStepIntegrator *osi)
{
  allOSI.insert(osi);
  isOSIAllocatedIn[osi] = false;
  osi->setStrategyPtr(this);
}


void Strategy::computeFreeState()
{
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->computeFreeState();
}

void Strategy::nextStep()
{
  // increment time step
  timeDiscretisation->increment();
  model->setCurrentT(model->getCurrentT() + timeDiscretisation->getH());

  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->nextStep();

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

  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->updateState();

  // compute output (y, ydot)
  if (nsProblem != NULL) nsProblem->updateOutput();
}

void Strategy::initialize()
{
  if (model == NULL)
    RuntimeException::selfThrow("Strategy initialization - model = NULL.");

  // Initializes NSDS topology (if required)
  // (ie takes into account modifications in the NonSmoothDynamicalSystem)
  Topology * topology = getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  if (!(topology->isUpToDate()))
    topology->updateTopology();

  // initialization of the OneStepIntegrators
  OSIIterator itOsi;
  for (itOsi = allOSI.begin(); itOsi != allOSI.end(); ++itOsi)
    (*itOsi)->initialize();

  // initialization of  OneStepNonSmoothProblem
  if (nsProblem != NULL)
    nsProblem->initialize();
}

OneStepIntegrator* Strategy::getIntegratorOfDSPtr(const int& numberDS) const
{
  ConstOSIIterator itOSI  = allOSI.begin();

  //  DSIterator itDS;
  DSSet dsList = (*itOSI)->getDynamicalSystems();
  DSIterator itDS = dsList.begin();
  bool found = false;
  while (!found && itOSI != allOSI.end())
  {
    while ((*itDS)->getNumber() != numberDS && itDS != dsList.end())
      itDS++;
    if (itDS == dsList.end()) // if not found ...
    {
      itOSI++;
      dsList = (*itOSI)->getDynamicalSystems();
      itDS = dsList.begin();
    }
    else
      found = true;
  }

  return (*itOSI);


  //   vector<OneStepIntegrator*>::const_iterator it  = allOSI.begin();
  //   while((*it)->getDynamicalSystemPtr()->getNumber() != numberDS && it!= allOSI.end())
  //     it++;

  //   if(it == allOSI.end())
  //     RuntimeException::selfThrow("Strategy::getIntegratorOfDSPtr(numberDS), no integrator corresponds to this dynamical sytem");

  //   return (*it);
}

OneStepIntegrator* Strategy::getIntegratorOfDSPtr(DynamicalSystem * ds) const
{

  ConstOSIIterator itOSI  = allOSI.begin();

  ConstDSIterator itDS = (*itOSI)->getDynamicalSystems().find(ds);

  // while ds is not in the set of the current osi, scan next osi ...
  while (itDS == (*itOSI)->getDynamicalSystems().end() && itOSI != allOSI.end())
  {
    itOSI++; // go to next osi in the list
    // check if ds is present in its set
    itDS = (*itOSI)->getDynamicalSystems().find(ds);
  }

  // if ds is not find in any of the osi of the vector -> exception
  if (itOSI == allOSI.end())
    RuntimeException::selfThrow("Strategy::getIntegratorOfDSPtr(ds), no integrator corresponds to this dynamical sytem");

  return (*itOSI);
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
  if (strategyxml != NULL)
  {
    string typeOSI;
    OSIIterator it;
    for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    {
      typeOSI = (*it)->getType();
      if (typeOSI == "Moreau")
        (static_cast<Moreau*>(*it))->saveIntegratorToXML();
      else if (typeOSI == "Lsodar")
        (static_cast<Lsodar*>(*it))->saveIntegratorToXML();
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
}

bool Strategy::hasDynamicalSystemIntegrator(DynamicalSystem* ds) const
{

  ConstOSIIterator itOSI;
  ConstDSIterator itDS;
  bool val = false; // true when ds found.

  // sweep list of OSI
  for (itOSI = allOSI.begin(); itOSI != allOSI.end(); ++itOSI)
  {
    // look for ds in each osi
    itDS = ((*itOSI)->getDynamicalSystems()).find(ds);
    if (itDS != ((*itOSI)->getDynamicalSystems()).end()) // if found ...
    {
      val = true;
      break;
    }
  }
  return val;
}

bool Strategy::hasDynamicalSystemIntegrator(OneStepIntegrator* osi) const
{
  DSIterator itDS;
  bool val = false; // true when ds found.

  // sweep list of ds of osi and check that none of its ds is already present in another osi of the strategy.
  for (itDS = (osi->getDynamicalSystems()).begin(); itDS != (osi->getDynamicalSystems()).end(); ++itDS)
  {
    if (hasDynamicalSystemIntegrator(*itDS))
    {
      val = true;
      break;
    }
  }
  return val;
}
