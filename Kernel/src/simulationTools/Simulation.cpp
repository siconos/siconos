/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
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
#include "Simulation.h"
// includes to be deleted thanks to factories:

// One Step Integrators
#include "Moreau.h"
#include "Lsodar.h"

// One Step Non Smooth Problems
#include "LCP.h"
#include "FrictionContact2D.h"
#include "FrictionContact3D.h"
#include "QP.h"
#include "Relay.h"

using namespace std;

// Warning: neither time discretisation nor OSI or OSNS are given in the constructors.
// The rule is that an object TimeDiscretisation/OSI/OSNS needs a Simulation to be constructed. Then, the right construction is:
//   Simulation * s = new Simulation(...);
//   TimeDiscretisation* t = new TimeDiscretisation(...,s,...);
//   OneStepNSProblem  * osns = new OneStepNSProblem(...,s,...);
//   OneStepIntegrator * osi = new OneStepIntegrator(...,s,...);

// --- Default constructor (protected) ---
Simulation::Simulation(const string type):
  name("unnamed"), simulationType(type), timeDiscretisation(NULL), isTimeDiscretisationAllocatedIn(false),
  simulationxml(NULL), model(NULL), levelMin(0), levelMax(0)
{}

// --- constructor with a Model and an id ---
Simulation::Simulation(Model * mainModel, const string id):
  name("unnamed"), simulationType(id), timeDiscretisation(NULL), isTimeDiscretisationAllocatedIn(false),
  simulationxml(NULL), model(mainModel), levelMin(0), levelMax(0)
{
  if (model == NULL)
    RuntimeException::selfThrow("Simulation constructor - model == NULL.");
  model->setSimulationPtr(this);
  // === indexSets will be updated during initialize() call ===
}

// --- xml constructor ---
Simulation::Simulation(SimulationXML* strxml, Model *newModel, const string id):
  name("unnamed"), simulationType(id), timeDiscretisation(NULL), isTimeDiscretisationAllocatedIn(true),
  simulationxml(strxml), model(newModel), levelMin(0), levelMax(0)
{
  if (simulationxml == NULL)
    RuntimeException::selfThrow("Simulation:: xml constructor - xml file = NULL");

  // === Model ===
  if (model == NULL)
    RuntimeException::selfThrow("Simulation constructor - model = NULL.");
  model->setSimulationPtr(this);

  // === Time discretisation ===
  timeDiscretisation = new TimeDiscretisation(simulationxml->getTimeDiscretisationXMLPtr(), this);

  // === OneStepIntegrators ===
  SetOfOSIXML OSIXMLList = simulationxml->getOneStepIntegratorsXML();
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

    else RuntimeException::selfThrow("Simulation::xml constructor - unknown one-step integrator type: " + typeOfOSI);

    // checkOSI.first is an iterator that points to the OSI inserted into the set.
    isOSIAllocatedIn[*(checkOSI.first)] = true ;
  }

  // === indexSets : computed during initialization ===

  // === OneStepNSProblems  ===
  // This depends on the type of strategy --> in derived class constructor
}

// --- Destructor ---
Simulation::~Simulation()
{
  if (isTimeDiscretisationAllocatedIn) delete timeDiscretisation;
  timeDiscretisation = NULL;

  // == delete OSI ==
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end(); ++it)
  {
    if (isOSIAllocatedIn[*it]) delete *it;
  }

  allOSI.clear();
  isOSIAllocatedIn.clear();

  // == delete OS NS Problems ==
  OSNSIterator itOSNS;
  for (itOSNS = allNSProblems.begin(); itOSNS != allNSProblems.end(); ++itOSNS)
  {
    if (isNSProblemAllocatedIn[itOSNS->second]) delete(itOSNS->second);
  }

  allNSProblems.clear();
  isNSProblemAllocatedIn.clear();

  indexSets.clear();
}

// Getters/setters

void Simulation::setTimeDiscretisationPtr(TimeDiscretisation* td)
{
  // Warning: this function may be used carefully because of the links between Model and TimeDiscretisation
  // The simulation of the td input MUST be the current simulation.
  //
  if (isTimeDiscretisationAllocatedIn) delete timeDiscretisation;
  timeDiscretisation = td;
  isTimeDiscretisationAllocatedIn = false;
}

void Simulation::setOneStepIntegrators(const OSISet& newVect)
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

OneStepIntegrator* Simulation::getIntegratorOfDSPtr(const int numberDS) const
{
  ConstOSIIterator itOSI  = allOSI.begin();

  //  DSIterator itDS;
  DynamicalSystemsSet dsList = (*itOSI)->getDynamicalSystems();
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
}

OneStepIntegrator* Simulation::getIntegratorOfDSPtr(DynamicalSystem * ds) const
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
    RuntimeException::selfThrow("Simulation::getIntegratorOfDSPtr(ds), no integrator corresponds to this dynamical sytem");

  return (*itOSI);
}

void Simulation::addOneStepIntegratorPtr(OneStepIntegrator *osi)
{
  if (hasCommonDSInIntegrators(osi))
    RuntimeException::selfThrow("Simulation::addOneStepIntegratorPtr(osi), one of the dynamical system of the new integrator has already an integrator defined in the simulation!");

  allOSI.insert(osi);
  isOSIAllocatedIn[osi] = false;
  osi->setSimulationPtr(this);
}

const bool Simulation::hasCommonDSInIntegrators(OneStepIntegrator* osi) const
{
  bool val = false;
  ConstOSIIterator itOSI;
  DynamicalSystemsSet commonDS;

  for (itOSI = allOSI.begin(); itOSI != allOSI.begin() ; ++itOSI)
  {
    // check intersection of set of dynamical system of new osi and all osi of the simulation
    val = (intersection((*itOSI)->getDynamicalSystems() , osi->getDynamicalSystems())).isEmpty();
    if (!val) break;
  }

  return val;
}

const UnitaryRelationsSet Simulation::getIndexSet(const unsigned int i) const
{
  if (i >= indexSets.size())
    RuntimeException::selfThrow("Simulation - getIndexSet(i) - index set(i) does not exist.");
  return indexSets[i];
}

OneStepNSProblem* Simulation::getOneStepNSProblemPtr(const std::string name)
{
  if (!hasOneStepNSProblem(name))
    RuntimeException::selfThrow("Simulation - getOneStepNSProblemPtr(name) - The One Step NS Problem is not in the simulation.");

  return allNSProblems[name];
}

void Simulation::setOneStepNSProblems(const OneStepNSProblems& mapOfOSNS)
{
  clearOneStepNSProblems();

  // Warning: pointers links between OneStepNSProblem of each map
  allNSProblems = mapOfOSNS;
  OSNSIterator itOSNS;
  for (itOSNS = allNSProblems.begin(); itOSNS != allNSProblems.end(); ++itOSNS)
    isNSProblemAllocatedIn[itOSNS->second] = false;
}

void Simulation::clearOneStepNSProblems()
{
  OSNSIterator itOSNS;
  for (itOSNS = allNSProblems.begin(); itOSNS != allNSProblems.end(); ++itOSNS)
  {
    if (isNSProblemAllocatedIn[itOSNS->second] && itOSNS->second != NULL)
      delete(itOSNS->second);
  }
  allNSProblems.clear();
  isNSProblemAllocatedIn.clear();
}

const bool Simulation::hasOneStepNSProblem(OneStepNSProblem* osns) const
{

  bool val = false; // true when osns found.

  ConstOSNSIterator it = allNSProblems.find(osns->getId());
  if (it != allNSProblems.end()) val = true;

  return val;
}

const bool Simulation::hasOneStepNSProblem(const string name) const
{
  bool val = false;
  ConstOSNSIterator it;
  for (it = allNSProblems.begin(); it != allNSProblems.end(); ++it)
    if ((it->second)->getId() == name)
    {
      val = true;
      break;
    }
  return val;
}

void Simulation::updateIndexSets()
{
  // Warning, I0 is not updated and must remain unchanged !
  if (indexSets.size() > 1)
  {
    for (unsigned int i = 1; i < indexSets.size() ; ++i)
      updateIndexSet(i);
  }
}

void Simulation::addOneStepNSProblemPtr(OneStepNSProblem* osns)
{
  if (hasOneStepNSProblem(osns))
    RuntimeException::selfThrow("Simulation - addOneStepNSProblemPtr(osns), the non smooth problem already exists in the Simulation. ");

  string name = osns->getId();
  allNSProblems[name] = osns;
  isNSProblemAllocatedIn[osns] = false;
}

void Simulation::initialize()
{
  if (model == NULL)
    RuntimeException::selfThrow("Simulation initialization - model = NULL.");

  // The number of indexSets is given by the maximum value of relative degrees of the unitary relations.

  // We first need to initialize the topology (computes UnitaryRelation sets, relative degrees ...)
  model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->initialize();
}

void Simulation::computeFreeState()
{
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->computeFreeState();
}

void Simulation::nextStep()
{
  // increment time step
  timeDiscretisation->increment();
  model->setCurrentT(model->getCurrentT() + timeDiscretisation->getH());

  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->nextStep();

  OSNSIterator itOsns;
  for (itOsns = allNSProblems.begin(); itOsns != allNSProblems.end(); ++itOsns)
    (itOsns->second)->nextStep();
}

void Simulation::computeOneStepNSProblem(const std::string name)
{
  if (!hasOneStepNSProblem(name))
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem does not exist in the simulation. Id:" + name);
  if (allNSProblems[name] == NULL)
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem == NULL, Id: " + name);

  allNSProblems[name]->compute(model->getCurrentT());
}

void Simulation::newtonSolve(const double criterion, const unsigned int maxStep)
{
  // At the time, only for time stepping
  if (simulationType != "TimeStepping")
    RuntimeException::selfThrow("Simulation::newtonSolve - Not yet implemented for simulation of type" + simulationType);

  bool isNewtonConverge = false;
  unsigned int nbNewtonStep = 0; // number of Newton iterations
  while ((!isNewtonConverge) && (nbNewtonStep <= maxStep))
  {
    nbNewtonStep++;
    computeFreeState();
    updateIndexSets();
    computeOneStepNSProblem("timeStepping");
    update(levelMin);
    isNewtonConverge = newtonCheckConvergence(criterion);
  }
  if (!isNewtonConverge)
    cout << "Newton process stopped: reach max step number" << endl ;

  // time step increment
  model->setCurrentT(model->getCurrentT() + timeDiscretisation->getH());
}

bool Simulation::newtonCheckConvergence(const double criterion)
{
  bool checkConvergence = false;
  // get the non smooth dynamical system
  NonSmoothDynamicalSystem* nsds = model-> getNonSmoothDynamicalSystemPtr();
  // get the nsds indicator of convergence
  double nsdsConverge = nsds -> nsdsConvergenceIndicator();
  if (nsdsConverge < criterion) checkConvergence = true ;

  return(checkConvergence);
}

void Simulation::saveSimulationToXML()
{
  if (simulationxml != NULL)
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
      else RuntimeException::selfThrow("Simulation::saveSimulationToXML - wrong type of OneStepIntegrator");
    }

    //       if( getSimulationXMLPtr()->hasOneStepNSProblemXML() )
    //  {
    //    string NSPType = nsProblem->getType();
    //    if( NSPType == "LCP" )
    //      (static_cast<LCP*>(nsProblem))->saveNSProblemToXML();
    //    else if( NSPType == "FrictionContact2D" || NSPType == "FrictionContact3D")
    //      (static_cast<FrictionContact*>(nsProblem))->saveNSProblemToXML();
    //    else if( NSPType == "QP")
    //      (static_cast<QP*>(nsProblem))->saveNSProblemToXML();
    //    else if( NSPType == "Relay" )
    //      (static_cast<Relay*>(nsProblem))->saveNSProblemToXML();
    //    else
    //      RuntimeException::selfThrow("Simulation::saveSimulationToXML - wrong type of OneStepNSProblem");
    //  }
  }
  else RuntimeException::selfThrow("Simulation::saveSimulationToXML - SimulationXML = NULL");
}

void Simulation::updateInput(int level)
{
  if (level == -1)
    level = levelMin; // We use this since it is impossible to set levelMin as defaultValue in OneStepNSProblem.h

  double time = model->getCurrentT();
  InteractionsSet allInter = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getInteractions();
  InteractionsIterator it;

  // First, r (or p) is set to zero in all DynamicalSystems.
  OSIIterator itOSI;
  for (itOSI = allOSI.begin(); itOSI != allOSI.end() ; ++itOSI)
    (*itOSI)->resetNonSmoothPart();

  // We compute inpute using lambda(levelMin).
  for (it = allInter.begin(); it != allInter.end(); it++)
    (*it)->getRelationPtr() -> computeInput(time, level);
}

void Simulation::updateOutput(const int level0, int level1)
{

  if (level1 == -1)
    level1 = levelMax;

  double time = model->getCurrentT();
  InteractionsSet allInter = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getInteractions();
  InteractionsIterator it;

  for (it = allInter.begin(); it != allInter.end(); it++)
  {
    for (int i = level0; i <= level1; ++i)
      (*it)->getRelationPtr()->computeOutput(time , i);
  }
}


