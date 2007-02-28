/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "SimulationXML.h"
#include "DynamicalSystem.h"
#include "NonSmoothDynamicalSystem.h"
#include "Topology.h"
#include "Interaction.h"
#include "Relation.h"
#include "OneStepIntegratorXML.h"
#include "TimeDiscretisation.h"
#include "Model.h"
#include "EventsManager.h"

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

// Warning: neither OSI nor OSNS are given in the constructors.
// The rule is that an object OSI/OSNS needs a Simulation to be constructed. Then, the right construction is:
//   Simulation * s = new Simulation(...);
//   OneStepNSProblem  * osns = new OneStepNSProblem(...,s,...);
//   OneStepIntegrator * osi = new OneStepIntegrator(...,s,...);

// --- Default constructor (protected) ---
Simulation::Simulation(const string& type):
  name("unnamed"), simulationType(type), timeDiscretisation(NULL), eventsManager(NULL),
  simulationxml(NULL), model(NULL), levelMin(0), levelMax(0)
{
  isAllocatedIn["eventsManager"] = false;
  isAllocatedIn["timeDiscretisation"] = false;
}

// --- Constructor with a TimeDiscretisation (and thus a Model) and an id ---
Simulation::Simulation(TimeDiscretisation* td, const string& id):
  name("unnamed"), simulationType(id), timeDiscretisation(td), eventsManager(NULL),
  simulationxml(NULL), model(td->getModelPtr()), levelMin(0), levelMax(0)
{
  if (timeDiscretisation == NULL)
    RuntimeException::selfThrow("Simulation constructor - timeDiscretisation == NULL.");
  isAllocatedIn["timeDiscretisation"] = false;

  // Note that the link to the Model is done through the TimeDiscretisation object.
  if (model == NULL)
    RuntimeException::selfThrow("Simulation constructor - model == NULL.");

  model->setSimulationPtr(this);
  // === indexSets will be updated during initialize() call ===

  // === Events manager creation ===
  eventsManager = new EventsManager(this); //
  isAllocatedIn["eventsManager"] = true;
}

// --- xml constructor ---
Simulation::Simulation(SimulationXML* strxml, Model *newModel, const string& id):
  name("unnamed"), simulationType(id), timeDiscretisation(NULL), eventsManager(NULL),
  simulationxml(strxml), model(newModel), levelMin(0), levelMax(0)
{
  if (simulationxml == NULL)
    RuntimeException::selfThrow("Simulation:: xml constructor - xml file = NULL");

  // === Model ===
  if (model == NULL)
    RuntimeException::selfThrow("Simulation:: xml constructor - model = NULL.");
  model->setSimulationPtr(this);

  // === Time discretisation ===
  timeDiscretisation = new TimeDiscretisation(simulationxml->getTimeDiscretisationXMLPtr(), model);
  isAllocatedIn["timeDiscretisation"] = true;

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
  // This depends on the type of simulation --> in derived class constructor

  // === Events manager creation ===
  eventsManager = new EventsManager(this); //
  isAllocatedIn["eventsManager"] = true;
}

// --- Destructor ---
Simulation::~Simulation()
{
  // == EventsManager ==
  if (isAllocatedIn["eventsManager"])
    delete eventsManager;
  eventsManager = NULL;

  // == Time discretisation ==
  if (isAllocatedIn["timeDiscretisation"]) delete timeDiscretisation;
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
  if (isAllocatedIn["timeDiscretisation"]) delete timeDiscretisation;
  timeDiscretisation = td;
  isAllocatedIn["timeDiscretisation"] = false;
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

OneStepIntegrator* Simulation::getIntegratorOfDSPtr(int numberDS) const
{

  DSOSIConstIterator it = osiMap.begin();
  bool found = false;

  while (!found || it != osiMap.end())
  {
    if ((it->first)->getNumber() == numberDS)
      found = true;
    else ++it;
  }

  return (it->second);
}

OneStepIntegrator* Simulation::getIntegratorOfDSPtr(DynamicalSystem * ds) const
{
  return (osiMap.find(ds))->second;
}

void Simulation::addOneStepIntegratorPtr(OneStepIntegrator *osi)
{
  if (hasCommonDSInIntegrators(osi))
    RuntimeException::selfThrow("Simulation::addOneStepIntegratorPtr(osi), one of the dynamical system of the new integrator has already an integrator defined in the simulation!");

  allOSI.insert(osi);
  isOSIAllocatedIn[osi] = false;
  osi->setSimulationPtr(this);
}

void Simulation::addInOSIMap(DynamicalSystem * ds, OneStepIntegrator * osi)
{
  osiMap[ds] = osi;
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

const UnitaryRelationsSet Simulation::getIndexSet(unsigned int i) const
{
  if (i >= indexSets.size())
    RuntimeException::selfThrow("Simulation - getIndexSet(i) - index set(i) does not exist.");
  return indexSets[i];
}

OneStepNSProblem* Simulation::getOneStepNSProblemPtr(const std::string& name)
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

const bool Simulation::hasOneStepNSProblem(const string& name) const
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

  // Initialize the user time discretisation.
  timeDiscretisation->initialize();

  // === Events manager initialization ===
  eventsManager->initialize();

  // We first need to initialize the topology (computes UnitaryRelation sets, relative degrees ...)
  model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->initialize();

  // initialization of the OneStepIntegrators
  OSIIterator itOsi;
  for (itOsi = allOSI.begin(); itOsi != allOSI.end(); ++itOsi)
    (*itOsi)->initialize();
  // The number of indexSets is given by the maximum value of relative degrees of the unitary relations.

  InteractionsSet allInteractions = model->getNonSmoothDynamicalSystemPtr()->getInteractions();
  double t0 = model->getT0();
  if (!allInteractions.isEmpty()) // ie if some Interactions have been declared
  {
    initLevelMax();

    InteractionsIterator it;
    for (it = allInteractions.begin(); it != allInteractions.end(); ++it)
      (*it)->initialize(t0, levelMax + 1);

    indexSets.resize(levelMax + 1);
    indexSets[0] = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getIndexSet0();
  }

  // Initialize OneStepNSProblem: in derived classes specific function.
  initOSNS();

  // == Call process functions of events (usefull to save initial values for example) ==
  eventsManager->process();
  // Set Model current time (warning: current time of the model corresponds to the time of the next event to be treated).
  model->setCurrentT(getNextTime());

  // End of initialize:
  //  - all OSI and OSNS (ie DS and Interactions) states are computed for time t0 and saved into memories.
  //  - Sensors or related objects are updated for t=t0.
  //  - current time of the model is equal to t1, time of the first event after t0.
  //  - currentEvent of the simu. corresponds to t0 and nextEvent to t1.
}

void Simulation::reset()
{
  // r (or p) is set to zero in all DynamicalSystems.
  OSIIterator itOSI;
  for (itOSI = allOSI.begin(); itOSI != allOSI.end() ; ++itOSI)
    (*itOSI)->resetNonSmoothPart();
}

void Simulation::saveInMemory()
{
  // Save OSI state (DynamicalSystems) in Memory.
  OSIIterator it;
  for (it = allOSI.begin(); it != allOSI.end() ; ++it)
    (*it)->saveInMemory();

  // Save OSNS state (Interactions) in Memory.
  OSNSIterator itOsns;
  for (itOsns = allNSProblems.begin(); itOsns != allNSProblems.end(); ++itOsns)
    (itOsns->second)->saveInMemory();
}

void Simulation::computeOneStepNSProblem(const std::string& name)
{
  if (!hasOneStepNSProblem(name))
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem does not exist in the simulation. Id:" + name);
  if (allNSProblems[name] == NULL)
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem == NULL, Id: " + name);

  allNSProblems[name]->compute(model->getCurrentT());
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

  //  double time = getNextTime();
  double time = model->getCurrentT();
  InteractionsSet allInter = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getInteractions();
  InteractionsIterator it;

  reset();

  // We compute input using lambda(levelMin).
  for (it = allInter.begin(); it != allInter.end(); it++)
    (*it)->computeInput(time, level);
}

void Simulation::updateOutput(int level0, int level1)
{

  if (level1 == -1)
    level1 = levelMax;

  double time = model->getCurrentT();
  InteractionsSet allInter = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getInteractions();
  InteractionsIterator it;

  for (it = allInter.begin(); it != allInter.end(); it++)
  {
    for (int i = level0; i <= level1; ++i)
      (*it)->computeOutput(time , i);
  }
}

void Simulation::run(const std::string&, double, unsigned int)
{
  // Note that input arg. are useless in general case. Only useful for timeStepping.

  unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of " << simulationType << " simulation - This may take a while ... ====" << endl;
  while (eventsManager->hasNextEvent())
  {
    advanceToEvent();
    eventsManager->processEvents();
    count++;
  }
  cout << "===== End of " << simulationType << "simulation. " << count << " events have been processed. ==== " << endl;
}

void Simulation::processEvents()
{
  eventsManager->processEvents();
}
