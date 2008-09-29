/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "Model.h"
#include "EventsManager.h"

// One Step Integrators
#include "Moreau.h"
#include "Lsodar.h"

// One Step Non Smooth Problems
#include "LCP.h"
#include "QP.h"
#include "Relay.h"

using namespace std;

// Warning: neither OSI nor OSNS are given in the constructors.
// The rule is that an object OSI/OSNS needs a Simulation to be constructed. Then, the right construction is:
//   SP::Simulation s = new Simulation(...);
//   OneStepNSProblem  * osns = new OneStepNSProblem(...,s,...);
//   SP::OneStepIntegrator  osi = new OneStepIntegrator(...,s,...);

// --- Constructor with a TimeDiscretisation (and thus a Model) and an id ---
Simulation::Simulation(SP::TimeDiscretisation td, const string& id):
  name("unnamed"), simulationType(id), timeDiscretisation(td), tinit(0.0), tend(0.0), tout(0.0), model(td->getModelPtr()), levelMin(0), levelMax(0), tolerance(DEFAULT_TOLERANCE), printStat(false)
{
  if (!timeDiscretisation)
    RuntimeException::selfThrow("Simulation constructor - timeDiscretisation == NULL.");

  // Note that the link to the Model is done through the TimeDiscretisation object.
  if (!model)
    RuntimeException::selfThrow("Simulation constructor - model == NULL.");

  // === indexSets will be updated during initialize() call ===

  allOSI.reset(new OSISet());
  allNSProblems.reset(new OneStepNSProblems());
}

// --- xml constructor ---
Simulation::Simulation(SP::SimulationXML strxml, SP::Model newModel, const string& id):
  name("unnamed"), simulationType(id), tinit(0.0), tend(0.0), tout(0.0),
  simulationxml(strxml), model(newModel), levelMin(0), levelMax(0), tolerance(DEFAULT_TOLERANCE), printStat(false)
{
  if (!simulationxml)
    RuntimeException::selfThrow("Simulation:: xml constructor - xml file = NULL");

  // === Model ===
  if (!model)
    RuntimeException::selfThrow("Simulation:: xml constructor - model = NULL.");
  model->setSimulationPtr(shared_from_this());

  // === Time discretisation ===
  timeDiscretisation.reset(new TimeDiscretisation(simulationxml->getTimeDiscretisationXMLPtr(), model));

  // === OneStepIntegrators ===
  SetOfOSIXML OSIXMLList = simulationxml->getOneStepIntegratorsXML();
  SetOfOSIXMLIt it;
  CheckInsertOSI checkOSI;
  string typeOfOSI;
  allOSI.reset(new OSISet());
  for (it = OSIXMLList.begin(); it != OSIXMLList.end(); ++it)
  {
    typeOfOSI = (*it)->getType();
    // if OSI is a Moreau
    if (typeOfOSI == MOREAU_TAG)
      checkOSI = allOSI->insert(SP::Moreau(new Moreau(*it, shared_from_this())));

    else if (typeOfOSI == LSODAR_TAG) // if OSI is a Lsodar-type
      checkOSI = allOSI->insert(SP::Lsodar(new Lsodar(*it, shared_from_this())));

    else RuntimeException::selfThrow("Simulation::xml constructor - unknown one-step integrator type: " + typeOfOSI);

    // checkOSI.first is an iterator that points to the OSI inserted into the set.

  }

  // === indexSets : computed during initialization ===

  // === OneStepNSProblems  ===
  // This depends on the type of simulation --> in derived class constructor

  // === Events manager creation ===
  eventsManager.reset(new EventsManager(shared_from_this())); //
  allNSProblems.reset(new OneStepNSProblems());
}

// --- Destructor ---
Simulation::~Simulation()
{
  allNSProblems->clear();
  // -> see shared ressources for this
  if (statOut.is_open()) statOut.close();
}

// Getters/setters

void Simulation::setTimeDiscretisationPtr(SP::TimeDiscretisation td)
{
  // Warning: this function may be used carefully because of the links between Model and TimeDiscretisation
  // td will replace the current timeDiscretisation, which may change the linked Model.
  //
  if (model != td->getModelPtr())
    RuntimeException::selfThrow("Simulation::setTimeDiscretisationPtr(td) - The model associated to td is different from the one associated to the simulation.");
  timeDiscretisation = td;
}

void Simulation::setOneStepIntegrators(const OSISet& newSet)
{
  allOSI->clear();
  OSIIterator it;
  // Warning: pointers links between osi of newSet and allOSI.
  for (it = newSet.begin(); it != newSet.end(); ++it)
  {
    allOSI->insert(*it);
  }
}

SP::OneStepIntegrator Simulation::getIntegratorOfDSPtr(int numberDS) const
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

SP::OneStepIntegrator Simulation::getIntegratorOfDSPtr(SP::DynamicalSystem ds) const
{
  DSOSIConstIterator it = osiMap.find(ds);
  if (it == osiMap.end())
    RuntimeException::selfThrow("Simulation::getIntegratorOfDSPtr(ds), ds not found in the integrator set.");
  return it->second;
}

void Simulation::addOneStepIntegratorPtr(SP::OneStepIntegrator osi)
{
  allOSI->insert(osi);
  osi->setSimulationPtr(shared_from_this());

  // Note: each (ds,osi) pair will be registered into the osiMap during initialize() call (in osi->initialize).
  // During this step, we will check that each ds belongs to one and only one osi.

}

void Simulation::addInOSIMap(SP::DynamicalSystem ds, SP::OneStepIntegrator  osi)
{
  if (osiMap.find(ds) != osiMap.end()) // ie if ds is already registered in the map with another integrator
    RuntimeException::selfThrow("Simulation::addInOSIMap(ds,osi), ds is already associated with another one-step integrator");

  osiMap[ds] = osi;
}

SP::UnitaryRelationsSet Simulation::getIndexSetPtr(unsigned int i)
{
  if (i >= indexSets.size())
    RuntimeException::selfThrow("Simulation - getIndexSetPtr(i) - index set(i) does not exist.");
  return (indexSets[i]);
}

SP::OneStepNSProblem Simulation::getOneStepNSProblemPtr(const std::string& name)
{
  if (!hasOneStepNSProblem(name))
    RuntimeException::selfThrow("Simulation - getOneStepNSProblemPtr(name) - The One Step NS Problem is not in the simulation.");

  return (*allNSProblems)[name];
}

void Simulation::setOneStepNSProblems(const OneStepNSProblems& mapOfOSNS)
{
  clearOneStepNSProblems();

  // Warning: pointers links between OneStepNSProblem of each map
  allNSProblems.reset(new OneStepNSProblems());

  ConstOSNSIterator itOSNS;
  for (itOSNS = mapOfOSNS.begin(); itOSNS != mapOfOSNS.end(); ++itOSNS)
  {
    (*allNSProblems)[itOSNS->first] = itOSNS->second;
  }
}


void Simulation::clearOneStepNSProblems()
{
  allNSProblems->clear();
}

const bool Simulation::hasOneStepNSProblem(SP::OneStepNSProblem osns) const
{

  bool val = false; // true when osns found.

  ConstOSNSIterator it = allNSProblems->find(osns->getId());
  if (it != allNSProblems->end()) val = true;

  return val;
}

const bool Simulation::hasOneStepNSProblem(const string& name) const
{
  bool val = false;
  ConstOSNSIterator it;
  for (it = allNSProblems->begin(); it != allNSProblems->end(); ++it)
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

void Simulation::addOneStepNSProblemPtr(SP::OneStepNSProblem osns)
{
  if (hasOneStepNSProblem(osns))
    RuntimeException::selfThrow("Simulation - addOneStepNSProblemPtr(osns), the non smooth problem already exists in the Simulation. ");

  string name = osns->getId();
  (*allNSProblems)[name] = osns;
}

void Simulation::initialize()
{

  assert(model || !"Simulation initialization - model = NULL.");

  eventsManager.reset(new EventsManager(shared_from_this())); //
  model->setSimulationPtr(shared_from_this());

  // === Events manager initialization ===
  tinit = eventsManager->getStartingTime();
  // We first need to initialize the topology (computes UnitaryRelation sets, relative degrees ...)
  model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->initialize();

  // === OneStepIntegrators initialization ===
  OSIIterator itOsi;
  for (itOsi = allOSI->begin(); itOsi != allOSI->end(); ++itOsi)
    (*itOsi)->initialize();

  // === IndexSets building ===
  // The number of indexSets is given by the maximum value of relative degrees of the unitary relations.
  SP::InteractionsSet allInteractions = model->getNonSmoothDynamicalSystemPtr()->getInteractions();
  if (!allInteractions->isEmpty()) // ie if some Interactions have been declared
  {
    initLevelMax();

    InteractionsIterator it;
    for (it = allInteractions->begin(); it != allInteractions->end(); ++it)
      (*it)->initialize(tinit, levelMax + 1);

    indexSets.resize(levelMax + 1);
    // Link with index0 of the Topology.
    indexSets[0] = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getIndexSet0Ptr();
    for (unsigned int i = 1; i < indexSets.size(); ++i)
      indexSets[i].reset(new UnitaryRelationsSet());
  }

  // Initialize OneStepNSProblem: in derived classes specific functions.
  initOSNS();

  // Process events at time tinit. Useful to save values in memories for example.
  // Warning: can not be called during eventsManager->initialize, because it needs
  // the initialization of OSI, OSNS ...
  eventsManager->initialize();
  tend =  eventsManager->getNextTime();

  // Set Model current time (warning: current time of the model
  // corresponds to the time of the next event to be treated).
  model->setCurrentTime(getNextTime());

  // End of initialize:

  //  - all OSI and OSNS (ie DS and Interactions) states are computed for time tinit and saved into memories.
  //  - Sensors or related objects are updated for t=tinit.
  //  - current time of the model is equal to t1, time of the first event after tinit.
  //  - currentEvent of the simu. corresponds to tinit and nextEvent to tend.


  // If printStat is true, open output file.
  if (printStat)
  {
    statOut.open("simulationStat.dat", std::ofstream::out);
    statOut << "============================================" << endl;
    statOut << " Siconos Simulation of type " << simulationType << "." << endl;
    statOut << endl;
    statOut << "The tolerance parameter is equal to: " << tolerance << endl;
    statOut << endl << endl;
  }
}

void Simulation::reset()
{
  // r (or p) is set to zero in all DynamicalSystems.
  OSIIterator itOSI;
  for (itOSI = allOSI->begin(); itOSI != allOSI->end() ; ++itOSI)
    (*itOSI)->resetNonSmoothPart();
}


void Simulation::saveInMemory()
{
  // Save OSI state (DynamicalSystems) in Memory.
  OSIIterator it;
  for (it = allOSI->begin(); it != allOSI->end() ; ++it)
    (*it)->saveInMemory();

  // Save OSNS state (Interactions) in Memory.
  OSNSIterator itOsns;
  for (itOsns = allNSProblems->begin(); itOsns != allNSProblems->end(); ++itOsns)
    (itOsns->second)->saveInMemory();
}

int Simulation::computeOneStepNSProblem(const std::string& name)
{
  if (!hasOneStepNSProblem(name))
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem does not exist in the simulation. Id:" + name);
  if (!(*allNSProblems)[name])
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem == NULL, Id: " + name);

  return (*allNSProblems)[name]->compute(model->getCurrentTime());
}

void Simulation::update()
{
  for (unsigned int i = 1; i < levelMax; ++i)
    update(i);
}

void Simulation::saveSimulationToXML()
{
  if (simulationxml)
  {
    string typeOSI;
    OSIIterator it;
    for (it = allOSI->begin(); it != allOSI->end() ; ++it)
    {
      typeOSI = (*it)->getType();
      if (typeOSI == "Moreau")
        (boost::static_pointer_cast<Moreau>(*it))->saveIntegratorToXML();
      else if (typeOSI == "Lsodar")
        (boost::static_pointer_cast<Lsodar>(*it))->saveIntegratorToXML();
      else RuntimeException::selfThrow("Simulation::saveSimulationToXML - wrong type of OneStepIntegrator");
    }

    //       if( getSimulationXMLPtr()->hasOneStepNSProblemXML() )
    //  {
    //    string NSPType = nsProblem->getType();
    //    if( NSPType == "LCP" )
    //      (static_cast<LCP*>(nsProblem))->saveNSProblemToXML();
    //    else if( NSPType == "FrictionContact")
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
  // To compute input(level) (ie with lambda[level]) for all Interactions.

  if (level == -1)
    level = levelMin; // We use this since it is impossible to set levelMin as defaultValue in OneStepNSProblem.h

  //  double time = getNextTime();
  double time = model->getCurrentTime();
  SP::Topology topology = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  InteractionsIterator it;

  // Set dynamical systems non-smooth part to zero.
  reset();

  // We compute input using lambda(level).
  for (it = topology->interactionsBegin(); it != topology->interactionsEnd(); it++)
    (*it)->computeInput(time, level);
}

void Simulation::updateOutput(int level0, int level1)
{
  // To compute output() for all levels between level0 and level1 (included), for all Interactions.
  if (level1 == -1)
    level1 = levelMax;

  double time = model->getCurrentTime();
  SP::Topology topology = model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  InteractionsIterator it;

  for (it = topology->interactionsBegin(); it != topology->interactionsEnd(); it++)
  {
    for (int i = level0; i <= level1; ++i)
      (*it)->computeOutput(time , i);
  }
}

void Simulation::run(const std::string&, double, unsigned int)
{
  // Note that input arg. are useless in general case. Only useful for timeStepping.

  unsigned int count = 0; // events counter.
  cout << " ==== Start of " << simulationType << " simulation - This may take a while ... ====" << endl;
  while (getNextTime() <= model->getFinalT())
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
