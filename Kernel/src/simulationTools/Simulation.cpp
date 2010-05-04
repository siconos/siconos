/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "Simulation.hpp"
#include "SimulationXML.hpp"
#include "DynamicalSystem.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"
#include "OneStepIntegratorXML.hpp"
#include "EventsManager.hpp"

// One Step Integrators
#include "Moreau.hpp"
#include "Lsodar.hpp"

// One Step Non Smooth Problems
#include "LCP.hpp"
#include "QP.hpp"
#include "Relay.hpp"

using namespace std;

// --- Constructor with a TimeDiscretisation (and thus a Model) and an
// --- id ---
Simulation::Simulation(SP::TimeDiscretisation td):
  _name("unnamed"), _timeDiscretisation(td),
  _tinit(0.0), _tend(0.0), _tout(0.0), _levelMin(0), _levelMax(0),
  _tolerance(DEFAULT_TOLERANCE), _printStat(false)
{
  if (!_timeDiscretisation)
    RuntimeException::selfThrow("Simulation constructor - timeDiscretisation == NULL.");
  mUseRelativeConvergenceCriterion = false;
  mRelativeConvergenceCriterionHeld = false;
  mRelativeConvergenceTol = 10e-3;

  // === indexSets will be updated during initialize() call ===

  _allOSI.reset(new OSISet());
  _allNSProblems.reset(new OneStepNSProblems());
  _eventsManager.reset(new EventsManager()); //
}

// --- xml constructor ---
Simulation::Simulation(SP::SimulationXML strxml, double t0, double T, SP::DynamicalSystemsSet dsList,
                       SP::InteractionsSet interactionsList):
  _name("unnamed"), _tinit(0.0), _tend(0.0), _tout(0.0),
  _simulationxml(strxml), _levelMin(0), _levelMax(0), _tolerance(DEFAULT_TOLERANCE), _printStat(false)
{
  if (!_simulationxml)
    RuntimeException::selfThrow("Simulation:: xml constructor - xml file = NULL");
  mUseRelativeConvergenceCriterion = false;
  mRelativeConvergenceCriterionHeld = false;
  mRelativeConvergenceTol = 10e-3;


  // === Model ===

  // === Time discretisation ===
  _timeDiscretisation.reset(new TimeDiscretisation(_simulationxml->timeDiscretisationXML(), t0, T));

  // === OneStepIntegrators ===
  SetOfOSIXML OSIXMLList = _simulationxml->getOneStepIntegratorsXML();
  SetOfOSIXMLIt it;
  string typeOfOSI;
  _allOSI.reset(new OSISet());
  for (it = OSIXMLList.begin(); it != OSIXMLList.end(); ++it)
  {
    typeOfOSI = (*it)->getType();
    // if OSI is a Moreau
    if (typeOfOSI == MOREAU_TAG)
      _allOSI->insert(SP::Moreau(new Moreau(*it, dsList)));

    else if (typeOfOSI == LSODAR_TAG) // if OSI is a Lsodar-type
      _allOSI->insert(SP::Lsodar(new Lsodar(*it, dsList, interactionsList)));

    else RuntimeException::selfThrow("Simulation::xml constructor - unknown one-step integrator type: " + typeOfOSI);
  }

  // === indexSets : computed during initialization ===

  // === OneStepNSProblems  ===
  // This depends on the type of simulation --> in derived class constructor

  // === Events manager creation ===
  _eventsManager.reset(new EventsManager()); //
  _allNSProblems.reset(new OneStepNSProblems());
}

// --- Destructor ---
Simulation::~Simulation()
{
  _allNSProblems->clear();
  _allOSI->clear();
  _osiMap.clear();

  _allNSProblems->clear();
  // -> see shared ressources for this
  if (statOut.is_open()) statOut.close();
}

// Getters/setters

void Simulation::setOneStepIntegrators(const OSISet& newSet)
{
  _allOSI->clear();
  _allOSI->insert(newSet.begin(), newSet.end());
}

SP::OneStepIntegrator Simulation::integratorOfDS(int numberDS) const
{

  DSOSIConstIterator it = _osiMap.begin();
  bool found = false;

  while (!found || it != _osiMap.end())
  {
    if ((it->first)->number() == numberDS)
      found = true;
    else ++it;
  }

  return (it->second);
}

SP::OneStepIntegrator Simulation::integratorOfDS(SP::DynamicalSystem ds) const
{
  DSOSIConstIterator it = _osiMap.find(ds);
  if (it == _osiMap.end())
    RuntimeException::selfThrow("Simulation::integratorOfDS(ds), ds not found in the integrator set.");
  return it->second;
}

void Simulation::insertIntegrator(SP::OneStepIntegrator osi)
{
  _allOSI->insert(osi);
  // Note: each (ds,osi) pair will be registered into the _osiMap
  // during initialize() call (in osi->initialize).  During this step,
  // we will check that each ds belongs to one and only one osi.
}

void Simulation::addInOSIMap(SP::DynamicalSystem ds, SP::OneStepIntegrator  osi)
{
  if (_osiMap.find(ds) != _osiMap.end()) // ie if ds is already registered
    // in the map with another
    // integrator
    ;/*RuntimeException::selfThrow("Simulation::addInOSIMap(ds,osi), ds is already associated with another one-step integrator");  */
  _osiMap[ds] = osi;
}


SP::OneStepNSProblem Simulation::oneStepNSProblem(int Id)
{
  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - oneStepNSProblem(Id) - The One Step NS Problem is not in the simulation.");

  return (*_allNSProblems)[Id];
}
/*
void Simulation::setOneStepNSProblems(const OneStepNSProblems& mapOfOSNS)
{
  clearOneStepNSProblems();

  // Warning: pointers links between OneStepNSProblem of each map
  _allNSProblems.reset(new OneStepNSProblems());
  for(ConstOSNSIterator itOSNS = mapOfOSNS.begin(); itOSNS != mapOfOSNS.end(); ++itOSNS)
    (*_allNSProblems)[itOSNS->first] = itOSNS->second;
}


void Simulation::clearOneStepNSProblems()
{
  _allNSProblems->clear();
}

const bool Simulation::hasOneStepNSProblem(SP::OneStepNSProblem osns) const
{

  bool val = false; // true when osns found.

  ConstOSNSIterator it = _allNSProblems->find(osns->getId());
  if (it!= _allNSProblems->end() ) val = true;

  return val;
}

const bool Simulation::hasOneStepNSProblem(const string& name) const
{
  bool val = false;
  ConstOSNSIterator it;
  for(it = _allNSProblems->begin(); it!= _allNSProblems->end(); ++it)
    if( (it->second)->getId() == name)
      {
  val = true;
  break;
      }
  return val;
}
*/
void Simulation::updateIndexSets()
{
  // Warning, I0 is not updated and must remain unchanged !
  unsigned int nindexsets = model()->nonSmoothDynamicalSystem()
                            ->topology()->indexSetsSize();

  if (nindexsets > 1)
  {
    for (unsigned int i = 1; i < nindexsets ; ++i)
      updateIndexSet(i);
  }
}

void Simulation::insertNonSmoothProblem(SP::OneStepNSProblem osns, int Id)
{
  if (((*_allNSProblems)[Id]))
    RuntimeException::selfThrow("Simulation - insertNonSmoothProblem(osns), trying to insert a OSNSP already existing. ");
  (*_allNSProblems)[Id] = osns;

}

void Simulation::updateInteractions()
{

  SP::InteractionsSet allInteractions =
    model()->nonSmoothDynamicalSystem()->interactions();
  initLevelMax();

  double time = model()->currentTime(); // init with current model time

  std::for_each(allInteractions->begin(), allInteractions->end(),
                boost::bind(&Interaction::initialize, _1, time, _levelMax + 1));

  initOSNS();

}

void Simulation::initialize(SP::Model m, bool withOSI)
{
  // === Connection with the model ===
  assert(m || !"Simulation::initialize(model) - model = NULL.");
  _model = boost::weak_ptr<Model>(m);

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  // === Events manager initialization ===
  _eventsManager->initialize(shared_from_this());
  _tinit = _eventsManager->startingTime();

  if (withOSI)
  {
    // === OneStepIntegrators initialization ===

    for (OSIIterator itosi = _allOSI->begin();
         itosi != _allOSI->end(); ++itosi)
    {


      for (DSIterator itds = (*itosi)->dynamicalSystems()->begin();
           itds != (*itosi)->dynamicalSystems()->end();
           ++itds)
      {
        (*itds)->initialize(Type::name(*shared_from_this()), model()->t0(),
                            (*itosi)->getSizeMem());
        addInOSIMap(*itds, *itosi);
      }

      (*itosi)->setSimulationPtr(shared_from_this());
      (*itosi)->initialize();

    }
  }

  // === IndexSets building ===

  // The number of indexSets is given by the maximum value of relative
  // degrees of the unitary relations.
  SP::InteractionsSet allInteractions =
    model()->nonSmoothDynamicalSystem()->interactions();
  //  if( !allInteractions->isEmpty() ) // ie if some Interactions
  //  have been declared
  {
    initLevelMax();
    topo->indexSetsResize(_levelMax + 1);

    std::for_each(allInteractions->begin(), allInteractions->end(),
                  boost::bind(&Interaction::initialize, _1, _tinit, _levelMax + 1));

    for (unsigned int i = 1; i < topo->indexSetsSize(); ++i)
      topo->resetIndexSetPtr(i);

    // Initialize OneStepNSProblem: in derived classes specific
    // functions.
    initOSNS();
  }


  // Process events at time _tinit. Useful to save values in memories
  // for example.  Warning: can not be called during
  // eventsManager->initialize, because it needs the initialization of
  // OSI, OSNS ...
  _eventsManager->preUpdate();
  _tend =  _eventsManager->nextTime();

  // Set Model current time (warning: current time of the model
  // corresponds to the time of the next event to be treated).
  model()->setCurrentTime(nextTime());

  // End of initialize:

  //  - all OSI and OSNS (ie DS and Interactions) states are computed
  //  - for time _tinit and saved into memories.
  //  - Sensors or related objects are updated for t=_tinit.
  //  - current time of the model is equal to t1, time of the first
  //  - event after _tinit.
  //  - currentEvent of the simu. corresponds to _tinit and nextEvent
  //  - to _tend.

  // If _printStat is true, open output file.
  if (_printStat)
  {
    statOut.open("simulationStat.dat", std::ofstream::out);
    statOut << "============================================" << endl;
    statOut << " Siconos Simulation of type " << typeName() << "." << endl;
    statOut << endl;
    statOut << "The tolerance parameter is equal to: " << _tolerance << endl;
    statOut << endl << endl;
  }
}

void Simulation::reset()
{
  // r (or p) is set to zero in all DynamicalSystems.
  OSIIterator itOSI;
  for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    (*itOSI)->resetNonSmoothPart();
}


void Simulation::saveInMemory()
{
  // Save OSI state (DynamicalSystems) in Memory.
  OSIIterator it;
  for (it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    (*it)->saveInMemory();

  // Save OSNS state (Interactions) in Memory.
  OSNSIterator itOsns;
  for (itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
    (*itOsns)->saveInMemory();
}

int Simulation::computeOneStepNSProblem(int Id)
{

  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem == NULL, Id: " + Id);

  return (*_allNSProblems)[Id]->compute(model()->currentTime());
}

void Simulation::update()
{
  for (unsigned int i = 1; i < _levelMax; ++i)
    update(i);
}

void Simulation::saveSimulationToXML()
{
  if (_simulationxml)
  {
    OSI::TYPES typeOSI;
    OSIIterator it;
    for (it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    {
      typeOSI = (*it)->getType();
      if (typeOSI == OSI::MOREAU)
        (boost::static_pointer_cast<Moreau>(*it))->saveIntegratorToXML();
      else if (typeOSI == OSI::LSODAR)
        (boost::static_pointer_cast<Lsodar>(*it))->saveIntegratorToXML();
      else RuntimeException::selfThrow("Simulation::saveSimulationToXML - wrong type of OneStepIntegrator");
    }

  }
  else RuntimeException::selfThrow("Simulation::saveSimulationToXML - SimulationXML = NULL");
}

void Simulation::updateInput(int level)
{
  // To compute input(level) (ie with lambda[level]) for all Interactions.

  if (level == -1)
    level = _levelMin; // We use this since it is impossible to set
  // _levelMin as defaultValue in
  // OneStepNSProblem.h

  //  double time = nextTime();
  double time = model()->currentTime();
  SP::Topology topology = model()->nonSmoothDynamicalSystem()->topology();
  InteractionsIterator it;

  // Set dynamical systems non-smooth part to zero.
  reset();

  // We compute input using lambda(level).
  for (it = topology->interactions()->begin();
       it != topology->interactions()->end(); it++)
    (*it)->computeInput(time, level);
}

void Simulation::updateOutput(int level0, int level1)
{
  // To compute output() for all levels between level0 and level1
  // (included), for all Interactions.
  if (level1 == -1)
    level1 = _levelMax;

  double time = model()->currentTime();
  SP::Topology topology = model()->nonSmoothDynamicalSystem()->topology();
  InteractionsIterator it;

  for (it = topology->interactions()->begin();
       it != topology->interactions()->end(); it++)
  {
    for (int i = level0; i <= level1; ++i)
      (*it)->computeOutput(time , i);
  }
}

void Simulation::run(const std::string&, double, unsigned int)
{
  // Note that input arg. are useless in general case. Only useful for
  // timeStepping.

  unsigned int count = 0; // events counter.
  cout << " ==== Start of " << typeName() << " simulation - This may take a while ... ====" << endl;
  while (nextTime() <= model()->finalT())
  {
    advanceToEvent();
    _eventsManager->processEvents();
    count++;
  }
  cout << "===== End of " << typeName() << "simulation. " << count << " events have been processed. ==== " << endl;
}

void Simulation::processEvents()
{
  _eventsManager->processEvents();

  /* should be evaluated only if needed */
  SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    dsGraph->bundle(*vi)->endStep();
  }
}

