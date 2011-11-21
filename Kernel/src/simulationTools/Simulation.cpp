/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "LagrangianDS.hpp"

// One Step Integrators
#include "Moreau.hpp"
#include "Lsodar.hpp"
#include "D1MinusLinear.hpp"
#include "SchatzmanPaoli.hpp"

// One Step Non Smooth Problems
#include "LCP.hpp"
#include "QP.hpp"
#include "Relay.hpp"

using namespace std;


// --- Constructor with a TimeDiscretisation (and thus a Model) and an
// --- id ---
Simulation::Simulation(SP::TimeDiscretisation td):
  _name("unnamed"), _timeDiscretisation(td),
  _tinit(0.0), _tend(0.0), _tout(0.0),
  _tolerance(DEFAULT_TOLERANCE), _printStat(false), _staticLevels(false), _levelsAreComputed(false)
{
  if (!_timeDiscretisation)
    RuntimeException::selfThrow("Simulation constructor - timeDiscretisation == NULL.");
  _useRelativeConvergenceCriterion = false;
  _relativeConvergenceCriterionHeld = false;
  _relativeConvergenceTol = 10e-3;

  // === indexSets will be updated during initialize() call ===

  _allOSI.reset(new OSISet());
  _allNSProblems.reset(new OneStepNSProblems());
  _eventsManager.reset(new EventsManager()); //
}

// --- xml constructor ---
Simulation::Simulation(SP::SimulationXML strxml, double t0, double T, SP::DynamicalSystemsSet dsList,
                       SP::InteractionsSet interactionsList):
  _name("unnamed"), _tinit(0.0), _tend(0.0), _tout(0.0),
  _simulationxml(strxml),
  _tolerance(DEFAULT_TOLERANCE), _printStat(false), _staticLevels(false), _levelsAreComputed(false)
{
  if (!_simulationxml)
    RuntimeException::selfThrow("Simulation:: xml constructor - xml file = NULL");
  _useRelativeConvergenceCriterion = false;
  _relativeConvergenceCriterionHeld = false;
  _relativeConvergenceTol = 10e-3;


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
  _interactionOsiMap.clear();

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

SP::OneStepIntegrator Simulation::integratorOfInteraction(SP::Interaction inter) const
{
  InteractionOSIConstIterator it = _interactionOsiMap.find(inter);
  if (it == _interactionOsiMap.end())
    RuntimeException::selfThrow("Simulation::integratorOfInteraction(inter), inter not found in the integrator set.");
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

void Simulation::addInteractionInOSIMap(SP::Interaction inter, SP::OneStepIntegrator  osi)
{
  if (_interactionOsiMap.find(inter) != _interactionOsiMap.end())
    // in the map with another
    // integrator
    ;/*RuntimeException::selfThrow("Simulation::addInOSIMap(ds,osi), ds is already associated with another one-step integrator");  */
  _interactionOsiMap[inter] = osi;
}

SP::OneStepNSProblem Simulation::oneStepNSProblem(int Id)
{
  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - oneStepNSProblem(Id) - The One Step NS Problem is not in the simulation.");

  return (*_allNSProblems)[Id];
}

void Simulation::updateIndexSets()
{
  // Warning, I0 is not updated and must remain unchanged !
  unsigned int nindexsets = model()->nonSmoothDynamicalSystem()
                            ->topology()->indexSetsSize();
  if (nindexsets > 1)
  {
    for (unsigned int i = 1; i < nindexsets ; ++i)
    {
      updateIndexSet(i);
      model()->nonSmoothDynamicalSystem()->topology()->indexSet(i)->update_vertices_indices();
      model()->nonSmoothDynamicalSystem()->topology()->indexSet(i)->update_edges_indices();

    }
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
  if (!_allNSProblems->empty())
  {
    double time = model()->currentTime(); // init with current model time

    ComputeLevelsForInputAndOutput();

    std::for_each(allInteractions->begin(), allInteractions->end(),
                  boost::bind(&Interaction::initialize, _1, time));

    initOSNS();
  }

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
  //===
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
        (*itds)->initialize(model()->t0(),
                            (*itosi)->getSizeMem());
        addInOSIMap(*itds, *itosi);
      }

      (*itosi)->setSimulationPtr(shared_from_this());
      (*itosi)->initialize();

    }
  }



  // === IndexSets building ===
  SP::InteractionsSet allInteractions =
    model()->nonSmoothDynamicalSystem()->interactions();
  //  if( !allInteractions->isEmpty() ) // ie if some Interactions
  //  have been declared
  {
    ComputeLevelsForInputAndOutput();
    if (_allNSProblems->empty())
    {
      topo->indexSetsResize(0);
    }
    else
    {
      if (_levelsAreComputed)
      {
        topo->indexSetsResize(_levelMaxForOutput + 1);
      }
      else
      {
        topo->indexSetsResize(LEVELMAX);
        // ComputeLevelsForInputAndOutput will resize the indexSets when some interactions appear
      }

      std::for_each(allInteractions->begin(), allInteractions->end(),
                    boost::bind(&Interaction::initialize, _1, _tinit));

      for (unsigned int i = 1; i < topo->indexSetsSize(); ++i)
        topo->resetIndexSetPtr(i);
      // Initialize OneStepNSProblem: in derived classes specific functions.
      initOSNS();
    }
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
  //std::cout << "     Simulation::reset()"  <<std::endl;
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
  {
    (*itOsns)->saveInMemory();
    (*itOsns)->saveTimeStepInMemory();
  }
}

int Simulation::computeOneStepNSProblem(int Id)
{

  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem == NULL, Id: " + Id);

  return (*_allNSProblems)[Id]->compute(model()->currentTime());
}

void Simulation::update()
{
  assert(0);
  // for(unsigned int i = 1; i<_levelMax; ++i)
  //   update(i);
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

void Simulation::updateInput(unsigned int level)
{
  // To compute input(level) (ie with lambda[level]) for all Interactions.
  assert(level >= 0);

  //  double time = nextTime();
  double time = model()->currentTime();
  SP::Topology topology = model()->nonSmoothDynamicalSystem()->topology();
  InteractionsIterator it;

  // Set dynamical systems non-smooth part to zero.
  reset();

  // We compute input using lambda(level).
  for (it = topology->interactions()->begin();
       it != topology->interactions()->end(); it++)
  {
    assert((*it)->lowerLevelForInput() <= level);
    assert((*it)->upperLevelForInput() >= level);
    (*it)->computeInput(time, level);
  }
}

void Simulation::updateOutput(unsigned int level)
{
  // To compute output(level) (ie with y[level]) for all Interactions.
  assert(level >= 0);

  double time = model()->currentTime();
  SP::Topology topology = model()->nonSmoothDynamicalSystem()->topology();
  InteractionsIterator it;

  for (it = topology->interactions()->begin();
       it != topology->interactions()->end(); it++)
  {
    assert((*it)->lowerLevelForOutput() <= level);
    assert((*it)->upperLevelForOutput() >= level);
    (*it)->computeOutput(time , level);
  }
}

void Simulation::run()
{
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



/** a visitor to set the level(s) parameter in Interaction
 *  and to compute the levelMin and levelMax
 */
// class NonSmoothLaw;
// class DynamicalSystem;
// class OneStepIntegrator;
// class Moreau;

struct Simulation::SetupLevels : public SiconosVisitor
{
  SP::Simulation _parent;
  SP::Interaction _interaction;
  SP::NonSmoothLaw _nonSmoothLaw;
  SP::DynamicalSystem _ds;
  SetupLevels(SP::Simulation s, SP::Interaction inter,
              SP::DynamicalSystem ds) :
    _parent(s), _interaction(inter), _ds(ds)
  {
    _nonSmoothLaw = inter->nonSmoothLaw();
  };

  void visit(const Moreau&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {

      if (Type::name(*_parent) == "TimeStepping")
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 1;
        upperLevelForInput = 1;
      }
      else if (Type::name(*_parent) == "TimeSteppingProjectOnConstraints")
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 0;
        upperLevelForInput = 1;
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit - unknown simulation type: " + Type::name(*_parent));
    }
    else if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      if (Type::name(*_parent) == "TimeStepping")
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 0;
        lowerLevelForInput = 0;
        upperLevelForInput = 0;
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit - unknown simulation type: " + Type::name(*_parent));
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);
    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };
  void visit(const SchatzmanPaoli&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {

      if (Type::name(*_parent) == "TimeStepping")
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 0;
        lowerLevelForInput = 0;
        upperLevelForInput = 0;
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit - unknown simulation type: " + Type::name(*_parent));
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);

    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(2);
  };
  void visit(const D1MinusLinear&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    /** there is only a test on the dstype and simulation since  we assume that
     * we implicitely the nonsmooth law when a DS type is chosen
     */

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {


      if (Type::name(*_parent) == "TimeSteppingD1Minus")
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 2 ;
        lowerLevelForInput = 1;
        upperLevelForInput = 2;
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit(const D1MinusLinear&) - unknown simulation type: " + Type::name(*_parent));
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit(const D1MinusLinear&) - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);

    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);
    _interaction->setSteps(1);
  };



  void visit(const Lsodar&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    /** there is only a test on the dstype and simulation since  we assume that
     * we implicitely the nonsmooth law when a DS type is chosen
     */

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {
      if (Type::name(*_parent) == "EventDriven")
      {
        Type::Siconos nslType = Type::value(*_nonSmoothLaw);

        if (nslType == Type::NewtonImpactNSL || nslType == Type::MultipleImpactNSL)
        {
          lowerLevelForOutput = 0;
          upperLevelForOutput = 2 ;
          lowerLevelForInput = 1;
          upperLevelForInput = 2;
        }
        else if (nslType ==  Type::NewtonImpactFrictionNSL)
        {
          lowerLevelForOutput = 0;
          upperLevelForOutput = 4;
          lowerLevelForInput = 1;
          upperLevelForInput = 2;
          RuntimeException::selfThrow("Simulation::SetupLevels::visit - simulation of type: " + Type::name(*_parent) + " not yet implemented for nonsmooth law of type NewtonImpactFrictionNSL");
        }
        else
        {
          RuntimeException::selfThrow("Simulation::SetupLevels::visit - simulation of type: " + Type::name(*_parent) + "not yet implemented  for nonsmooth of type");
        }
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit - unknown simulation type: " + Type::name(*_parent));
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);



    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };


};

void Simulation::ComputeLevelsForInputAndOutput(SP::Interaction inter)
{
  /** \warning. We test only for the first Dynamical of the interaction.
   * we assume that the osi(s) are consistent for one interaction
   */
  SP::DynamicalSystem ds = *(inter->dynamicalSystemsBegin());
  SP::OneStepIntegrator Osi =  integratorOfDS(ds);
  addInteractionInOSIMap(inter, Osi);

  boost::shared_ptr<SetupLevels> setupLevels;
  setupLevels.reset(new SetupLevels(shared_from_this(), inter, ds));
  Osi->accept(*(setupLevels.get()));
}

void Simulation::ComputeLevelsForInputAndOutput()
{
  SP::InteractionsSet allInteractions =
    model()->nonSmoothDynamicalSystem()->interactions();
  if (_allNSProblems->empty())
  {
    _levelsAreComputed = true;
    model()->nonSmoothDynamicalSystem()->topology()->indexSetsResize(_levelMaxForOutput);
  }
  else if (not _staticLevels or not _levelsAreComputed)
  {
    _levelMinForInput = LEVELMAX;
    _levelMaxForInput = 0;
    _levelMinForOutput = LEVELMAX;
    _levelMaxForOutput = 0;

    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      ComputeLevelsForInputAndOutput(*it);
      _levelsAreComputed = true;
    }
    if (_levelsAreComputed)
    {
      if (model()->nonSmoothDynamicalSystem()->topology()->indexSetsSize() == LEVELMAX
          or  model()->nonSmoothDynamicalSystem()->topology()->indexSetsSize() < _levelMaxForOutput + 1)
      {
        model()->nonSmoothDynamicalSystem()->topology()->indexSetsResize(_levelMaxForOutput + 1);
      }

    }
  }
#include <debug.h>
  // #define DEBUG_MESSAGES 1
  DEBUG_PRINTF("_levelMinForInput =%d\n", _levelMinForInput);
  DEBUG_PRINTF("_levelMaxForInput =%d\n", _levelMaxForInput);
  DEBUG_PRINTF("_levelMinForOutput =%d\n", _levelMinForInput);
  DEBUG_PRINTF("_levelMaxForOutput =%d\n", _levelMaxForInput);
}

