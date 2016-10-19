/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "Simulation.hpp"
#include "EventDriven.hpp"
#include "DynamicalSystem.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"
#include "EventsManager.hpp"
#include "LagrangianDS.hpp"

// One Step Integrators
#include "EulerMoreauOSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "LsodarOSI.hpp"
#include "Hem5OSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "D1MinusLinearOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"
#include "ZeroOrderHoldOSI.hpp"
// One Step Non Smooth Problems
#include "LCP.hpp"
#include "QP.hpp"
#include "Relay.hpp"
#include "NonSmoothLaw.hpp"
#include "TypeName.hpp"
// for Debug
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <debug.h>
#include <fstream>

#include "Model.hpp"



// --- Constructor with a TimeDiscretisation (and thus a Model) and an
// --- id ---
Simulation::Simulation(SP::TimeDiscretisation td):
  _name("unnamed"), _tinit(0.0), _tend(0.0), _tout(0.0),
  _numberOfIndexSets(0),
  _tolerance(DEFAULT_TOLERANCE), _printStat(false),
  _staticLevels(false)
{
  if (!td)
    RuntimeException::selfThrow("Simulation constructor - timeDiscretisation == NULL.");
  _useRelativeConvergenceCriterion = false;
  _relativeConvergenceCriterionHeld = false;
  _relativeConvergenceTol = 10e-3;

  // === indexSets will be updated during initialize() call ===

  _allOSI.reset(new OSISet());
  _allNSProblems.reset(new OneStepNSProblems());
  _eventsManager.reset(new EventsManager(td)); //
}

// --- Destructor ---
Simulation::~Simulation()
{
  clear();
  // -> see shared ressources for this
  if (statOut.is_open()) statOut.close();
}

void Simulation::setTimeDiscretisationPtr(SP::TimeDiscretisation td)
{
  _eventsManager->setTimeDiscretisationPtr(td);
}

double Simulation::getTk() const
{
  return _eventsManager->getTk();
}

double Simulation::getTkp1() const
{
  return _eventsManager->getTkp1();
}

double Simulation::getTkp2() const
{
  return _eventsManager->getTkp2();
}

double Simulation::currentTimeStep() const
{
  return _eventsManager->currentTimeStep();
}

double Simulation::startingTime() const
{
  return _eventsManager->startingTime();
}

double Simulation::nextTime() const
{
  return _eventsManager->nextTime();
}

bool Simulation::hasNextEvent() const
{
  return _eventsManager->hasNextEvent();
}


// clear all maps to break shared_ptr cycle
void Simulation::clear()
{
  if (_allOSI)
  {
    _allOSI->clear();
  }
  if (_allNSProblems)
  {
    _allNSProblems->clear();
  }
}

// Getters/setters

void Simulation::insertIntegrator(SP::OneStepIntegrator osi)
{
  _allOSI->insert(osi);
}

SP::InteractionsGraph Simulation::indexSet(unsigned int i)
{
  return _nsds->topology()->indexSet(i) ;
}

SP::OneStepNSProblem Simulation::oneStepNSProblem(int Id)
{
  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - oneStepNSProblem(Id) - The One Step NS Problem is not in the simulation.");

  return (*_allNSProblems)[Id];
}

void Simulation::updateIndexSets()
{

  DEBUG_BEGIN("Simulation::updateIndexSets()\n");
  // update I0 indices
  unsigned int nindexsets = _nsds->topology()->indexSetsSize();

  DEBUG_PRINTF("  nindexsets = %d\n", nindexsets   );
  if (nindexsets > 1)
  {
    for (unsigned int i = 1; i < nindexsets ; ++i)
    {
      updateIndexSet(i);
      _nsds->topology()->indexSet(i)->update_vertices_indices();
      _nsds->topology()->indexSet(i)->update_edges_indices();
    }
  }
  DEBUG_END("Simulation::updateIndexSets()\n");

}

void Simulation::insertNonSmoothProblem(SP::OneStepNSProblem osns, int Id)
{
  if (_allNSProblems->size() > 0 && ((*_allNSProblems)[Id]))
    RuntimeException::selfThrow("Simulation - insertNonSmoothProblem(osns), trying to insert a OSNSP already existing. ");
  (*_allNSProblems)[Id] = osns;

}

void Simulation::initialize(SP::Model m, bool withOSI)
{
  // === Connection with the model ===
  assert(m && "Simulation::initialize(model) - model = NULL.");

  _T = m->finalT();

  _nsds =  m->nonSmoothDynamicalSystem();

  // === Events manager initialization ===
  _eventsManager->initialize(_T);
  _tinit = _eventsManager->startingTime();
  //===


  if (withOSI)
  {


    if (numberOfOSI() == 0)
      RuntimeException::selfThrow("Simulation::initialize No OSI !");


    DynamicalSystemsGraph::VIterator dsi, dsend;
    SP::DynamicalSystemsGraph DSG = _nsds->topology()->dSG(0);
    for (std11::tie(dsi, dsend) = DSG->vertices(); dsi != dsend; ++dsi)
    {
      SP::OneStepIntegrator osi = DSG->properties(*dsi).osi;
      SP::DynamicalSystem ds = DSG->bundle(*dsi);
      if (!osi)
      {
        // By default, if the user has not set the OSI, we assign the first OSI to all DS
        _nsds->topology()->setOSI(ds,*_allOSI->begin());
        //std::cout << "By default, if the user has not set the OSI, we assign the first OSI to all DS"<<std::endl;
      }

      osi = DSG->properties(*dsi).osi;
      ds->initialize(m->t0(), osi->getSizeMem());
    }


    // === OneStepIntegrators initialization ===
    for (OSIIterator itosi = _allOSI->begin();
         itosi != _allOSI->end(); ++itosi)
    {
      // for (DSIterator itds = (*itosi)->dynamicalSystems()->begin();
      //      itds != (*itosi)->dynamicalSystems()->end();
      //      ++itds)
      // {
      //   (*itds)->initialize(startingTime(),
      //                       (*itosi)->getSizeMem());
      //   addInOSIMap(*itds, *itosi);
      // }

      (*itosi)->setSimulationPtr(shared_from_this());
      (*itosi)->initialize(*m);

    }
  }

  // This is the default
  _levelMinForInput = LEVELMAX;
  _levelMaxForInput = 0;
  _levelMinForOutput = LEVELMAX;
  _levelMaxForOutput = 0;

  computeLevelsForInputAndOutput();

  // Loop over all DS in the graph, to reset NS part of each DS.
  // Note FP : this was formerly done in inter->initialize call with local levels values
  // but I think it's ok (better?) to do it with the simulation levels values.
  DynamicalSystemsGraph::VIterator dsi, dsend;
  SP::DynamicalSystemsGraph DSG = _nsds->topology()->dSG(0);
  for (std11::tie(dsi, dsend) = DSG->vertices(); dsi != dsend; ++dsi)
  {
    //assert(_levelMinForInput <= _levelMaxForInput);
    for (unsigned int k = _levelMinForInput ; k < _levelMaxForInput + 1; k++)
    {
      DSG->bundle(*dsi)->initializeNonSmoothInput(k);
    }
  }

  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet0->bundle(*ui);
    inter.initialize(_tinit, indexSet0->properties(*ui));
  }

  // Initialize OneStepNSProblem(s). Depends on the type of simulation.
  // Warning FP : must be done in any case, even if the interactions set
  // is empty.
  initOSNS();

  // Process events at time _tinit. Useful to save values in memories
  // for example.  Warning: can not be called during
  // eventsManager->initialize, because it needs the initialization of
  // OSI, OSNS ...
  _eventsManager->preUpdate(*this);

  _tend =  _eventsManager->nextTime();

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
    statOut.open("simulationStat.dat", std::ios::out | std::ios::trunc);
    if (!statOut.is_open())
      SiconosVectorException::selfThrow("writing error : Fail to open file simulationStat.dat ");
    statOut << "============================================" <<std::endl;
    statOut << " Siconos Simulation of type " << Type::name(*this) << "." <<std::endl;
    statOut <<std::endl;
    statOut << "The tolerance parameter is equal to: " << _tolerance <<std::endl;
    statOut <<std::endl <<std::endl;
  }
}

void Simulation::initializeInteraction(double time, SP::Interaction inter)
{
  // determine which (lower and upper) levels are required for this Interaction
  // in this Simulation.
  computeLevelsForInputAndOutput(inter);

  // Get the interaction properties from the topology for initialization.
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();
  InteractionsGraph::VDescriptor ui = indexSet0->descriptor(inter);

  // This calls computeOutput() and initializes qMemory and q_k.
  inter->initialize(time, indexSet0->properties(ui));
}



int Simulation::computeOneStepNSProblem(int Id)
{
  DEBUG_BEGIN("Simulation::computeOneStepNSProblem(int Id)\n");
  DEBUG_PRINTF("with Id = %i\n", Id);

  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem == NULL, Id: " + Id);

  DEBUG_END("Simulation::computeOneStepNSProblem(int Id)\n");
  return (*_allNSProblems)[Id]->compute(nextTime());


}


SP::SiconosVector Simulation::y(unsigned int level, unsigned int coor)
{
  // return output(level) (ie with y[level]) for all Interactions.
  // assert(level>=0);

  DEBUG_BEGIN("Simulation::output(unsigned int level, unsigned int coor)\n");
  DEBUG_PRINTF("with level = %i and coor = %i \n", level,coor);

  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();

  SP::SiconosVector y (new SiconosVector (_nsds->topology()->indexSet0()->size() ));
  int i=0;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    y->setValue(i,inter->y(level)->getValue(coor));
    i++;
  }
  DEBUG_END("Simulation::output(unsigned int level, unsigned int coor)\n");
  return y;
}

SP::SiconosVector Simulation::lambda(unsigned int level, unsigned int coor)
{
  // return input(level) (ie with lambda[level]) for all Interactions.
  // assert(level>=0);

  DEBUG_BEGIN("Simulation::input(unsigned int level, unsigned int coor)\n");
  DEBUG_PRINTF("with level = %i and coor = %i \n", level,coor);

  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();

  SP::SiconosVector lambda (new SiconosVector (_nsds->topology()->indexSet0()->size() ));
  int i=0;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    lambda->setValue(i,inter->lambda(level)->getValue(coor));
    i++;
  }
  DEBUG_END("Simulation::input(unsigned int level, unsigned int coor)\n");
  return lambda;
}


void Simulation::run()
{
  unsigned int count = 0; // events counter.

  std::cout << " ==== Start of " << Type::name(*this) << " simulation - This may take a while ... ====" <<std::endl;
  while (hasNextEvent())
  {
    advanceToEvent();
    processEvents();
    count++;
  }
  std::cout << "===== End of " << Type::name(*this) << " simulation. " << count << " events have been processed. ==== " <<std::endl;
}

void Simulation::processEvents()
{
  _eventsManager->processEvents(*this);

  if (_eventsManager->hasNextEvent())
  {
    // For TimeStepping Scheme, need to update IndexSets, but not for EventDriven scheme
    if (Type::value(*this) != Type::EventDriven)
    {
      updateIndexSets();
    }
  }

  /* should be evaluated only if needed */
  SP::DynamicalSystemsGraph dsGraph = _nsds->dynamicalSystems();
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
// class MoreauJeanOSI;
struct Simulation::SetupLevels : public SiconosVisitor
{
  using SiconosVisitor::visit;

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

  void visit(const MoreauJeanOSI&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {

      if (Type::value(*_parent) == Type::TimeStepping)
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 1;
        upperLevelForInput = 1;
      }
      else if (Type::value(*_parent) == Type::TimeSteppingDirectProjection)
      {
        // Warning : we never enter this case !!!
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 0;
        upperLevelForInput = 1;
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit - unknown simulation type: " + Type::name(*_parent));
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);

    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);

    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };

  void visit(const MoreauJeanGOSI&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {

      if (Type::value(*_parent) == Type::TimeStepping)
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 1;
        upperLevelForInput = 1;
      }
      else if (Type::value(*_parent) == Type::TimeSteppingDirectProjection)
      {
        // Warning : we never enter this case !!!
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 0;
        upperLevelForInput = 1;
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit - unknown simulation type: " + Type::name(*_parent));
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);

    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);

    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };

  void visit(const EulerMoreauOSI&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {


      if (Type::value(*_parent) == Type::TimeStepping)
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
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);
    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };




  void visit(const MoreauJeanDirectProjectionOSI& moreauCPOSI)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {
      if (Type::value(*_parent) == Type::TimeSteppingDirectProjection)
      {
        // Warning : we never enter this case !!!
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 0;
        upperLevelForInput = 1;
      }
      else
      {
        RuntimeException::selfThrow("Simulation::SetupLevels::visit(const MoreauJeanCombinedProjectionOSI) - unknown simulation type: " + Type::name(*_parent));
      }
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit(const MoreauJeanCombinedProjectionOSI) - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);
    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);

  };


  void visit(const MoreauJeanCombinedProjectionOSI& moreauCPOSI)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {
      if (Type::value(*_parent) == Type::TimeStepping)
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1;
        lowerLevelForInput = 1;
        upperLevelForInput = 1;
      }
      else if (Type::value(*_parent) == Type::TimeSteppingCombinedProjection)
      {
        // Warning : we never enter this case !!!
        lowerLevelForOutput = 0;
        upperLevelForOutput = 1 ;
        lowerLevelForInput = 0;
        upperLevelForInput = 1;
      }
      else
      {
        RuntimeException::selfThrow("Simulation::SetupLevels::visit(const MoreauJeanCombinedProjectionOSI) - unknown simulation type: " + Type::name(*_parent));
      }
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit(const MoreauJeanCombinedProjectionOSI) - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);


    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);

  };



  void visit(const SchatzmanPaoliOSI&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
    {

      if (Type::value(*_parent) == Type::TimeStepping)
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
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);

    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(2);
  };
  void visit(const D1MinusLinearOSI& d1OSI)
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


      if (Type::value(*_parent) == Type::TimeSteppingD1Minus)
      {
        lowerLevelForOutput = 0;
        upperLevelForOutput = 2 ;
        lowerLevelForInput = 1;
        upperLevelForInput = 2;
      }
      else
        RuntimeException::selfThrow("Simulation::SetupLevels::visit(const D1MinusLinearOSI&) - unknown simulation type: " + Type::name(*_parent));
    }
    else RuntimeException::selfThrow("Simulation::SetupLevels::visit(const D1MinusLinearOSI&) - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);

    /* Get the number of required index sets. */

    unsigned int nbIndexSets =  d1OSI.numberOfIndexSets();

    _parent->_numberOfIndexSets = std::max<int>(nbIndexSets, _parent->_numberOfIndexSets);

    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(2); // Two evaluations of lambda(2) are made for each time--step
  };



  void visit(const LsodarOSI&)
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
      if (Type::value(*_parent) == Type::EventDriven)
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
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);
    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };



  void visit(const Hem5OSI&)
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
      if (Type::value(*_parent) == Type::EventDriven)
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
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);
    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };

  void visit(const NewMarkAlphaOSI&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    /** there is only a test on the dstype and simulation since  we assume that
     * we implicitely the nonsmooth law when a DS type is chosen
     */

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      if (Type::value(*_parent) == Type::EventDriven)
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
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);
    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };




  void visit(const ZeroOrderHoldOSI&)
  {
    unsigned int lowerLevelForOutput = LEVELMAX;
    unsigned int upperLevelForOutput = 0;
    unsigned int lowerLevelForInput = LEVELMAX;
    unsigned int upperLevelForInput = 0;

    Type::Siconos dsType = Type::value(*_ds);

    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
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
    else
      RuntimeException::selfThrow("Simulation::SetupLevels::visit - not yet implemented for Dynamical system type :" + dsType);

    _parent->_levelMinForInput = std::min<int>(lowerLevelForInput, _parent->_levelMinForInput);
    _parent->_levelMaxForInput = std::max<int>(upperLevelForInput, _parent->_levelMaxForInput);
    _parent->_levelMinForOutput = std::min<int>(lowerLevelForOutput, _parent->_levelMinForInput);
    _parent->_levelMaxForOutput = std::max<int>(upperLevelForOutput, _parent->_levelMaxForInput);
    _parent->_numberOfIndexSets = std::max<int>(_parent->_levelMaxForOutput + 1, _parent->_numberOfIndexSets);
    _interaction->setLowerLevelForOutput(lowerLevelForOutput);
    _interaction->setUpperLevelForOutput(upperLevelForOutput);

    _interaction->setLowerLevelForInput(lowerLevelForInput);
    _interaction->setUpperLevelForInput(upperLevelForInput);

    _interaction->setSteps(1);
  };


};

void Simulation::computeLevelsForInputAndOutput(SP::Interaction inter, bool init)
{
  DEBUG_PRINT("Simulation::computeLevelsForInputAndOutput(SP::Interaction inter, bool init)\n");

 /** \warning. We test only for the first Dynamical of the interaction.
   * we assume that the osi(s) are consistent for one interaction
   */
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet(0);
  SP::DynamicalSystem ds = indexSet0->properties(indexSet0->descriptor(inter)).source;

  SP::DynamicalSystemsGraph DSG0 = _nsds->topology()->dSG(0);
  SP::OneStepIntegrator osi = DSG0->properties(DSG0->descriptor(ds)).osi;

  if (!osi)
    RuntimeException::selfThrow("Simulation::computeLevelsForInputAndOutput osi does not exists");
  indexSet0->properties(indexSet0->descriptor(inter)).osi = osi;
  std11::shared_ptr<SetupLevels> setupLevels;
  setupLevels.reset(new SetupLevels(shared_from_this(), inter, ds));
  osi->accept(*(setupLevels.get()));
  if (!init) // We are not computing the levels at the initialization
  {
    SP::Topology topo = _nsds->topology();
    unsigned int indxSize = topo->indexSetsSize();
    assert (_numberOfIndexSets >0);
    if ((indxSize == LEVELMAX) || (indxSize < _numberOfIndexSets ))
    {
      topo->indexSetsResize(_numberOfIndexSets);
      // Init if the size has changed
      for (unsigned int i = indxSize; i < topo->indexSetsSize(); i++) // ++i ???
        topo->resetIndexSetPtr(i);
    }
  }
}

void Simulation::computeLevelsForInputAndOutput()
{
  DEBUG_PRINT("Simulation::computeLevelsForInputAndOutput()\n");

  SP::Topology topo = _nsds->topology();

  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet0 = topo->indexSet0();
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    computeLevelsForInputAndOutput(indexSet0->bundle(*ui), true);
  }

  unsigned int indxSize = topo->indexSetsSize();
  if ((indxSize == LEVELMAX) || (indxSize < _numberOfIndexSets ))
  {
    topo->indexSetsResize(_numberOfIndexSets );
    // Init if the size has changed
    for (unsigned int i = indxSize; i < topo->indexSetsSize(); i++) // ++i ???
      topo->resetIndexSetPtr(i);
  }
  DEBUG_PRINTF("_numberOfIndexSets =%d\n", _numberOfIndexSets);
  DEBUG_PRINTF("_levelMinForInput =%d\n", _levelMinForInput);
  DEBUG_PRINTF("_levelMaxForInput =%d\n", _levelMaxForInput);
  DEBUG_PRINTF("_levelMinForOutput =%d\n", _levelMinForInput);
  DEBUG_PRINTF("_levelMaxForOutput =%d\n", _levelMaxForInput);
}

void Simulation::updateT(double T)
{
  _T = T;
  _eventsManager->updateT(T);
}


void Simulation::link(SP::Interaction inter,
                      SP::DynamicalSystem ds1,
                      SP::DynamicalSystem ds2)
{
  DEBUG_PRINTF("link interaction : %d\n", inter->number());

  nonSmoothDynamicalSystem()->link(inter, ds1, ds2);

  initializeInteraction(nextTime(), inter);

  // Note FP : ds init should probably be done once and only once for
  // all ds (like in simulation->initialize()) but where/when?
  // Note SS : in InteractionManager::buildGraph()?
  unsigned int levelMinForInput = inter->lowerLevelForInput();
  unsigned int levelMaxForInput = inter->upperLevelForInput();
  bool has2DS = inter->has2Bodies();
  for (unsigned int k = levelMinForInput ; k < levelMaxForInput + 1; k++)
  {
    ds1->initializeNonSmoothInput(k);
    if(has2DS)
      ds2->initializeNonSmoothInput(k);
  }

  _linkOrUnlink = true;
}

void Simulation::unlink(SP::Interaction inter)
{
  nonSmoothDynamicalSystem()->removeInteraction(inter);
  _linkOrUnlink = true;
}

void Simulation::updateInteractions()
{
  // Update interactions if a manager was provided
  if (_interman) {
    _linkOrUnlink = false;
    _interman->updateInteractions(shared_from_this());

    if (_linkOrUnlink)
      initOSNS();
  }
}
