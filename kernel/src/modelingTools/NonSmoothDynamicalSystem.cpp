/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "NonSmoothDynamicalSystem.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"

#include <SiconosConfig.h>
#include <functional>
using namespace std::placeholders;

#include <limits>

// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"


using namespace RELATION;

//  constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(double t0, double T):
  _t0(t0), _T(T)
{
  // === Builds an empty topology ===
  _topology.reset(new Topology());
  // we push a first element in the list to avoid acces to null when
  // we call --_changeLog.end();
  _changeLog.push_back(Change(clearTopology));
  DEBUG_EXPR((--_changeLog.end())->display());

  // see Simulation::initialize() for an explanation of why we
  // implement this changelog
};

NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  clear();
}

// changelog
void NonSmoothDynamicalSystem::Change::display() const
{
  std::cout << "Changes display   " << this <<std::endl;
  if(typeOfChange == addDynamicalSystem)
  {
    std::cout << "typeOfChange : " << typeOfChange << " : addDynamicalSystem" << std::endl;
  }
  else if(typeOfChange == rmDynamicalSystem)
  {
    std::cout << "typeOfChange : " << typeOfChange << " : rmDynamicalSystem" << std::endl;
  }
  else if(typeOfChange == addInteraction)
  {
    std::cout << "typeOfChange : " << typeOfChange << " : addInteraction" << std::endl;
  }
  else if(typeOfChange == rmInteraction)
  {
    std::cout << "typeOfChange : " << typeOfChange << " : rmInteraction" << std::endl;
  }
  else if(typeOfChange == clearTopology)
  {
    std::cout << "typeOfChange : " << typeOfChange << " : clearTopology" << std::endl;
  }
}

void NonSmoothDynamicalSystem::clearChangeLogTo(const ChangeLogIter& it)
{
  /* Given an interator into the changelog list, clear everything that
   * comes before it. User must be careful calling this if he has two
   * simulations, but in the one-simulation case (currently 100% of
   * cases), calling this will prevent changelog from building up
   * forever. Important especially for simulations using an
   * InteractionManager, e.g. mechanics_run.py. */
  while(_changeLog.begin() != it.it)
  {
    _changeLog.pop_front();
    assert((_changeLog.end() != it.it) && (_changeLog.begin() != _changeLog.end())
           && "NSDS::clearChangeLogTo: iterator not in list!");
  }
}


// === DynamicalSystems management ===

void NonSmoothDynamicalSystem::display() const
{
  std::cout << " ===== Non Smooth Dynamical System display ===== " <<std::endl;
  std::cout << "---> isBVP = " << _BVP <<std::endl;
  dynamicalSystems()->begin();
  _topology->indexSet0()->display();
  std::cout << "---> last change : " <<std::endl;
  (--_changeLog.end())->display();
  std::cout << "===================================================" <<std::endl;
}

void  NonSmoothDynamicalSystem::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  // some checks here ...
  if(!ds)
  {
    THROW_EXCEPTION("NonSmoothDynamicalSystem::insertDynamicalSystem :: DS is nul");
  }

  // Do not insert the same ds several times : results in errors in initialisation process.
  if(! _topology->hasDynamicalSystem(ds))
  {
    _topology->insertDynamicalSystem(ds);
    _changeLog.push_back(Change(addDynamicalSystem,ds));
    _mIsLinear = ((ds)->isLinear() && _mIsLinear);
  }
}

void  NonSmoothDynamicalSystem::removeDynamicalSystem(SP::DynamicalSystem ds)
{
  _topology->removeDynamicalSystem(ds);
  _changeLog.push_back(Change(rmDynamicalSystem,ds));
}
void  NonSmoothDynamicalSystem::removeInteraction(SP::Interaction inter)
{
  _topology->removeInteraction(inter);
  _changeLog.push_back(Change(rmInteraction,inter));
}

void NonSmoothDynamicalSystem::link(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2)
{
  _mIsLinear = (inter->relation()->isLinear() && _mIsLinear);
  _topology->link(inter, ds1, ds2);
  _changeLog.push_back(Change(addInteraction,inter));
};


void NonSmoothDynamicalSystem::clear()
{
  _topology->clear();
  _changeLog.push_back(Change(clearTopology));
}

void NonSmoothDynamicalSystem::setSymmetric(bool val)
{
  _topology->setSymmetric(val);
}


void NonSmoothDynamicalSystem::reset()
{
  DynamicalSystemsGraph::VIterator vi;
  for(vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dynamicalSystems()->bundle(*vi)->resetNonSmoothPart(1);
  }
}

void NonSmoothDynamicalSystem::reset(unsigned int level)
{
  DynamicalSystemsGraph::VIterator vi;
  for(vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dynamicalSystems()->bundle(*vi)->resetNonSmoothPart(level);
  }
}

void NonSmoothDynamicalSystem::swapInMemory()
{
  //could be better to call bind method
  DynamicalSystemsGraph::VIterator vi;
  for(vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dynamicalSystems()->bundle(*vi)->swapInMemory();
  }
}
void NonSmoothDynamicalSystem::pushInteractionsInMemory()
{
  // Save Interactions state into Memory.

  if(_topology->indexSet0()->size() > 0)
  {
    // Temp FP : saveInOldVar was called for each osns and each osns call
    // swapInOldVar for all interactions in the nsds.
    // ==> let's do it only once, by the simu.

    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet0 = _topology->indexSet0();
    for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      indexSet0->bundle(*ui)->swapInMemory();
    }
  }
}
void NonSmoothDynamicalSystem::updateDSPlugins(double time)
{
  //could be better to call bind method
  DynamicalSystemsGraph::VIterator vi;
  for(vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dynamicalSystems()->bundle(*vi)->updatePlugins(time);
  }
}
void NonSmoothDynamicalSystem::updateInput(double time, unsigned int level)
{

  DEBUG_BEGIN("Nonsmoothdynamicalsystem::updateInput(double time, unsigned int level)\n");
  DEBUG_PRINTF("with level = %i\n", level);


  // To compute input(level) (ie with lambda[level]) for all Interactions.
  //  assert(level>=0);

  // Set dynamical systems non-smooth part to zero.
  reset(level);

  // We compute input using lambda(level).
  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _topology->indexSet0();
  for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForInput() <= level);
    assert(inter->upperLevelForInput() >= level);
    inter->computeInput(time, level);
  }

  DEBUG_END("Nonsmoothdynamicalsystem::updateInput(double time, unsigned int level)\n");

}


void NonSmoothDynamicalSystem::updateOutput(double time, unsigned int level)
{

  // To compute output(level) (ie with y[level]) for all Interactions.
  //  assert(level>=0);

  DEBUG_BEGIN("NonSmoothDynamicalSystem::updateOutput(unsigned int level)\n");
  DEBUG_PRINTF("with level = %i\n", level);
  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _topology->indexSet0();
  for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    inter->computeOutput(time, level);
  }
  DEBUG_END("NonSmoothDynamicalSystem::updateOutput(unsigned int level)\n");
}


void NonSmoothDynamicalSystem::updateOutput(double time, unsigned int level_min, unsigned int level_max)
{

  // To compute output(level) (ie with y[level]) for all Interactions in I0
  // and for a range of levels in a single pass through I0.
  //  assert(level>=0);

  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _topology->indexSet0();
  for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level_max);
    assert(inter->upperLevelForOutput() >= level_min);
    for(unsigned int level = level_min; level<=level_max; ++level)
      inter->computeOutput(time, level);
  }
}

void NonSmoothDynamicalSystem::computeInteractionJacobians(double time)
{

  DEBUG_BEGIN("NonSmoothDynamicalSystem::computeInteractionJacobians(double time)\n");
  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _topology->indexSet0();
  for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    inter->relation()->computeJach(time, *inter);
    inter->relation()->computeJacg(time, *inter);
  }
  DEBUG_END("NonSmoothDynamicalSystem::computeInteractionJacobians(double time)\n");
}

void NonSmoothDynamicalSystem::computeInteractionJacobians(double time, InteractionsGraph& indexSet)
{
  DEBUG_BEGIN("NonSmoothDynamicalSystem::computeInteractionJacobians(double time)\n");
  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
  {
    inter = indexSet.bundle(*ui);
    inter->relation()->computeJach(time, *inter);
    inter->relation()->computeJacg(time, *inter);
  }
  DEBUG_END("NonSmoothDynamicalSystem::computeInteractionJacobians(double time)\n");
}

void NonSmoothDynamicalSystem::visitDynamicalSystems(SP::SiconosVisitor visitor)
{
  DynamicalSystemsGraph &dsg = *dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg.vertices();
  for(; dsi != dsiend; ++dsi)
  {
    dsg.bundle(*dsi)->acceptSP(visitor);
  }
}
std::vector<SP::DynamicalSystem> NonSmoothDynamicalSystem::dynamicalSystemsVector() const
{
  std::vector<SP::DynamicalSystem> dynamicalSystemsVector;
  DynamicalSystemsGraph &dsg = *dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg.vertices();
  for(; dsi != dsiend; ++dsi)
  {
    dynamicalSystemsVector.push_back(dsg.bundle(*dsi));
  }

  return dynamicalSystemsVector;
}
std::vector<SP::Interaction> NonSmoothDynamicalSystem::InteractionsVector() const
{
  std::vector<SP::Interaction> interactionsVector;
  SP::InteractionsGraph indexSet0 = _topology->indexSet0();
  InteractionsGraph::VIterator ui, uiend;
  for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    interactionsVector.push_back(indexSet0->bundle(*ui));
  }

  return interactionsVector;
}
