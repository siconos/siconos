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
#include "NonSmoothDynamicalSystem.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"

#include <SiconosConfig.h>
#if defined(SICONOS_STD_FUNCTIONAL) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#include <functional>
using namespace std::placeholders;
#else
#include <boost/bind.hpp>
#include <boost/weak_ptr.hpp>
#endif

#include <limits>

// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "debug.h"


using namespace RELATION;

// --- CONSTRUCTORS/DESTRUCTOR ---

// Default constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(): _BVP(false), _mIsLinear(true)
{
  // === Builds an empty topology ===
  _topology.reset(new Topology());
};


NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  clear();
}

// === DynamicalSystems management ===

void NonSmoothDynamicalSystem::display() const
{
  std::cout << " ===== Non Smooth Dynamical System display ===== " <<std::endl;
  std::cout << "---> isBVP = " << _BVP <<std::endl;
  dynamicalSystems()->begin();
  _topology->indexSet0()->display();
  std::cout << "===================================================" <<std::endl;
}

void NonSmoothDynamicalSystem::link(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2)
{
  _mIsLinear = (inter->relation()->isLinear() && _mIsLinear);
  _topology->link(inter, ds1, ds2);
};


void NonSmoothDynamicalSystem::clear()
{
  _topology->clear();
}

void NonSmoothDynamicalSystem::setSymmetric(bool val)
{
  _topology->setSymmetric(val);
}


void NonSmoothDynamicalSystem::reset()
{
  DynamicalSystemsGraph::VIterator vi;
  for (vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dynamicalSystems()->bundle(*vi)->resetNonSmoothPart(1);
  }
}

void NonSmoothDynamicalSystem::reset(unsigned int level)
{
  DynamicalSystemsGraph::VIterator vi;
  for (vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dynamicalSystems()->bundle(*vi)->resetNonSmoothPart(level);
  }
}

void NonSmoothDynamicalSystem::swapInMemory()
{
  //could be better to call bind method
  DynamicalSystemsGraph::VIterator vi;
  for (vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dynamicalSystems()->bundle(*vi)->swapInMemory();
  }
}
void NonSmoothDynamicalSystem::pushInteractionsInMemory()
{
  // Save Interactions state into Memory.

  if (_topology->indexSet0()->size() > 0)
  {
    // Temp FP : saveInOldVar was called for each osns and each osns call
    // swapInOldVar for all interactions in the nsds.
    // ==> let's do it only once, by the simu.

    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet0 = _topology->indexSet0();
    for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      indexSet0->bundle(*ui)->swapInOldVariables();
      indexSet0->bundle(*ui)->swapInMemory();
    }
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
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForInput() <= level);
    assert(inter->upperLevelForInput() >= level);
    inter->computeInput(time, indexSet0->properties(*ui), level);
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
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    inter->computeOutput(time, indexSet0->properties(*ui), level);
  }
  DEBUG_END("NonSmoothDynamicalSystem::updateOutput(unsigned int level)\n");

}


void NonSmoothDynamicalSystem::visitDynamicalSystems(SP::SiconosVisitor visitor)
{
  DynamicalSystemsGraph &dsg = *dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std11::tie(dsi, dsiend) = dsg.vertices();
  for (; dsi != dsiend; ++dsi)
  {
    dsg.bundle(*dsi)->acceptSP(visitor);
  }
}
