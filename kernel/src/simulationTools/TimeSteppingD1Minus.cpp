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

#include "TimeSteppingD1Minus.hpp"
#include "D1MinusLinearOSI.hpp"
#include "TimeDiscretisation.hpp"
#include "Topology.hpp"
//#include "Interaction.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "TypeName.hpp"
#include "NonSmoothLaw.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "EventsManager.hpp"
#include "OneStepNSProblem.hpp"

#include <ciso646>

#include <SiconosConfig.h>
#if defined(SICONOS_STD_FUNCTIONAL) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#include <functional>
using namespace std::placeholders;
#else
#include <boost/bind.hpp>
#include <boost/weak_ptr.hpp>
#endif


using namespace RELATION;

void TimeSteppingD1Minus::initOSNS()
{
  // initialize OSNS for InteractionsGraph from Topology
  assert(_nsds->topology()->isUpToDate());
  SP::Topology topo =  _nsds->topology();

  // there is at least one OSNP
  if (!_allNSProblems->empty())
  {
    if (_allNSProblems->size() != 2)
      RuntimeException::selfThrow("TimeSteppingD1Minus::initOSNS, TimeSteppingD1Minus simulation must have two OneStepNonsmoothProblems.");

    //update all index sets
    updateIndexSets();

    // update output
    updateOutput();
  }
}

TimeSteppingD1Minus::TimeSteppingD1Minus(SP::TimeDiscretisation td, int nb) : Simulation(td)
{
  (*_allNSProblems).resize(nb);
}

TimeSteppingD1Minus::~TimeSteppingD1Minus()
{
}

void TimeSteppingD1Minus::updateIndexSet(unsigned int i)
{
  DEBUG_PRINTF("\nTimeSteppingD1Minus::updateIndexSet(unsigned int i) for i = %i\n", i);
  // To update IndexSet i: add or remove Interactions from
  // this set, depending on y values.

  assert(_nsds);
  assert(_nsds->topology());

  SP::Topology topo = _nsds->topology();

  assert(i < topo->indexSetsSize() &&
         "TimeSteppingD1Minus::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to Topology.
  assert(i > 0 &&
         "TimeSteppingD1Minus::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  SP::InteractionsGraph indexSet0 = topo->indexSet(0); // ALL Interactions : formula (8.30) of Acary2008
  SP::InteractionsGraph indexSetCurrent = topo->indexSet(i); // ACTIVE Interactions for IMPACTS
  assert(indexSet0);
  assert(indexSetCurrent);
  DynamicalSystemsGraph& DSG0= *nonSmoothDynamicalSystem()->dynamicalSystems();
  topo->setHasChanged(false); // only with changed topology, OSNS will be forced to update themselves


  DEBUG_PRINTF("\nINDEXSETS BEFORE UPDATE for level i = %i\n", i);
  DEBUG_PRINTF(" indexSet0 size : %ld\n", indexSet0->size());
  DEBUG_PRINTF(" indexSet(%i) size : %ld\n", i, topo->indexSet(i)->size());

  InteractionsGraph::VIterator uipend, uip;
  for (std11::tie(uip, uipend) = indexSet0->vertices(); uip != uipend; ++uip)
    /* loop over ALL vertices in indexSet0 */
  {

    SP::Interaction inter = indexSet0->bundle(*uip);

    // We assume that the integrator of the ds1 drive the update of the index set
    //SP::OneStepIntegrator Osi = indexSetCurrent->properties(*uip).osi;
    SP::DynamicalSystem ds1 = indexSetCurrent->properties(*uip).source;
    OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;
    unsigned int levelMaxForInput = osi.levelMaxForInput();
    if ((!indexSetCurrent->is_vertex(inter))
        and (osi.addInteractionInIndexSet(inter, i)))
    {
      indexSetCurrent->copy_vertex(inter, *indexSet0);
      topo->setHasChanged(true);
    }
    else if ((indexSetCurrent->is_vertex(inter))
             and !(osi.addInteractionInIndexSet(inter, i)))
    {
      indexSetCurrent->remove_vertex(inter);
      topo->setHasChanged(true);
      if (i <= levelMaxForInput)
      {
        DEBUG_PRINTF("Reset to zero inter->lambda(%i)", i);
        inter->lambda(i)->zero();
      }
    }

    if (!indexSetCurrent->is_vertex(inter))
    {
      DEBUG_PRINTF("The current interaction is not in the indexSet(%i)\n",(int)i);
      if (i <= levelMaxForInput)
      {
        DEBUG_EXPR(inter->lambda(i)->display());
        inter->lambda(i)->zero();
      }
    }
    else
    {
      DEBUG_PRINTF("The current interaction is in the indexSet(%i)\n",(int)i);
      DEBUG_EXPR(if (i <= levelMaxForInput) inter->lambda(i)->display());
    }




  }
  DEBUG_PRINTF("\nINDEXSETS AFTER UPDATE for level i = %i\n", i);
  DEBUG_PRINTF(" indexSet0 size : %ld\n", indexSet0->size());
  DEBUG_PRINTF(" indexSet(%i) size : %ld\n", i, topo->indexSet(i)->size());
}

void TimeSteppingD1Minus::run()
{
  unsigned int count = 0;
  std::cout << " ==== Start of " << Type::name(*this) << " simulation - This may take a while ... ====" <<std::endl;
  while (_eventsManager->hasNextEvent())
  {
    advanceToEvent();

    processEvents();
    count++;
  }
  std::cout << "===== End of " << Type::name(*this) << "simulation. " << count << " events have been processed. ==== " <<std::endl;
}

void TimeSteppingD1Minus::advanceToEvent()
{
  // Update interactions if a manager was provided
  updateInteractions();

  // we start after initialization (initOSNS) with
  // * initial state (q_0, v_0^+)
  // * updated indexset (I_0^+)
  // * updated  gaps and gap velocities (g_0^+)
  //
  // hence we end this procedure with
  // * state (q_{k+1}, v_{k+1}^+)
  // * updated gaps and gap velocities (g_{k+1}^+)
  // * indexset (I_{k+1}^+)

  // Initialize lambdas of all interactions.
  SP::InteractionsGraph indexSet0 = _nsds->
                                    topology()->indexSet(0);
  InteractionsGraph::VIterator ui, uiend, vnext;
  std11::tie(ui, uiend) = indexSet0->vertices();
  for (vnext = ui; ui != uiend; ui = vnext)
  {
    ++vnext;
    indexSet0->bundle(*ui)->resetAllLambda();
  }

  // calculate residu without nonsmooth event with OSI
  // * calculate position q_{k+1} in ds->q()
  // * calculate velocity v_{k+1}^- and not free velocity in ds->velocity()
  // * calculate free residu in ds->free()
  computeResidu();

  // calculate state without nonsmooth event with OSI
  // * calculate free velocity and not v_{k+1}^- in ds->velocity
  computeFreeState();

  // event (impulse) calculation only when there has been a topology change (here: closing contact)
  // * calculate gap velocity using free velocity with OSI
  // * calculate local impulse (Lambda_{k+1}^+)
  //
  // Maurice Bremond: indices must be recomputed
  // as we deal with dynamic graphs, vertices and edges are stored
  // in lists for fast add/remove during updateIndexSet(i)
  // we need indices of list elements to build the OSNS Matrix so we
  // need an update if graph has changed
  // this should be done in updateIndexSet(i) for all integrators only
  // if a graph has changed
  //updateIndexSet(1);
  //_nsds->topology()->indexSet(1)->update_vertices_indices();
  //_nsds->topology()->indexSet(1)->update_edges_indices();

  //if(_nsds->topology()->hasChanged())
  //{
  //  for(OSNSIterator itOsns = _allNSProblems->begin(); itOsns != _allNSProblems->end(); ++itOsns)
  //  {
  //    (*itOsns)->setHasBeenUpdated(false);
  //  }
  //}

  if (!_allNSProblems->empty())
    computeOneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY);

  DEBUG_EXPR(
    if (_nsds->topology()->indexSet(1)->size() >0)
      (*_allNSProblems)[SICONOS_OSNSP_TS_VELOCITY]->display();
    );

  // update on impulse level
  // * calculate global impulse (p_{k+1}^+)
  // * update velocity (v_{k+1}^+) with OSI in ds->velocity
  // * calculate local gaps (g_{k+1}^+)
  update(1);

  // indexset (I_{k+1}^+) is calculated in Simulation::processEvent
}


void TimeSteppingD1Minus::updateInput(unsigned int level)
{
  DEBUG_BEGIN("TimeSteppingD1Minus::updateInput(unsigned int level)\n");
  OSIIterator itOSI;
  // 1 - compute input (lambda -> r)
  if (!_allNSProblems->empty())
  {
    for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
      (*itOSI)->updateInput(nextTime(),level);
    //_nsds->updateInput(nextTime(),levelInput);
  }
  DEBUG_END("TimeSteppingD1Minus::updateInput(unsigned int level)\n");

}




// void TimeSteppingD1Minus::updateInput(unsigned int level)
// {
//   //  assert(level>=0);

//   double time = nextTime();
//   SP::Topology topology = _nsds->topology();
//   InteractionsIterator it;

//   // // set dynamical systems non-smooth part to zero.
//   // for (OSIIterator itOSI = _allOSI->begin(); itOSI != _allOSI->end(); ++itOSI)
//   // {
//   //   for (DSIterator itDS = (*itOSI)->dynamicalSystems()->begin(); itDS != (*itOSI)->dynamicalSystems()->end(); ++itDS)
//   //   {
//   //     Type::Siconos dsType = Type::value(**itDS);
//   //     if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
//   //       RuntimeException::selfThrow("TimeSteppingD1Minus::updateInput - not implemented for Dynamical system type: " + dsType);
//   //     else
//   //     {
//   //       SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*itDS);
//   //       if (d->p(level)) d->p(level)->zero();
//   //     }
//   //   }
//   // }

//   // Set dynamical systems non-smooth part to zero.
//   reset(level);


//   // we compute input using lambda(level).
//   for (it = topology->interactions()->begin(); it != topology->interactions()->end(); it++)
//   {
//     assert((*it)->lowerLevelForInput() <= level);
//     assert((*it)->upperLevelForInput() >= level);
//     (*it)->computeInput(time, level);
//   }
// }

void TimeSteppingD1Minus::computeResidu()
{
  for (OSIIterator it = _allOSI->begin(); it != _allOSI->end() ; ++it)
    (*it)->computeResidu();
}

void TimeSteppingD1Minus::computeFreeState()
{
  std::for_each(_allOSI->begin(), _allOSI->end(), std11::bind(&OneStepIntegrator::computeFreeState, _1));
}
