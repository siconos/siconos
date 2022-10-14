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
#include "NonSmoothEvent.hpp"

#include "Event.hpp"
#include "EventDriven.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SiconosException.hpp"
#include "SimulationGraphs.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #include "siconos_debug.h"

void siconos::simulation::NonSmoothEvent::process(Simulation& simulation)
{
  auto eventDriven = dynamic_cast<EventDriven*>(&simulation);
  if (!eventDriven)
    THROW_EXCEPTION("NonSmoothEvent::process failed; Simulation is not of EventDriven type.");

  if (!(simulation.oneStepNSProblems()->empty())) {
    // Compute y[0], y[1] and update index sets. => already done
    // during advance to event ...

    // Get the required index sets ...
    auto indexSet0 = simulation.indexSet(0);

    // Update all the index sets ...
    eventDriven->updateIndexSets();
    auto indexSet1 = simulation.indexSet(1);
    auto indexSet2 = simulation.indexSet(2);
    bool found = true;
    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
      found = indexSet2->is_vertex(indexSet1->bundle(*ui));
      if (!found) break;
    }
    // ---> solve impact LCP if IndexSet[1]\IndexSet[2] is not empty.
    if (!found) {
      // For Event-Driven algo., memories vectors are of size 2
      // (ie 2 interactionBlocks).  First interactionBlock (pos 0, last
      // in) for post-event values and last interactionBlock (pos 1,
      // first in) for pre-event values.

      simulation.nonSmoothDynamicalSystem()->swapInMemory();  // To save pre-impact values
      simulation.nonSmoothDynamicalSystem()
          ->pushInteractionsInMemory();  // To save pre-impact values

      // solve the LCP-impact => y[1],lambda[1]
      eventDriven->computeOneStepNSProblem(SICONOS_OSNSP_ED_IMPACT);  // solveLCPImpact();
      // compute p[1], post-impact velocity, y[1] and indexSet[2]
      simulation.update(1);
      // Update the corresponding index set ...
      eventDriven->updateIndexSets();
    }

    //---> solve acceleration LCP if IndexSet[2] is not empty
    if (indexSet2->size() > 0) {
      // solve LCP-acceleration
      eventDriven->computeOneStepNSProblem(
          SICONOS_OSNSP_ED_SMOOTH_ACC);  // solveLCPAcceleration();
      // update input of level 2, acceleration and output of level 2
      simulation.update(2);
      // for all index in IndexSets[2], update the index set according to y[2] and/or lambda[2]
      // sign.
      eventDriven->updateIndexSetsWithDoubleCondition();
    }

    // Save results in memory
    simulation.nonSmoothDynamicalSystem()->swapInMemory();  // To save pre-impact values
    simulation.nonSmoothDynamicalSystem()
        ->pushInteractionsInMemory();  // To save pre-impact values
  }
}
