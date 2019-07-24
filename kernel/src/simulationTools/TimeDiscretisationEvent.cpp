/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "TimeDiscretisationEvent.hpp"
#include "EventFactory.hpp"
#include "Simulation.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
using namespace EventFactory;

// Default constructor
TimeDiscretisationEvent::TimeDiscretisationEvent(): Event(0.0, TD_EVENT)
{}

TimeDiscretisationEvent::TimeDiscretisationEvent(double time, int notUsed): Event(time, TD_EVENT)
{}

TimeDiscretisationEvent::~TimeDiscretisationEvent()
{}

void TimeDiscretisationEvent::process(Simulation& simulation)
{
  // Update y[i] values in Interactions with new DS states.
  //simulation->updateOutput(0, 1);
  // Save state(s) in Memories (DS and Interactions, through OSI and OSNS).

  simulation.nonSmoothDynamicalSystem()->swapInMemory();  // To save pre-impact values
  simulation.nonSmoothDynamicalSystem()->pushInteractionsInMemory();  // To save pre-impact values

}

void TimeDiscretisationEvent::update(unsigned int k)
{
  assert(k > _k && "TimeDiscretisationEvent::update - next step has to be greater than the current one");
  if (_td) // if no TimeDiscretisation, then do nothing
  {
    if (_td->hGmp())
      incrementTime(k-_k);
    else
      setTime(_td->getTk(k));

    _k = k;
  }
}

AUTO_REGISTER_EVENT(TD_EVENT, TimeDiscretisationEvent)
