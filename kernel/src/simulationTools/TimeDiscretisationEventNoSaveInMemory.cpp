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
#include "TimeDiscretisationEventNoSaveInMemory.hpp"
#include "TimeDiscretisation.hpp"
#include "EventFactory.hpp"

using namespace EventFactory;

// Default constructor
TimeDiscretisationEventNoSaveInMemory::TimeDiscretisationEventNoSaveInMemory(): Event(0.0, TD_EVENT)
{}

TimeDiscretisationEventNoSaveInMemory::TimeDiscretisationEventNoSaveInMemory(double time, int notUsed): Event(time, TD_EVENT)
{}

TimeDiscretisationEventNoSaveInMemory::~TimeDiscretisationEventNoSaveInMemory()
{}

void TimeDiscretisationEventNoSaveInMemory::process(Simulation& simulation)
{}

void TimeDiscretisationEventNoSaveInMemory::update(unsigned int k)
{
  assert(k > _k && "TimeDiscretisationEvent::update - next step has to be greater than the current one");
  if(_td)  // if no TimeDiscretisation, then do nothing
  {
    if(_td->hGmp())
      incrementTime(k-_k);
    else
      setTime(_td->getTk(k));

    _k = k;
  }
}
//AUTO_REGISTER_EVENT(TD_EVENT, TimeDiscretisationEventNoSaveInMemory)
