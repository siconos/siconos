/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
/*! \file
 Non-Smooth Events
*/
#ifndef NONSMOOTHEVENT_H
#define NONSMOOTHEVENT_H

/** Events due to non smooth behavior (contact occurence...)
 *
 * Those events are detected during Simulation process (integration of the smooth part with a roots-finding algorithm)
 * and scheduled into the EventsManager.
 *
 */

#include "Event.hpp"

class NonSmoothEvent : public Event
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NonSmoothEvent);


  /** Default constructor */
  NonSmoothEvent();

public:

  /** constructor with time value as a parameter
  *  \param time the time of the first event (a double)
  *  \param notUsed unused parameter (an int)
  */
  NonSmoothEvent(double time, int notUsed);

  /** destructor
  */
  ~NonSmoothEvent();

  /** OSNS solving and IndexSets updating
  *  \param simulation the simulation that owns this Event (through the EventsManager)
  */
  void process(Simulation& simulation);
};

#endif // NonSmoothEvent_H
