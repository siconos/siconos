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
/*! \file
Time Discretisation Events
*/
#ifndef TIMEDISCRETISATIONEVENTNOSAVEINMEMORY_H
#define TIMEDISCRETISATIONEVENTNOSAVEINMEMORY_H

#include "Event.hpp"

/** Event that corresponds to user-defined time discretisation points
 *  This Event does not automatically save in memory some variables.
 *  Use it at your own risk
 *
 */
class TimeDiscretisationEventNoSaveInMemory : public Event
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(TimeDiscretisationEventNoSaveInMemory);


  /** Default constructor */
  TimeDiscretisationEventNoSaveInMemory();

public:

  /** constructor with time value as a parameter
  *  \param time starting time (a double)
  *  \param notUsed unused int
  */
  TimeDiscretisationEventNoSaveInMemory(double time, int notUsed);

  /** destructor
  */
  ~TimeDiscretisationEventNoSaveInMemory();

  /** increment the TimeDiscretisation and to change the time of the Event
   * \param k the next index for this event
   */
  void update(unsigned int k);

  /**
  *  \param simulation the simulation that owns this Event (through the EventsManager)
  */
  void process(Simulation& simulation);
};

#endif // TimeDiscretisationEventNoSaveInMemory_H
