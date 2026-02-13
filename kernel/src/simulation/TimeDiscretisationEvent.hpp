/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef TIMEDISCREVENT_H
#define TIMEDISCREVENT_H

#include "Event.hpp"

namespace siconos::simulation {

/** Event that corresponds to user-defined time discretisation points
 *
 */
class TimeDiscretisationEvent : public Event {
 private:
  ACCEPT_SERIALIZATION(TimeDiscretisationEvent);

  /** Turn this on to limit memory print of event (no swap of the nsds during process ...*/
  bool noSaveInMemory_{false};

 public:
  /** constructor with time value as a parameter
   *  \param time starting time (a double)
   */
  TimeDiscretisationEvent(double time) : Event(time, EventType::TD){};

  /** destructor
   */
  ~TimeDiscretisationEvent() noexcept = default;

  /**
   *  \param simulation the simulation that owns this Event (through the EventsManager)
   */
  void process(Simulation& simulation);

  /** increment the TimeDiscretisation and to change the time of the Event
   * \param k the next index for this event
   */
  void update(unsigned int k);

  /** If called, this event won't save nsds values (state ...) into memory
      use this at your own risk, many integrators needs previous values
   *  to integrate properly
   */
  void noSaveInMemory() { noSaveInMemory_ = true; };
};

// Register the event into the factory
static EventRegistration<TimeDiscretisationEvent> reg_TD(EventType::TD);

}  // namespace siconos::simulation
#endif  // TimeDiscretisationEvent_H
