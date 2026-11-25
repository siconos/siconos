/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
  TimeDiscretisationEvent(double time) : Event(time, EventType::TD) {};

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
