// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file
 Non-Smooth Events
*/
#ifndef NONSMOOTHEVENT_H
#define NONSMOOTHEVENT_H

#include "Event.hpp"

namespace siconos::simulation {

/** Events due to non smooth behavior (contact occurence...)
 *
 * Those events are detected during Simulation process (integration of the smooth part with a
 * roots-finding algorithm) and scheduled into the EventsManager.
 *
 */
class NonSmoothEvent : public Event {
 private:
  ACCEPT_SERIALIZATION(NonSmoothEvent);

 public:
  /** constructor with time value as a parameter
   *  \param time the time of the first event (a double)
   */
  NonSmoothEvent(double time) : Event(time, EventType::NS) {};

  /** destructor
   */
  ~NonSmoothEvent() noexcept = default;

  /** OSNS solving and IndexSets updating
   *  \param simulation the simulation that owns this Event (through the EventsManager)
   */
  void process(Simulation& simulation);
};

// Register the event into the factory
static EventRegistration<NonSmoothEvent> reg_NS(EventType::NS);

}  // namespace siconos::simulation
#endif  // NonSmoothEvent_H
