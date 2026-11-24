// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file ObserverEvent.hpp
  Observer Events
*/
#ifndef ObserverEvent_H
#define ObserverEvent_H

#include "Event.hpp"
namespace siconos::control {
class Observer;

/** Events when the observer updates the state estimate
 *
 *
 */
class ObserverEvent : public siconos::simulation::Event {
 private:
  ACCEPT_SERIALIZATION(ObserverEvent);

  using EventType = siconos::simulation::EventType;

  /** The observer linked to the present event */
  std::shared_ptr<Observer> _observer{nullptr};

  // /** Default constructor */
  // ObserverEvent(): Event(0.0, EventType::Observer, true) {};

 public:
  /** constructor with time value as a parameter
   *  \param time the starting time of the Event
   */
  ObserverEvent(double time) : Event(time, EventType::Observer, true) {};

  /** destructor
   */
  ~ObserverEvent() noexcept = default;

  /** get the Observer linked to this Event
   *  \return a std::shared_ptr<Observer> to the Observer
   */
  inline std::shared_ptr<Observer> observer() const { return _observer; };

  /** set the Observer linked to this Event
   *  \param newObserver the std::shared_ptr<Observer>
   */
  void setObserverPtr(std::shared_ptr<Observer> newObserver) { _observer = newObserver; };

  /** Call the capture method of the linked Observer
   *  \param sim a std::shared_ptr<siconos::simulation::Simulation> (ignored).
   */
  void process(siconos::simulation::Simulation& sim);
};

}  // namespace siconos::control

namespace siconos::simulation {
// Register the event into the factory
static EventRegistration<siconos::control::ObserverEvent> reg_OB(EventType::Observer);

}  // namespace siconos::simulation

#endif  // ObserverEvent_H
