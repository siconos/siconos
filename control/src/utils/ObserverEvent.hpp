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
  ObserverEvent(double time) : Event(time, EventType::Observer, true){};

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
  void setObserverPtr(std::shared_ptr<Observer> newObserver)
  {
    _observer = newObserver;
  };

  /** Call the capture method of the linked Observer
   *  \param sim a std::shared_ptr<siconos::simulation::Simulation> (ignored).
   */
  void process(siconos::simulation::Simulation& sim);
};


}  // namespace siconos::control

namespace siconos::simulation{
// Register the event into the factory
  static EventRegistration<siconos::control::ObserverEvent> reg_OB(EventType::Observer);


}

#endif  // ObserverEvent_H
