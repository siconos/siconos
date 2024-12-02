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
/*! \file ActuatorEvent.hpp
  \brief Actuator Events
*/
#ifndef ActuatorEvent_H
#define ActuatorEvent_H

#include "Event.hpp"

namespace siconos::control {

class Actuator;

/** Events when sensor data capture is done.
 *
 */
class ActuatorEvent : public siconos::simulation::Event {
 private:
  ACCEPT_SERIALIZATION(ActuatorEvent);

  using EventType = siconos::simulation::EventType;

  /** The actuator linked to the present event */
  std::shared_ptr<Actuator> _actuator{nullptr};

  // /** Default constructor */
  // ActuatorEvent(): Event(0.0, EventType::Actuator, true) {};

 public:
  /** constructor with time value as a parameter
   *  \param time the time of the Event
   */
  ActuatorEvent(double time) : Event(time, EventType::Actuator, true){};

  /** destructor
   */
  ~ActuatorEvent() noexcept = default;

  /** get the Actuator linked to this Event
   *  \return a pointer to Actuator
   */
  inline std::shared_ptr<Actuator> actuator() const { return _actuator; };

  /** set the Actuator linked to this ActuatorEvent
   *  \param newActuator the Actuator associated with this Event.
   */
  void setActuatorPtr(std::shared_ptr<Actuator> newActuator) { _actuator = newActuator; };

  /** Call the actuate method of the Actuator
   *  \param sim ignored argument.
   */
  void process(siconos::simulation::Simulation& sim);

};
}  // namespace siconos::control

namespace siconos::simulation{
    // Register the event into the factory
  static EventRegistration<siconos::control::ActuatorEvent> reg_AC(EventType::Actuator);

}
#endif  // ActuatorEvent_H
