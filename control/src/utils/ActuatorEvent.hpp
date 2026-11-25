/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
  ActuatorEvent(double time) : Event(time, EventType::Actuator, true) {};

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

namespace siconos::simulation {
// Register the event into the factory
static EventRegistration<siconos::control::ActuatorEvent> reg_AC(EventType::Actuator);

}  // namespace siconos::simulation
#endif  // ActuatorEvent_H
