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
/*! \file SensorEvent.hpp
  Sensor Events
*/
#ifndef SensorEvent_H
#define SensorEvent_H

#include "Event.hpp"

namespace siconos::control {

class Sensor;

/** Events when sensor data capture is done.
 *
 */
class SensorEvent : public siconos::simulation::Event {
 private:
  ACCEPT_SERIALIZATION(SensorEvent);

  using EventType = siconos::simulation::EventType;
  
  /** The sensor linked to the present event */
  std::shared_ptr<Sensor> _sensor{nullptr};

  // /** Default constructor */
  // SensorEvent() : Event(0.0, EventType::Sensor, true){};

 public:
  /** constructor with time value as a parameter
   *  \param time the starting time of the Event
   */
  SensorEvent(double time) : Event(time, EventType::Sensor, true){};

  /** destructor
   */
  ~SensorEvent() noexcept = default;

  /** get the Sensor linked to this Event
   *  \return a pointer to the Sensor
   */
  inline std::shared_ptr<Sensor> sensor() const { return _sensor; };

  /** set the Sensor linked to this Event
   *  \param newSensor the std::shared_ptr<Sensor>
   */
  void setSensorPtr(std::shared_ptr<Sensor> newSensor) { _sensor = newSensor; };

  /** Call the capture method of the linked Sensor
   *  \param sim a std::shared_ptr<siconos::simulation::Simulation> (ignored).
   */
  void process(siconos::simulation::Simulation& sim);
};


}  // namespace siconos::control

namespace siconos::simulation{

    // Register the event into the factory
  static EventRegistration<siconos::control::SensorEvent> reg_SE(EventType::Sensor);
}

#endif  // SensorEvent_H
