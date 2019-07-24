/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "SiconosControlFwd.hpp"
#include "ControlTypeDef.hpp"

/** Events when sensor data capture is done.
 *
 */
class ActuatorEvent : public Event
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(ActuatorEvent);

  /** The actuator linked to the present event */
  SP::Actuator _actuator;

  /** Default constructor */
  ActuatorEvent(): Event(0.0, ACTUATOR_EVENT, true) {};

public:

  /** constructor with time value as a parameter
   *  \param time the time of the Event
   *  \param name the type of Event
   */
  ActuatorEvent(double time, int name): Event(time, name, true) {};

  /** destructor
   */
  ~ActuatorEvent() {};

  /** get the Actuator linked to this Event
   *  \return a pointer to Actuator
   */
  inline SP::Actuator actuator() const
  {
    return _actuator;
  };

  /** set the Actuator linked to this ActuatorEvent
   *  \param newActuator the Actuator associated with this Event.
   */
  void setActuatorPtr(SP::Actuator newActuator)
  {
    _actuator = newActuator;
  };

  /** Call the actuate method of the Actuator
   *  \param sim ignored argument.
   */
  void process(Simulation& sim);

};

#endif // ActuatorEvent_H
