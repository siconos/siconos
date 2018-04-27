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
/*! \file SensorEvent.hpp
  Sensor Events
*/
#ifndef SensorEvent_H
#define SensorEvent_H

#include "Event.hpp"
#include "SiconosControlFwd.hpp"
#include "ControlTypeDef.hpp"

/** Events when sensor data capture is done.
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 01, 2007
 *
 *
 */
class SensorEvent : public Event
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SensorEvent);


  /** The sensor linked to the present event */
  SP::Sensor _sensor;

  /** Default constructor */
  SensorEvent(): Event(0.0, SENSOR_EVENT, true) {};

public:

  /** constructor with time value as a parameter
   *  \param time the starting time of the Event
   *  \param name the type of the Event
   */
  SensorEvent(double time, int name): Event(time, name, true) {};

  /** destructor
   */
  ~SensorEvent() {};

  /** get the Sensor linked to this Event
   *  \return a pointer to the Sensor
   */
  inline SP::Sensor sensor() const
  {
    return _sensor;
  };

  /** set the Sensor linked to this Event
   *  \param newSensor the SP::Sensor
   */
  void setSensorPtr(SP::Sensor newSensor)
  {
    _sensor = newSensor;
  };

  /** Call the capture method of the linked Sensor
   *  \param sim a SP::Simulation (ignored).
   */
  void process(Simulation& sim);

};

#endif // SensorEvent_H
