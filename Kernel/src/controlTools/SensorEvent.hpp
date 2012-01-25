/* Siconos-Kernel, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
/*! \file SensorEvent.hpp
  Sensor Events
*/
#ifndef SensorEvent_H
#define SensorEvent_H

#include "Event.hpp"
#include "Sensor.hpp"

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
  SensorEvent(): Event(0.0, 3) {};

public:

  /** constructor with time value as a parameter
   *  \param time the starting time of the Event
   *  \param name the type of the Event
   */
  SensorEvent(double time, int name): Event(time, name) {};

  /** destructor
   */
  ~SensorEvent() {};

  /** get the Sensor linked to this Event
   *  \return a SP::Sensor to the Sensor
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
  void process(SP::Simulation sim);

  /** Increment time of the present event according to
      the time discretisation of the linked Actuator
  */
  void update();
};

#endif // SensorEvent_H
