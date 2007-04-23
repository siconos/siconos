/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
/*! \file SensorEvent.h
Sensor Events
*/
#ifndef SensorEvent_H
#define SensorEvent_H

#include "Event.h"

class Sensor;

/** Events when sensor data capture is done.
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) February 01, 2007
 *
 *
 */
class SensorEvent : public Event
{

private:

  /** The sensor linked to the present event */
  Sensor* sensor;

  /** Default constructor */
  SensorEvent();

public:

  /** constructor with time value as a parameter
   *  \param a double
   *  \param a string, type of Event
   */
  SensorEvent(double, const std::string&);

  /** destructor
  */
  ~SensorEvent();

  /** get the Sensor linked to this Event
   *  \return a pointer to Sensor
   */
  inline Sensor* getSensorPtr() const
  {
    return sensor;
  };

  /** set the TimeDiscretisation linked to this Sensor
   *  \param a pointer to TimeDiscretisation.
   */
  void setSensorPtr(Sensor * newSensor)
  {
    sensor = newSensor;
  };

  /**
   *  \param Simulation*, the simulation that owns this Event (through the EventsManager)
   */
  void process(Simulation*);
};

#endif // SensorEvent_H
