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
/*! \file ActuatorEvent.hpp
  \brief Actuator Events
*/
#ifndef ActuatorEvent_H
#define ActuatorEvent_H

#include "Event.hpp"
#include "Actuator.hpp"

/** Events when sensor data capture is done.
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 09, 2007
 *
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
  ActuatorEvent(): Event(0.0, 4) {};

public:

  /** constructor with time value as a parameter
   *  \param time the time of the Event
   *  \param name the type of Event
   */
  ActuatorEvent(double time, int name): Event(time, name) {};

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
  void process(SP::Simulation sim);

  /** Increment the time discretisation of the linked Actuator
  */
  void update();
};

#endif // ActuatorEvent_H
