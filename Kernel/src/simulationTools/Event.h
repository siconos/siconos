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
/*! \file
 General Event
*/

#ifndef EVENT_H
#define EVENT_H

#include "SiconosConst.h"
#include<string>
#include<iostream>

class Simulation;

const unsigned long int DEFAULT_EVENT_TIME = 0;
const std::string DEFAULT_EVENT_TYPE = "undefined";

/** virtual class that represents generic time events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) February 21, 2006
 *
 *  This base class simply records the time at which the event will take place. A pure virtual function named process
 *  will be invoked to execute the event.
 *  The time is represented with an unsigned long int.
 *
 * Derived classes:
 * - TimeDiscretisationEvent: events that corresponds to user-defined time-discretisation points
 * - NonSmoothEvent: specific events, detected during simulation, when constraints are violated (thanks to roots-finding algorithm)
 * - SensorEvent: event dedicated to data capture through user-defined sensors.
 *
 */

class Event
{
protected:

  /** Date of the present event,
   *  represented with a long int */
  const unsigned long int timeOfEvent;

  /** Id or type of the Event */
  const std::string type;

  /** Default constructor */
  Event();

  /** copy constructor ; private => no copy nor pass-by-value.
  */
  Event(const Event&);

public:

  /** constructor with time value and type as input
  *  \param an unsigned int
  *  \param a string
  */
  Event(unsigned long int, const std::string & = DEFAULT_EVENT_TYPE);

  /** destructor
   */
  virtual ~Event();

  /** get the time of the present event (unsigned int format)
   *  \return an unsigned long int
   */
  inline const unsigned long int getTimeOfEvent() const
  {
    return timeOfEvent ;
  };

  /** get a type of the present event
   *  \return an std::string
   */
  inline const std::string getType() const
  {
    return type ;
  };

  /** display Event data
   */
  void display() const ;

  /** virtual function which actions depends on event type
   * \param Simulation*, the simulation that owns this Event (through the EventsManager)
   */
  virtual void process(Simulation*) = 0;
};

#endif // Event_H
