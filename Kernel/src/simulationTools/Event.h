/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include<gmp.h>
#include<string>
#include<iostream>

class Simulation;

const std::string DEFAULT_EVENT_TYPE = "undefined";

// tick default value
const double DEFAULT_TICK = 1e-16;

/** virtual class that represents generic time events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 21, 2006
 *
 *  This base class simply records the time at which the event will take place. A pure virtual function named process
 *  will be invoked to execute the event.
 *  The time is represented with a mpz_t, from gmp library. See http://gmplib.org.
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
   *  represented with a mpz_t */
  mpz_t timeOfEvent;

  /** Id or type of the Event */
  const std::string type;

  /** Date of the present event,
   *  represented with a double */
  const double dTime;

  /** confidence interval used to convert double time value to mpz_t
   */
  static double tick;

  /** Default constructor */
  Event();

  /** copy constructor ; private => no copy nor pass-by-value.
  */
  Event(const Event&);

public:

  /** constructor with time value and type as input
  *  \param double
  *  \param a string
  */
  Event(double, const std::string & = DEFAULT_EVENT_TYPE);

  /** destructor
   */
  virtual ~Event();

  /** get tick value
  *  \return a double
  */
  inline const double getTick() const
  {
    return tick ;
  };

  /** set tick value
  *  \param a double
  */
  inline void setTick(double newTick)
  {
    std::cout << "Warning: you change tick value for EventsManager -> a new initialization of the object is required. " << std::endl;
    tick = newTick;
  };

  /** get the time of the present event (mpz_t format)
   *  \return a mpz_t
   */
  inline const mpz_t * getTimeOfEvent() const
  {
    return &timeOfEvent ;
  };

  /** get the time of the present event (double format)
   *  \return a double
   */
  inline const double getDoubleTimeOfEvent() const
  {
    return dTime;
  }

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
