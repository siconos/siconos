/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
/*! \file Event.hpp
  General interface for Events
*/

#ifndef Event_H
#define Event_H

#include <cmath>
#include <iostream> // Warning (FP): iostream must be included before gmp
#include <gmp.h>
#include "SiconosConst.hpp"
#include "SimulationTypeDef.hpp"
#include "ControlTypeDef.hpp"
#include "SiconosPointers.hpp"
#include "SiconosSerialization.hpp"

class Simulation;

// tick default value
const double DEFAULT_TICK = 1e-16;

/** Abstract class that represents generic time events.
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
 *
 * Existing types of events:
 *  0 -> undef
 *  1 -> TimeDiscretisation
 *  2 -> NonSmooth
 *  3 -> Sensor
 *  4 -> Observer
 *  5 -> Actuator
 */

class Event
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Event);


  /** Date of the present event,
   *  represented with a mpz_t */
  mpz_t timeOfEvent;

  /** Id or type of the Event */
  int type;

  /** Date of the present event,
   *  represented with a double */
  double dTime;

  /** confidence interval used to convert double time value to mpz_t
   */
  static double tick;

  /** TimeDiscretisation for the Event (unused only in the NonSmoothEvent) */
  SP::TimeDiscretisation _td;

  /** Default constructor */
  Event(): type(0), dTime(0.0)
  {
    mpz_init(timeOfEvent);
  };

  /** copy constructor ; private => no copy nor pass-by-value.
   */
  // Event(const Event&); pb python link

  /** assignment operator private => no assign allowed
   */
//  Event& operator = (const Event&);

public:

  /** constructor with time value and type as input
   *  \param time the starting type (a double)
   *  \param newType the Event type (an int)
   */
  Event(double time, int newType = 0);

  /** destructor
   */
  virtual ~Event();

  /** get tick value
   *  \return a double
   */
  inline double getTick() const
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
  inline double getDoubleTimeOfEvent() const
  {
    return dTime;
  }

  /** set the time of the present event (double format)
   *  \param a double
   */
  inline void setTime(double time)
  {
    dTime = time;
    mpz_clear(timeOfEvent);
    mpz_init_set_d(timeOfEvent, rint(dTime / tick));
  };

  /** get a type of the present event
   *  \return an std::string
   */
  inline int getType() const
  {
    return type ;
  };

  /** set a new type for the present event
   *  \param a std::string
   */
  inline void setType(int newType)
  {
    type = newType;
  };

  /** Set the TimeDiscretisation
   * \param td a TimeDiscretisation for this Event
   */
  inline void setTimeDiscretisation(SP::TimeDiscretisation td)
  {
    _td = td;
  }

  /** display Event data
   */
  void display() const ;

  /** virtual function which actions depends on event type
   * \param sim the simulation that owns this Event (through the EventsManager)
   */
  virtual void process(Simulation& sim) = 0;

  /** virtual function which actions depends on event type.
   * The generic implementation present in this object is to increment the
   * TimeDiscretisation and to chamge the time of the current Event */
  virtual void update();

};
#endif // Event_H
