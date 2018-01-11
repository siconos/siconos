/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
/*! \file Event.hpp
  General interface for Events
*/

#ifndef Event_H
#define Event_H

#include <cmath>
#include <gmp.h>
#include <cstddef>
#include "SiconosConst.hpp"
#include "SimulationTypeDef.hpp"
#include "SiconosPointers.hpp"
#include "SiconosSerialization.hpp"

// As always, MSVC miss C99
#if defined(_MSC_VER) && _MSC_VER < 1800
extern "C" double rint(double x);
#endif

// tick default value
// it has to be greater than DBL_EPSILON ...
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
  mpz_t _timeOfEvent;

  /** Number of ticks corresponding to a timestep */
  mpz_t _tickIncrement;

  /** Id or type of the Event */
  int _type;

  /** Date of the present event,
   *  represented with a double */
  double _dTime;

  /** confidence interval used to convert double time value to mpz_t
   */
  static double _tick;

  /** has one Event object been instanciated. Use to detect in setTick potentially dangerous cases*/
  static bool _eventCreated;

  /** index for the current Event*/
  unsigned int _k;

  /** TimeDiscretisation for the Event (unused only in the NonSmoothEvent) */
  SP::TimeDiscretisation _td;

  /** For automatic rescheduling */
  bool _reschedule;

  /** Default constructor */
  Event(): _type(0), _dTime(0.0), _k(0), _reschedule(false)
  {
    mpz_init(_timeOfEvent);
    mpz_init(_tickIncrement);
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
   *  \param reschedule set this to true if the event has to be rescheduled
   */
  Event(double time, int newType = 0, bool reschedule = false);

  /** destructor
   */
  virtual ~Event();

  /** get tick value
   *  \return a double
   */
  inline double getTick() const
  {
    return _tick;
  };

  /** set tick value
   *  \param newTick the new tick value
   */
  static void setTick(double newTick);

  /** get the time of the present event (mpz_t format)
   *  \return a mpz_t
   */
  inline const mpz_t * getTimeOfEvent() const
  {
    return &_timeOfEvent ;
  };

  /** get the time of the present event (double format)
   *  \return a double
   */
  inline double getDoubleTimeOfEvent() const
  {
    return _dTime;
  }

  inline void incrementTime(unsigned int step = 1)
  {
    for (unsigned int i = 0; i < step; i++)
      mpz_add(_timeOfEvent, _timeOfEvent, _tickIncrement);
    _dTime = mpz_get_d(_timeOfEvent)*_tick;
  }

  /** set the time of the present event (double format)
   *  \param time the new time
   */
  inline void setTime(double time)
  {
    _dTime = time;
    mpz_set_d(_timeOfEvent, rint(_dTime / _tick));
  };

  /** get a type of the present event
   *  \return an std::string
   */
  inline int getType() const
  {
    return _type ;
  };

  /** set a new type for the present Event
   *  \param newType the new Event type
   */
  inline void setType(int newType)
  {
    _type = newType;
  };

  /** Set the current step k
   * \param newK the new value of _k
   */
  inline void setK(unsigned int newK) { _k = newK; };

  /** Set the TimeDiscretisation
   * \param td a TimeDiscretisation for this Event
   */
  void setTimeDiscretisation(SP::TimeDiscretisation td);

  /** Get the TimeDiscretisation
   * \return the TimeDiscretisation used in this Event
   */
  inline SP::TimeDiscretisation getTimeDiscretisation() const { return _td; };

  /** display Event data
   */
  void display() const ;

  /** virtual function which actions depends on event type
   * \param sim the simulation that owns this Event (through the EventsManager)
   */
  virtual void process(Simulation& sim) = 0;

  /** virtual function which actions depends on event type.
   * The generic implementation present in this object is to increment the
   * TimeDiscretisation and to chamge the time of the current Event 
   \param k meaning depends on the type of event. See derived class.
  */
  virtual void update(unsigned int k = 0);

  inline bool reschedule() const { return _reschedule; };
};
#endif // Event_H
