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
/*! \file EventsManager.hpp
  \brief Management of a Set (STL) of Events for the Simulation process
*/

#ifndef EventsManager_H
#define EventsManager_H

#include "SiconosConst.hpp"
#include "Event.hpp"
#include <gmp.h>
#include <iostream>
#include <set>

class Simulation;
class TimeDiscretisation;
const unsigned long int GAPLIMIT_DEFAULT = 100;

/** set of events, with an ordering based on Event time value (mpz_t) to compare Events */
typedef std::vector<SP::Event> EventsContainer; // Event are already sorted

/** Tools to handle a set of Events for the Simulation

   \author SICONOS Development Team - copyright INRIA
   \version 3.6.0.
   \date (Creation) December 22, 2012

   The EventsManager handles a set of events (from user time-discretisation, sensors, non-smooth ...),
   and is supposed to provide to the simulation the values of "current" and "next" events to define
   the time-integration interval.

   Events:
   - currentEvent: starting time for integration. Initialized with t0 of the simulation time-discretisation.
   - ETD: corresponds to the instant t[k+1] of the Simulation (user) TimeDiscretisation.
   - ENonSmooth: for EventDriven simulation only. Present only if one or more non-smooth events have been
   detected between currentEvent and the next event.
   - Sensor or Actuators Events. To each Sensor or Actuator declared in the ControlManager, corresponds
   an Event in the manager. When this event is processed, its time value is increased to next instant
   in the time-discretisation of the sensor/actuator.

   Examples:
   - for a TimeStepping, with one Sensor, the EventsManager looks like:\n
   {currentEvent(tk), ESensor(tsensor), ETD(tk+1)}
   - for an EventDriven, with one actuator, an non-smooth event detected at time tns: \n
   {currentEvent(tk), EActuator(tact), ENonSmooth(tns), ETD(tk+1)}.

   After each process, the time values of each event are updated and nextEvent points to the first event after currentEvent.

   \section EMMfunc Main functions
   - initialize(): process all events which have the same time as currentEvent
   - processEvents(): process all events simultaneous to nextEvent, increment them to next step, update index sets,
   increment currentEvent.


*/
class EventsManager
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(EventsManager);

  /** list of events
   * The first element is the last processed event. All the others
   * are unprocessed events.
   * This list is not fixed and can be updated at any time
   * depending on the simulation, user add ...
   */
  EventsContainer _events;

  /** Placeholder for the non smooth event */
  SP::Event _eNonSmooth;

  /** unsigned long int variable to check if two events are too close */
  static unsigned long int _GapLimit2Events;

  /** Insert an event in the event stack
   * \param e the event to insert
   * \return the position of the inserted event in the stack
   */
  unsigned int insertEv(SP::Event e);

  /** Update the set of events
   */
  void update(Simulation& sim);

  /** copy constructor => private: no copy nor pass-by-value.
   *  \param the eventsManager to be copied
   */
//  EventsManager(const EventsManager&);

  /** default constructor */
  EventsManager() {};

public:

  /**  default constructor
   * \param sim the Simulation that owns this EventsManager
   */
  EventsManager(Simulation& sim);

  /** destructor
   */
  virtual ~EventsManager() {};

  /** Set the gap limit between two events
   * \param var the new _GapLimit2Events
   */
  inline void setGapLimitEvents(unsigned long int var)
  {
    _GapLimit2Events = var;
  };

  /** Get the gap limit between two events  */
  inline unsigned long int getGapLimitEvents() const
  {
    return _GapLimit2Events;
  };

  /** initialize current, next events and the events stack.
      \param sim the simulation that owns the present eventsManager
   */
//  void initialize(const Simulation& sim);

  /** Change TimeDiscretisationEvent to TimeDiscretisationEventNoSaveInMemory
   * \warning use this at your own risk, many integrators needs previous values
   * to integrate properly
   * \param sim the Simulation that owns this EventsManager
   */
  void noSaveInMemory(const Simulation& sim);

  /** get the current event
   *  \return a pointer to Event
   */
  inline SP::Event currentEvent() const
  {
    return _events[0];
  };

  /** get the next event to be processed.
   *  \return a pointer to Event
   */
  inline SP::Event nextEvent() const
  {
    return _events[1];
  };

  /** return all the events
   * \return a reference to the events set
   */
  inline EventsContainer& events()
  {
    return _events;
  };

  /** check if there are some unprocessed events
   *  \return true if there are unprocessed events
   */
  inline bool hasNextEvent() const
  {
    return _events.size() > 1;
  };

  /** get the time of current event, in double format
   *  \return the time of the last processed events
   */
  double startingTime() const ;

  /** get the time of next event, in double format
   *  \return the time of the next events
   */
  double nextTime() const ;

  /** display EventsManager data
   */
  void display() const ;

  /** add a new Event in the allEvents list and update nextEvent value
   * \param sim the simulation that owns this EventsManager
   * \param time the time (double format) of occurence of the event
   * \param yes_update indicator to update or not the next event (default value is true)
   */
  void scheduleNonSmoothEvent(Simulation& sim, double time, bool yes_update = true);

  /** Process the next event, update the indexSets if necessary
   * \param sim the simulation that owns this EventsManager
   */
  void processEvents(Simulation& sim);

  /** function to be called once after initialization
   * \param sim the simulation that owns this EventsManager
   */
  void preUpdate(Simulation& sim);

  /** insert an event of a certain type. The event is created on the fly.
   * \param type the type of the event
   * \param time the time of the event
   * \return a reference to the Event
   */
  Event& insertEvent(const int type, const double& time);

  /** insert an event of a certain type. The event is created on the fly,
   * and the SP::TimeDiscretisation given in argument is stored inside
   * \param type the type of the event
   * \param td a TimeDiscretisation for the Event
   * \return a reference to the Event
   */
  Event& insertEvent(const int type, SP::TimeDiscretisation td);

};


DEFINE_SPTR(EventsManager)

#endif // EventsManager_H
