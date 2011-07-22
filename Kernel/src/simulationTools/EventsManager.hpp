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
/*! \file EventsManager.h
  \brief Management of a Set (STL) of Events for the Simulation process
*/

#ifndef EVENTSMANAGER_H
#define EVENTSMANAGER_H

#include "SiconosConst.hpp"
#include "Event.hpp"
#include <gmp.h>
#include<iostream>
#include<set>

class Simulation;
class TimeDiscretisation;
const unsigned long int GAPLIMIT_DEFAULT = 100;

/** Structure used for Events sorting. The time of occurence (in mpz_t format!!!)  is used to compare two Events. */
struct compareEvent
{
  bool operator()(const SP::Event e1, const SP::Event e2) const
  {
    const mpz_t *t1 = e1->getTimeOfEvent();
    const mpz_t *t2 = e2->getTimeOfEvent();
    int res = mpz_cmp(*t1, *t2); // res>0 if t1>t2, 0 if t1=t2 else res<0.
    return (res < 0);
  }
};

/** set of events, with an ordering based on Event time value (mpz_t) to compare Events
 *  A stl container of type "multiset" is used at the time
 *  \warning This may be not the best choice => review all possibi lities */
typedef std::multiset<SP::Event, compareEvent > EventsContainer; // sort in a chronological way

/** Iterator through a set of Events */
typedef EventsContainer::iterator EventsContainerIterator;

/** Tools to handle a set of Events for the Simulation

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) February 22, 2006

   The EventsManager handles a set of events (from user time-discretisation, sensors, non-smooth ...),
   and is supposed to provide to the simulation the values of "current" and "next" events to define
   the time-integration interval.

   Events:
   - currentEvent: starting time for integration. Initialized with t0 of the simulation time-discretisation.
   - ETD: corresponds to the instant t[k+1] of the Simulation (user) TimeDiscretisation.
   - ENonSmooth: for EventDriven simulation only. Present only if one or more non-smooth events have been
   detected between currentEvent and ETD.
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


  /** list of future, not processed, events.
   *  This list is not fixed and can be updated at any time
   *  depending on the simulation, user add ...
   *  The first event of this set is currentEvent, and the second is nextEvent.
   * ("first" and "second" defined according to event comparison operator)
   */
  EventsContainer _allEvents;

  /** Pointer to currentEvent, ie the simulation starting point.
   * It correponds to the first object in allEvents.
   */
  SP::Event _currentEvent;

  /** Pointer to nextEvent, ie the simulation ending point.
   * It correponds to the event following currentEvent and so
   * to the second object in allEvents.
   */
  SP::Event _nextEvent;

  /** Event which corresponds to time tk+h of the simulation time discretisation */
  SP::Event _ETD;

  /** Non Smooth Event: corresponds to the last non-smooth event detected*/
  SP::Event _ENonSmooth;

  /* link to the simulation that owns this manager*/
  SP::Simulation _simulation;

  /** bool to check if some NonSmooth events has already been inserted */
  bool _hasNS;

  /** bool to check if a ControlManager is associated to the simulation */
  bool _hasCM;

  /** unsigned long int variable to check if two events are so close or not */
  static unsigned long int GapLimit2Events;

  /** Creates and adds a new Event in the allEvents list
   *  \return false if Event already exists
   */
  bool createAndInsertEvent(int, double);

  /** copy constructor => private: no copy nor pass-by-value.
   *  \param the eventsManager to be copied
   */
  EventsManager(const EventsManager&);

public:

  /**  default constructor
   *  \param the simulation that owns this manager
   */
  EventsManager();

  /** destructor
   */
  ~EventsManager() {};

  /** set the gap limit between two events */
  inline void setGapLimitEvents(unsigned long int var)
  {
    GapLimit2Events = var;
  };
  /**Set the gap limit between two events  */
  inline unsigned long int getGapLimitEvents()
  {
    return GapLimit2Events;
  };

  /** initialize current, next events and the events stack.
      \param simulation, owner of the present eventsManager.
   */
  void initialize(SP::Simulation);

  /** Run process of all events simultaneous to currentEvent + synchr. with actuators/sensors
      Must be run at the end of simulation->initialize()
   */
  void preUpdate();

  /** add a set of existing Events into allEvents list
   *  \param an EventsContainer
   */
  void insertEvents(const EventsContainer&);

  /** get the list of all Events in the set
   *  \return a set of Events*
   */
  inline const EventsContainer events() const
  {
    return _allEvents ;
  };

  /** get the event that occurs at time inputTime
   *  \param a mpz_t
   *  \return a pointer to Event
   */
  SP::Event event(const mpz_t& inputTime) const;

  /** get the current event
   *  \return a pointer to Event
   */
  inline SP::Event currentEvent() const
  {
    return _currentEvent;
  };

  /** get the next event to be processed.
   *  \return a pointer to Event
   */
  inline SP::Event nextEvent() const
  {
    return _nextEvent;
  };

  /** get the event following inputEvent ("following" defined with
   *   operator(s) comparison of events)
   *  \param a pointer to Event
   *  \return a pointer to Events
   */
  SP::Event followingEvent(SP::Event) const;

  /** get the event that follows the event at time inputTime
   *   ("following" defined with operator(s) comparison of events)
   *  \param a mpz_t
   *  \return a pointer to Event
   */
  SP::Event followingEvent(const mpz_t& inputTime) const;

  /** get the Simulation
   *  \return a pointer to Simulation
   */
  inline SP::Simulation simulation() const
  {
    return _simulation;
  }

  /** set the Simulation of the OneStepNSProblem
   *  \param: a pointer on Simulation
   */
  inline void setSimulationPtr(SP::Simulation str)
  {
    _simulation = str;
  }

  /** check if event is present in allEvents list
   *  \return a bool
   */
  bool hasEvent(SP::Event) const ;

  /** check if some events remain in allEvents list
   *  \return a bool
   */
  bool hasNextEvent() const ;

  /** get the time (double format) of an event
   *  \return a double
   */
  double getTimeOfEvent(SP::Event) const;

  /** get the time of current event, in double format
   *  \return a double
   */
  double startingTime() const ;

  /** get the time of next event, in double format
   *  \return a double
   */
  double nextTime() const ;

  /** display EventsManager data
   */
  void display() const ;

  /** Insert an existing event into the set
      \param Event*, the event to be nserted
   */
  void insertEvent(SP::Event);

  /** Update nextEvent (useful in case of a new insertion)
   */
  void update();

  /** add a new Event in the allEvents list and update nextEvent value
   * \param the time (double format) of occurence of the event
   */
  void scheduleNonSmoothEvent(double);

  /** remove an Event from the unProcessed events list
   */
  void removeEvent(SP::Event);

  /** update current and next event positions and run process functions of current event
   */
  void processEvents();

  /** processEvents when neither non-smooth events nor control-related events are present.
   */
  void OptimizedProcessEvents();

  /** processEvents for general case
   */
  void GeneralProcessEvents();

  /** used when event's time has changed to resort it properly in allEvents set
      \param the event to be sorted
   */
  void SortEvent(SP::Event);
};

DEFINE_SPTR(EventsManager);

#endif // EventsManager_H
