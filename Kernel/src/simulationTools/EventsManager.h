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
 Management of a Set (STL) of Events
*/

#ifndef EVENTSMANAGER_H
#define EVENTSMANAGER_H

#include "SiconosConst.h"
#include "RuntimeCmp.h"
#include<iostream>
#include<set>

class Simulation;
class Event;
class TimeDiscretisation;

// tick default value
const double DEFAULT_TICK = 1e-07;

/** set of events, with a RuntimeCmp based on Event time value (unsigned int) to compare Events
 *  A stl container of type "multiset" is used at the time
 *  \Warning This may be not the best choice => review all possibi lities */
typedef std::multiset<Event*, RuntimeCmp<Event> > EventsContainer; // sort in a chronological way

/** Iterator through a set of Events */
typedef EventsContainer::iterator EventsContainerIterator;

/** Tools to handle a set of Events
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) February 22, 2006
 *
 * Some rules:
 *     - Events can be removed or insert anywhere in the list
 *     - An event can not be modified (but removed!) => no setEvent(Event) function
 *     - All Events are created inside this class
 *       That means user can not create an Event and then insert it. He is forced to
 *       scheduleEvent method.
 *     - Tick must be set at construction, and can not be change after.
 *
 */
class EventsManager
{
protected:

  /** list of events already processed.
   *  At the end of the process, currentEvent is inserted in this set.
   */
  EventsContainer pastEvents;

  /** list of future, not processed, events.
   *  This list is not fixed and can be updated at any time
   *  depending on the simulation, user add ...
   *  The first event of this set is currentEvent, and the second is nextEvent.
   * ("first" and "second" defined according to event comparison operator)
   */
  EventsContainer unProcessedEvents;

  /** Pointer to currentEvent, ie the simulation starting point.
    * It correponds to the first object in unProcessedEvents.
    */
  Event * currentEvent;

  /** Pointer to nextEvent, ie the simulation ending point.
    * It correponds to the event following currentEvent and so
    * to the second object in unProcessedEvents.
    */
  Event * nextEvent;

  /** confidence interval used to convert double time value to long unsigned int
   *  See doubleToIntTime function for details. Note that max unsigned long int value is
   * given by constant ULONG_MAX, from limits.h.
   */
  static double tick;

  /* link to the simulation that owns this manager*/
  Simulation * simulation;

  /** add a new Event in the unProcessedEvents list
  *  \return false if Event already exists
  */
  const bool insertEvent(const std::string&, double);

  /** copy constructor => private: no copy nor pass-by-value.
   *  \param the eventsManager to be copied
   */
  EventsManager(const EventsManager&);

public:

  /** default constructor
  *  \param the simulation that owns this manager
  */
  EventsManager(Simulation* = NULL);

  /** destructor
  */
  ~EventsManager();

  /** convert time from double to unsigned int according to tick.
  */
  static const unsigned long int doubleToIntTime(double);

  /** convert time from unsigned int to double according to tick.
  */
  static const double intToDoubleTime(unsigned long int);

  /** manager initialization function
   */
  void initialize();

  /** insert time discretisation into unProcessedEvents
   * \param a pointer to the TimeDiscretisation object to schedule.
   * \param a string, the type of Event associated to this discretisation.
   */
  void scheduleTimeDiscretisation(TimeDiscretisation*, const std::string&);

  /** add a set of Events into unProcessedEvents list
   *  \param an EventsContainer
   *  \return a bool, false if insertion failed.
   */
  const bool insertEvents(const EventsContainer&);

  // GETTERS/SETTERS

  /** get the list of past Events
  *  \return a set of Events*
  */
  inline const EventsContainer getPastEvents() const
  {
    return pastEvents ;
  };

  /* No setter for member pastEvents or unProcessedEvents-> depends on simulation and can not be set in another way.
   * Moreover, only EventsManager can create new Events.
   */

  /** get the list of unProcessed Events
  *  \return a set of Events*
  */
  inline const EventsContainer getUnProcessedEvents() const
  {
    return unProcessedEvents ;
  };

  /** get the event that occurs at time inputTime
  *  \param an unsigned long int
  *  \return a pointer to Event
  */
  Event* getEventPtr(unsigned long int inputTime) const;

  /** get the current event
  *  \return a pointer to Event
  */
  inline Event* getCurrentEventPtr() const
  {
    return currentEvent;
  };

  /** get the next event to be processed.
  *  \return a pointer to Event
  */
  inline Event* getNextEventPtr() const
  {
    return nextEvent;
  };

  /** get the event following inputEvent  ("following" defined with operator(s) comparison of events)
  *  \param a pointer to Event
  *  \return a pointer to Events
  */
  Event* getFollowingEventPtr(Event*) const;

  /** get the event that follows the event at time inputTime  ("following" defined with operator(s) comparison of events)
  *  \param an unsigned long int
  *  \return a pointer to Event
  */
  Event* getFollowingEventPtr(unsigned long int inputTime) const;

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

  /** get the Simulation
  *  \return a pointer to Simulation
  */
  inline Simulation* getSimulationPtr() const
  {
    return simulation;
  }

  /** set the Simulation of the OneStepNSProblem
  *  \param: a pointer on Simulation
  */
  inline void setSimulationPtr(Simulation* str)
  {
    simulation = str;
  }

  /** check if event is present in past of unProcessedEvents list
  *  \return a bool
  */
  const bool hasEvent(Event*) const ;

  /** check if some events remain in unProcessedEvents list
  *  \return a bool
  */
  const bool hasNextEvent() const ;

  /** get the time (double format) of an event
  *  \return a double
  */
  const double getTimeOfEvent(Event*) const;

  /** get the time (int format) of an event
  *  \return an unsigned long int
  */
  const unsigned long int getIntTimeOfEvent(Event*) const;

  /** get the time of current event, in double format
  *  \return a double
  */
  const double getCurrentTime() const ;

  /** get the time of next event, in double format
  *  \return a double
  */
  const double getNextTime() const ;

  /** display EventsManager data
  */
  void display() const ;

  /** add a new Event in the unProcessedEvents list and update nextEvent value
  *  \return false if Event already exists
  */
  const bool scheduleEvent(const std::string&, double);

  /** remove an Event from the unProcessed events list
  */
  void removeEvent(Event*);

  /** update current and next event positions and run process functions of current event
  */
  void processEvents();

  /** run process functions of current event and all simultaneous events.
  */
  void process();

  /** update current and next event positions.
  */
  void shiftEvents();
};

#endif // EventsManager_H
