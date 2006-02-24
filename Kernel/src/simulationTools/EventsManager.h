/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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

/** \class EventsManager
 *  \brief class that provides tools to handles events
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.1.
 *  \date (Creation) February 22, 2006
 *
 * Some rules:
 *     - Events can be removed or insert anywhere in the list
 *     - An event can not be modified (but removed!) => no setEvent(Event) function
 *     - All Events are created inside this class
 *     - Tick must be set at construction, and can not be change after.
 */

//#include "EventsManagerXML.h"
#ifndef EVENTSMANAGER_H
#define EVENTSMANAGER_H

#include "SiconosConst.h"
#include "RuntimeException.h"
#include "Strategy.h"
#include "TimeDiscrEvent.h"
#include "NonSmoothEvent.h"
#include "EventsComparison.h"
#include<set>
#include<deque>
#include<string>
#include<iostream>

//class EventsManagerXML;

class Strategy;

// tick default value
const double DEFAULT_TICK = 1e-09;

// The list of events, sort using operator defined in class EventsComparison.
//typedef std::set<Event*,RuntimeCmp<Event*> > eventsContainer;
/*  A stl container of type "set" is used at the time
 *  \Warning This may be not the best choice => review all possibilities */
typedef std::set<Event*, EventsComparison > eventsContainer; // sort in a chronological way
//typedef std::set<Event*,EventsComparison(reverse) > reverseEventsContainer; // reverse sorting.

// return value for insert and erase in set -> checkEventSet.second is a bool
typedef std::pair<eventsContainer::iterator, bool> checkEventSet;

class EventsManager
{
protected:

  /** list of events already processed.
   *  The last element of this list (ie the more "recent" according to EventsComparison)
   *  is currentEvent, ie the event from which simulation is supposed to start.
   */
  eventsContainer pastEvents;

  /** list of future, not processed, events.
   *  This list is not fixed and can be updated at any time
   *  depending on the simulation, user add ...
   */
  eventsContainer futureEvents;

  /** Pointer to currentEvent, ie the simulation starting point.
    * It correponds to the more recent object in pastEvents.
    */
  Event * currentEvent;

  /** Pointer to nextEvent, ie the simulation ending point.
    * It correponds to the event following currentEvent and to
    * the first object in futureEvents.
    */
  Event * nextEvent;

  /** confidence interval used to convert double time value to long unsigned int
   *  See doubleToIntTime function for details. Note that max unsigned long int value is
   * given by constant ULONG_MAX, from limits.h.
   */
  const double tick;

  /* link to the strategy that owns this manager*/
  Strategy * strategy;

  /** \fn unsigned long int doubleToIntTime(const double&)
   *  \brief convert time from double to unsigned int according to tick.
   */
  const unsigned long int doubleToIntTime(const double&) const;

  /** \fn double intToDoubleTime(const unsigned long int&) const;
   *  \brief convert time from unsigned int to double according to tick.
   */
  const double intToDoubleTime(const unsigned long int&) const;

public:

  /** \fn EventsManager(const EventsManager&)
   *  \brief copy constructor
   *  \param the eventsManager to be copied
   */
  EventsManager(const EventsManager&);

  /** \fn EventsManager(const double& = DEFAULT_TICK)
   *  \brief default constructor, with tick value as optional input
   *  \param an unsigned int
   *  \param a string
   */
  EventsManager(const double& = DEFAULT_TICK, Strategy* = NULL);

  /** \fn EventsManager(EventsManagerXML*)
   *  \brief constructor with XML object of the EventsManager
   *  \param a pointer to EventsManagerXML
   */
  // EventsManager(EventsManagerXML*);

  /** \fn ~EventsManager()
   *  \brief destructor
   */
  ~EventsManager();

  /** \fn initialize()
  *  \brief manager initialization function
  */
  void initialize();

  // GETTERS/SETTERS

  /** \fn inline const eventsContainer getPastEvents() const
   *  \brief get the list of past Events
   *  \return a set of Events*
   */
  inline const eventsContainer getPastEvents() const
  {
    return pastEvents ;
  };

  /* No setter for member pastEvents or futureEvents-> depends on simulation and can not be set in another way.
   * Moreover, only EventsManager can create new Events.
   */

  /** \fn inline const eventsContainer getFutureEvents() const
   *  \brief get the list of future Events
   *  \return a set of Events*
   */
  inline const eventsContainer getFutureEvents() const
  {
    return futureEvents ;
  };

  /** \fn Events* getNextEventPtr(Event* inputEvent) const
   *  \brief get the event following inputEvent  ("following" defined with operator(s) comparison of events)
   *  \return a pointer to Events
   */
  Event* getNextEventPtr(Event*) const;

  /** \fn inline const double getTick() const
   *  \brief get tick value
   *  \return a double
   */
  inline const double getTick() const
  {
    return tick ;
  };

  /** \fn Strategy* getStrategyPtr()
   *  \brief get the Strategy
   *  \return a pointer to Strategy
   */
  inline Strategy* getStrategyPtr() const
  {
    return strategy;
  }

  /** \fn void setStrategyPtr(Strategy*)
   *  \brief set the Strategy of the OneStepNSProblem
   *  \param: a pointer on Strategy
   */
  inline void setStrategyPtr(Strategy* str)
  {
    strategy = str;
  }

  // Events managements functions

  /** \fn const double getTimeOfEvent(Event*) const
   *  \brief get the time (double format) of an event
   *  \return a double
   */
  const double getTimeOfEvent(Event*) const;

  /** \fn const unsigned long int getIntTimeOfEvent(Event*) const
   *  \brief get the time (int format) of an event
   *  \return an unsigned long int
   */
  const unsigned long int getIntTimeOfEvent(Event*) const;

  /** \fn const double getCurrentTime() const
   *  \brief get the time of current event, in double format
   *  \return a double
   */
  const double getCurrentTime() const ;

  /** \fn const double getNextTime() const
   *  \brief get the time of next event, in double format
   *  \return a double
   */
  const double getNextTime() const ;

  /** \fn const bool hasEvent(Event* event) const
   *  \brief check if event is present in past of futureEvents list
   *  \return a bool
   */
  const bool hasEvent(Event*) const ;

  /** \fn const bool hasNextEvent() const
   *  \brief check if some events remain in futureEvents list
   *  \return a bool
   */
  const bool hasNextEvent() const ;

  /** \fn void display()
   *  \brief display EventsManager data
   */
  void display() const ;

  /** \fn const bool insertEvent(Event*)
   *  \brief add a new Event in the futureEvents list
   *  \return false if Event already exists
   */
  const bool insertEvent(const std::string&, const double&);

  /** \fn const bool scheduleEvent(Event*)
   *  \brief add a new Event in the futureEvents list and update nextEvent value
   *  \return false if Event already exists
   */
  const bool scheduleEvent(const std::string&, const double&);

  /** \fn void removeEvent(Event*)
   *  \brief remove an Event from the future events list
   */
  void removeEvent(Event*);

  /** \fn void processEvents()
   *  \brief run process functions of current event
   */
  void processEvents();
};

#endif // EventsManager_H
