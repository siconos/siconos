/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
#include "EventsManager.h"
using namespace std;

const unsigned long int EventsManager::doubleToIntTime(const double& doubleTime) const
{
  return (unsigned long int)ceil(doubleTime / tick);
}

const double EventsManager::intToDoubleTime(const unsigned long int& intTime) const
{
  return tick * intTime;
}

// Default/from data constructor
EventsManager::EventsManager(const double& inputTick, Strategy * newStrat):
  currentEvent(NULL), nextEvent(NULL), tick(inputTick), strategy(newStrat)
{}

// copy constructor
EventsManager::EventsManager(const EventsManager& newManager):
  currentEvent(NULL), nextEvent(NULL), tick(newManager.getTick()), strategy(newManager.getStrategyPtr())
{
  // ?? allow copy of Events ?? => no.
  RuntimeException::selfThrow("EventsManager copy constructor, not yet implemented, please avoid copy!");
}

EventsManager::~EventsManager()
{
  if (currentEvent != NULL) delete currentEvent;
  currentEvent = NULL;
  nextEvent = NULL;
  eventsContainer::iterator it;
  for (it = pastEvents.begin(); it != pastEvents.end(); ++it)
  {
    if (*it != NULL) delete(*it);
  }
  pastEvents.clear();
  for (it = futureEvents.begin(); it != futureEvents.end(); ++it)
  {
    if (*it != NULL) delete(*it);
  }
  futureEvents.clear();
}

void EventsManager::initialize()
{
  if (strategy == NULL)
    RuntimeException::selfThrow("EventsManager initialize, no strategy linked to the manager.");

  // get original, user, time discretisation.
  TimeDiscretisation * td = strategy->getTimeDiscretisationPtr();
  // get tk
  SimpleVector * tk = td->getTkPtr();
  unsigned int nSteps = tk->size();

  // get initial time, initialize currentEvent and insert it in future events.
  double t0 = td->getT0();
  currentEvent = new TimeDiscrEvent(doubleToIntTime(t0));
  futureEvents.clear();
  pastEvents.clear();
  checkEventSet checkSchedule = futureEvents.insert(currentEvent);
  if (!checkSchedule.second)
    RuntimeException::selfThrow("EventsManager initialize, currentEvent insertion failed.");

  // Create future events list, by adding all time discretisation events.
  bool isInsertOk;

  for (unsigned int i = 1 ; i < nSteps ; ++i)
    isInsertOk = insertEvent("TimeDiscretisationEvent", (*tk)(i));

  // Todo: allow external events (provided by user) addition.

  // Set nextEvent
  nextEvent = getNextEventPtr(currentEvent);
}

Event* EventsManager::getNextEventPtr(Event* inputEvent) const
{
  // look for inputEvent in the future events list ...
  eventsContainer::iterator next, current = futureEvents.find(inputEvent);
  // get iterator pointing to next event ...
  if (current != futureEvents.end())
    next = futureEvents.upper_bound(*current);
  else
    RuntimeException::selfThrow("EventsManager getNextEventPtr(input), Event input is not present in the set ");

  if (next == futureEvents.end())
    RuntimeException::selfThrow("EventsManager getNextEventPtr(input), no next event, input is the last one in the list.");

  return (*next);
}

const bool EventsManager::hasEvent(Event* event) const
{
  eventsContainer::iterator it = pastEvents.find(event);
  eventsContainer::iterator it2 = futureEvents.find(event);
  return ((it != pastEvents.end()) || (it2 != futureEvents.end()));
}

const bool EventsManager::hasNextEvent() const
{
  return (futureEvents.size() < 2); // minimum size of futureEvents is one, since it contains currentEvent.
}

const double EventsManager::getTimeOfEvent(Event* event) const
{
  if (!hasEvent(event))
    RuntimeException::selfThrow("EventsManager getTimeOfEvent, Event not present in the set ");
  return intToDoubleTime(event->getTimeOfEvent());
}

const unsigned long int EventsManager::getIntTimeOfEvent(Event* event) const
{
  if (!hasEvent(event))
    RuntimeException::selfThrow("EventsManager getTimeOfEvent, Event not present in the set ");
  return event->getTimeOfEvent();
}

const double EventsManager::getCurrentTime() const
{
  if (currentEvent == NULL)
    RuntimeException::selfThrow("EventsManager getCurrentTime, current event is NULL");
  return intToDoubleTime(currentEvent->getTimeOfEvent());
}

const double EventsManager::getNextTime() const
{
  if (nextEvent == NULL)
    RuntimeException::selfThrow("EventsManager getNextTime, next event is NULL");
  return intToDoubleTime(nextEvent->getTimeOfEvent());
}

void EventsManager::display() const
{
  cout << "=== EventsManager data display ===" << endl;
  cout << " - tick: " << tick << endl;
  eventsContainer::iterator it;
  cout << " ----- Already processed events: ----- " << endl;
  for (it = pastEvents.begin(); it != pastEvents.end(); ++it)
    (*it)->display();
  cout << " ----- Future events: ----- " << endl;
  for (it = futureEvents.begin(); it != futureEvents.end(); ++it)
    (*it)->display();
  cout << "===== End of EventsManager display =====" << endl;
}

// insert event functions.
// schedule is required by strategy.
//  -> necessary to insert a new event AND to update current/past event
// There are two different functions, to avoid multiple nextEvent update during initialize calls of insertEvent.
const bool EventsManager::insertEvent(const string& type, const double& time)
{
  unsigned long int intTime;

  // convert input double time to unsigned int
  intTime = doubleToIntTime(time);

  checkEventSet checkSchedule; // to check if insertion succeed or not.

  // scan various possible types for events
  if (type == "TimeDiscretisationEvent")
    checkSchedule = futureEvents.insert(new TimeDiscrEvent(intTime));
  else if (type == "NonSmoothEvent")
    checkSchedule = futureEvents.insert(new NonSmoothEvent(intTime));
  else
    RuntimeException::selfThrow("EventsManager scheduleEvent, unknown event type" + type);
  return checkSchedule.second;  // true if insert succeed.
}

const bool EventsManager::scheduleEvent(const string& type, const double& time)
{
  // Insert the event into the list
  bool isInsertOk = insertEvent(type, time);
  // update nextEvent value.(may have change because of insertion)
  nextEvent = getNextEventPtr(currentEvent);
  return isInsertOk;  // true if insert succeed.
}

void EventsManager::removeEvent(Event* event)
{
  eventsContainer::iterator it = futureEvents.find(event);
  if (it != futureEvents.end())
  {
    if ((*it) != NULL) delete(*it);
    futureEvents.erase(event);
  }
  else
    RuntimeException::selfThrow("EventsManager removeEvent, Event not present in the set ");

  // update nextEvent value (may have change because of removal)
  nextEvent = getNextEventPtr(currentEvent);
}

void EventsManager::processEvents()
{
  checkEventSet check; // to check if insertion succeed or not.

  // insert current event into pastEvents and remove it from futureEvents
  check = pastEvents.insert(currentEvent);
  if (!check.second) RuntimeException::selfThrow("EventsManager, processEvents, event insertion in pastEvents failed.");

  futureEvents.erase(currentEvent);

  // set new current event
  currentEvent = nextEvent;

  // get new next event
  nextEvent = getNextEventPtr(currentEvent);

  // process events
  eventsContainer::iterator it;
  for (it = futureEvents.begin(); it != futureEvents.end(); ++it)
    (*it)->process();
}

//void EventsManager::saveEventsManagerToXML()
//{
//  RuntimeException::selfThrow("saveEventsManagerToXML: not yet implemented");
//}
