/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include <limits> // for ULONG_MAX
using namespace std;

// PRIVATE METHODS

const unsigned long int EventsManager::doubleToIntTime(const double& doubleTime) const
{
  double res = ceil(doubleTime / tick);
  if (res > ULONG_MAX) // check if res value can be converted to unsigned long int.
    RuntimeException::selfThrow("EventsManager doubleToIntTime, conversion results in an overflow value > ULONG_MAX, max value for unsigned long int. Try to change tick value. ULONG_MAX=" + ULONG_MAX);
  return (unsigned long int)res;
}

const double EventsManager::intToDoubleTime(const unsigned long int& intTime) const
{
  return tick * intTime;
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
    checkSchedule = unProcessedEvents.insert(new TimeDiscrEvent(intTime));
  else if (type == "NonSmoothEvent")
    checkSchedule = unProcessedEvents.insert(new NonSmoothEvent(intTime));
  else
    RuntimeException::selfThrow("EventsManager insertEvent, unknown event type" + type);
  return checkSchedule.second;  // true if insert succeed.
}

// PUBLIC METHODS

// RuntimeCmp object used to sort Events in the sets (pastEvents and unProcessedEvents.
// !!! \todo Find a way to avoid global variable ... !!!
RuntimeCmp<Event> compareEvents(&Event::getTimeOfEvent);

// Default/from data constructor
EventsManager::EventsManager(const double& inputTick, Strategy * newStrat):
  pastEvents(compareEvents), unProcessedEvents(compareEvents), currentEvent(NULL), nextEvent(NULL), tick(inputTick), strategy(newStrat)
{}

// copy constructor
EventsManager::EventsManager(const EventsManager& newManager):
  pastEvents(compareEvents), unProcessedEvents(compareEvents), currentEvent(NULL), nextEvent(NULL), tick(newManager.getTick()), strategy(newManager.getStrategyPtr())
{
  // ?? allow copy of Events ?? => no.
  RuntimeException::selfThrow("EventsManager copy constructor, not yet implemented, please avoid copy!");
}

EventsManager::~EventsManager()
{
  currentEvent = NULL;
  nextEvent = NULL;
  eventsContainer::iterator it;
  for (it = pastEvents.begin(); it != pastEvents.end(); ++it)
  {
    if (*it != NULL) delete(*it);
  }
  pastEvents.clear();
  for (it = unProcessedEvents.begin(); it != unProcessedEvents.end(); ++it)
  {
    if (*it != NULL) delete(*it);
  }
  unProcessedEvents.clear();
}

void EventsManager::initialize()
{
  if (strategy == NULL)
    RuntimeException::selfThrow("EventsManager initialize, no strategy linked to the manager.");

  // === get original, user, time discretisation. ===
  TimeDiscretisation * td = strategy->getTimeDiscretisationPtr();

  // === insert TimeDiscretisation into the unProcessedEvents set. ===
  scheduleTimeDiscretisation(td);

  // Todo: allow external events (provided by user -> xml reading or anything else) addition.
  // Mind to check that new time-events are superior to t0/currentEvent ...

  // === Set current and nextEvent ===
  double t0 = td->getT0();
  currentEvent = getEventPtr(doubleToIntTime(t0));
  if (hasNextEvent())
    nextEvent = getNextEventPtr(doubleToIntTime(t0));
  else
    RuntimeException::selfThrow("EventsManager initialize, can not find next event since there is only one event in the set!");
}

void EventsManager::scheduleTimeDiscretisation(TimeDiscretisation* td)
{
  // === Clear unProcessedEvents set ===
  eventsContainer::iterator it;
  for (it = unProcessedEvents.begin(); it != unProcessedEvents.end(); ++it)
  {
    if (*it != NULL) delete(*it);
  }
  unProcessedEvents.clear();

  // === get tk ===
  SimpleVector * tk = td->getTkPtr();
  unsigned int nSteps = tk->size(); // number of time steps

  // === insert tk into the unProcessedEvents set ===
  // Create unProcessed events list, by adding time discretisation events.
  bool isInsertOk;
  // isInsertOk is not really usefull: no need to test, since unProcessedEvents is cleared at the beginning of this function
  for (unsigned int i = 0 ; i < nSteps ; ++i)
    isInsertOk = insertEvent("TimeDiscretisationEvent", (*tk)(i));
}

Event* EventsManager::getEventPtr(const unsigned long int& inputTime) const
{
  eventsContainer::iterator current;
  Event * searchedEvent = NULL;
  // look for the event following the one which time is inputTime
  for (current = unProcessedEvents.begin(); current != unProcessedEvents.end(); ++current)
  {
    if ((*current)->getTimeOfEvent() == inputTime)
    {
      searchedEvent = *current;
      break;
    }
  }
  if (searchedEvent == NULL)
    RuntimeException::selfThrow("EventsManager getEventPtr(inputTime), no Event corresponding to that time in the set.");

  return searchedEvent;
}

Event* EventsManager::getNextEventPtr(Event* inputEvent) const
{
  // look for inputEvent in the unProcessed events list ...
  eventsContainer::iterator next, current = unProcessedEvents.find(inputEvent);
  // get iterator pointing to next event ...
  if (current != unProcessedEvents.end())
    next = unProcessedEvents.upper_bound(*current);
  else
    RuntimeException::selfThrow("EventsManager getNextEventPtr(inputEvent), Event input is not present in the set ");

  if (next == unProcessedEvents.end())
    return NULL; //RuntimeException::selfThrow("EventsManager getNextEventPtr(inputEvent), no next event, input is the last one in the list.");
  else
    return (*next);
}

Event* EventsManager::getNextEventPtr(const unsigned long int& inputTime) const
{
  eventsContainer::iterator next = unProcessedEvents.upper_bound(getEventPtr(inputTime));

  if (next == unProcessedEvents.end())
    return NULL; //RuntimeException::selfThrow("EventsManager getNextEventPtr(inputTime), no next event, the one corresponding to inputTime is the last one in the list.");
  else
    return (*next);
}

const bool EventsManager::hasEvent(Event* event) const
{
  if (event == NULL) return false;
  eventsContainer::iterator it = pastEvents.find(event);
  eventsContainer::iterator it2 = unProcessedEvents.find(event);
  return ((it != pastEvents.end()) || (it2 != unProcessedEvents.end()));
}

const bool EventsManager::hasNextEvent() const
{
  return (!(unProcessedEvents.size() < 2)); // minimum size of unProcessedEvents is one, since it contains currentEvent.
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
  if (strategy != NULL)
    cout << "- This manager belongs to the strategy named \" " << strategy->getName() << "\", of type " << strategy->getType() << "." << endl;
  else
    cout << "- No strategy linked to this manager." << endl;
  cout << " - Tick: " << tick << endl;
  eventsContainer::iterator it;
  cout << " - The number of already processed events is: " << pastEvents.size() << endl;
  for (it = pastEvents.begin(); it != pastEvents.end(); ++it)
    (*it)->display();
  cout << " - The number of unprocessed events (including current one) is: " << unProcessedEvents.size() << endl;
  for (it = unProcessedEvents.begin(); it != unProcessedEvents.end(); ++it)
    (*it)->display();
  cout << "===== End of EventsManager display =====" << endl;
}

const bool EventsManager::scheduleEvent(const string& type, const double& time)
{
  // === get original, user, time discretisation. ===
  //     and check that t0 < time < T
  TimeDiscretisation * td = strategy->getTimeDiscretisationPtr();
  if (time < td->getT0() || time > td->getT())
    RuntimeException::selfThrow("EventsManager scheduleEvent(..., time), time out of bounds ([t0,T]).");

  // === Insert the event into the list ===
  bool isInsertOk = insertEvent(type, time);
  // update nextEvent value.(may have change because of insertion)
  nextEvent = getNextEventPtr(currentEvent);
  return isInsertOk;  // true if insert succeed.
}

void EventsManager::removeEvent(Event* event)
{
  if (event == currentEvent)
    RuntimeException::selfThrow("EventsManager removeEvent(input), input = currentEvent, you can not remove it!");

  if (unProcessedEvents.size() < 3) // check that at least 2 events remains in the set.
    RuntimeException::selfThrow("EventsManager removeEvent(input), can not remove input, else only currentEvent will remain in the set!");

  eventsContainer::iterator it = unProcessedEvents.find(event);
  if (it != unProcessedEvents.end())
  {
    //if ((*it)!=NULL) delete (*it);  // causes exception. Is erase enough to free memory?
    unProcessedEvents.erase(event);
  }
  else
    RuntimeException::selfThrow("EventsManager removeEvent(input), input is not present in the set.");

  // update nextEvent value (may have change because of removal)
  //unsigned long int t = currentEvent->getTimeOfEvent();
  nextEvent = getNextEventPtr(currentEvent);
}

void EventsManager::processEvents()
{
  // === process events ===

  eventsContainer::iterator it;
  for (it = unProcessedEvents.begin(); it != unProcessedEvents.end(); ++it)
    (*it)->process();

  // === update current and next event pointed values ===

  checkEventSet check; // to check if insertion succeed or not.
  // Get time of current event
  unsigned long int t = currentEvent->getTimeOfEvent();
  // Save currentEvent into pastEvents
  // Warning: do not directly insert or remove currentEvent. Mind the pointer links!!
  check = pastEvents.insert(getEventPtr(t));
  if (!check.second) RuntimeException::selfThrow("EventsManager, processEvents, event insertion in pastEvents failed.");

  // Get new currentEvent
  currentEvent = getNextEventPtr(t);

  // Remove event that occurs at time t (ie old current event) from unProcessedEvents set
  unProcessedEvents.erase(getEventPtr(t));

  // Get new nextEvent (return NULL if currentEvent is the last one)
  nextEvent = getNextEventPtr(currentEvent->getTimeOfEvent());
}

//void EventsManager::saveEventsManagerToXML()
//{
//  RuntimeException::selfThrow("saveEventsManagerToXML: not yet implemented");
//}
