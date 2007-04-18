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
#include "EventsManager.h"
#include "RuntimeCmp.h"
#include "EventFactory.h"
#include "TimeDiscretisation.h"
#include "Model.h"
#include "Simulation.h"
#include <math.h>
#include <limits> // for ULONG_MAX
using namespace std;

// insert event functions (PRIVATE).
// schedule is required by simulation.
//  -> necessary to insert a new event AND to update current/past event
// There are two different functions, to avoid multiple nextEvent update during initialize calls of insertEvent.
const bool EventsManager::insertEvent(const string& type, double time)
{
  unsigned long int intTime;

  // convert input double time to unsigned int

  intTime = doubleToIntTime(time);

  EventsContainerIterator it; // to check if insertion succeed or not.

  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  it = unProcessedEvents.insert(regEvent.instantiate(intTime, type));

  return (it != unProcessedEvents.end());
}

// PUBLIC METHODS

// RuntimeCmp object used to sort Events in the sets (pastEvents and unProcessedEvents.
// !!! \todo Find a way to avoid global variable ... !!!
RuntimeCmp<Event> compareEvents(&Event::getTimeOfEvent);

double EventsManager::tick = DEFAULT_TICK;

// Default/from data constructor
EventsManager::EventsManager(Simulation * newSimu):
  pastEvents(compareEvents), unProcessedEvents(compareEvents), currentEvent(NULL), nextEvent(NULL), simulation(newSimu)
{}

EventsManager::~EventsManager()
{
  currentEvent = NULL;
  nextEvent = NULL;
  EventsContainer::iterator it;
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

const unsigned long int EventsManager::doubleToIntTime(double doubleTime)
{
  double res = ceil(doubleTime / tick);
  if (res > ULONG_MAX) // check if res value can be converted to unsigned long int.
  {
    string report = "EventsManager doubleToIntTime, conversion results in an overflow value > ULONG_MAX, max value for unsigned long int. Try to change tick value. ULONG_MAX=";
    std::ostringstream oss;
    oss << report << ULONG_MAX;
    RuntimeException::selfThrow(oss.str());
  }
  return (unsigned long int)res;
}

const double EventsManager::intToDoubleTime(unsigned long int intTime)
{
  return tick * intTime;
}

void EventsManager::initialize()
{
  if (simulation == NULL)
    RuntimeException::selfThrow("EventsManager initialize, no simulation linked to the manager.");

  // === get original, user, time discretisation from Simulation object. ===
  TimeDiscretisation * td = simulation->getTimeDiscretisationPtr();

  // === insert TimeDiscretisation into the unProcessedEvents set. ===
  scheduleTimeDiscretisation(td, "TimeDiscretisationEvent");

  // Todo: allow external events (provided by user -> xml reading or anything else) addition.
  // Mind to check that new time-events are superior to t0/currentEvent ...

  // === Set current and nextEvent ===
  double t0 = td->getModelPtr()->getT0();
  currentEvent = getEventPtr(doubleToIntTime(t0));
  nextEvent = getFollowingEventPtr(doubleToIntTime(t0));
  if (nextEvent == NULL)
    RuntimeException::selfThrow("EventsManager initialize, can not find next event since there is only one event in the set!");
}

void EventsManager::scheduleTimeDiscretisation(TimeDiscretisation* td, const std::string& type)
{
  // === get tk ===
  SiconosVector * tk = td->getTkPtr();
  unsigned int nSteps = tk->size(); // number of time steps

  // === insert tk into the unProcessedEvents set ===
  // Create unProcessed events list, by adding time discretisation events.
  bool isInsertOk;
  // isInsertOk is not really usefull: no test at the time.
  for (unsigned int i = 0 ; i < nSteps ; ++i)
    isInsertOk = insertEvent(type, (*tk)(i));
}

const bool EventsManager::insertEvents(const EventsContainer& e)
{
  EventsContainerIterator it;
  unProcessedEvents.insert(e.begin(), e.end());
  return (it != unProcessedEvents.end());
}

Event* EventsManager::getEventPtr(unsigned long int inputTime) const
{
  EventsContainer::iterator current;
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

Event* EventsManager::getFollowingEventPtr(Event* inputEvent) const
{
  // look for inputEvent in the unProcessed events list ...
  EventsContainer::iterator next, current = unProcessedEvents.find(inputEvent);
  // get iterator pointing to next event ...
  if (current != unProcessedEvents.end())
    next = unProcessedEvents.upper_bound(*current);
  else
    RuntimeException::selfThrow("EventsManager getFollowingEventPtr(inputEvent), Event input is not present in the set ");

  if (next == unProcessedEvents.end())
    return NULL; //RuntimeException::selfThrow("EventsManager getFollowingEventPtr(inputEvent), no next event, input is the last one in the list.");
  else
    return (*next);
}

Event* EventsManager::getFollowingEventPtr(unsigned long int inputTime) const
{
  EventsContainer::iterator next = unProcessedEvents.upper_bound(getEventPtr(inputTime));

  if (next == unProcessedEvents.end())
    return NULL; //RuntimeException::selfThrow("EventsManager getFollowingEventPtr(inputTime), no next event, the one corresponding to inputTime is the last one in the list.");
  else
    return (*next);
}

const bool EventsManager::hasEvent(Event* event) const
{
  if (event == NULL) return false;
  EventsContainer::iterator it = pastEvents.find(event);
  EventsContainer::iterator it2 = unProcessedEvents.find(event);
  return ((it != pastEvents.end()) || (it2 != unProcessedEvents.end()));
}

const bool EventsManager::hasNextEvent() const
{
  return (nextEvent != NULL);
}

const double EventsManager::getTimeOfEvent(Event* event) const
{
  //  if(!hasEvent(event))
  if (event == NULL)
    RuntimeException::selfThrow("EventsManager getTimeOfEvent, Event == NULL (not present in the set?) ");
  return intToDoubleTime(event->getTimeOfEvent());
}

const unsigned long int EventsManager::getIntTimeOfEvent(Event* event) const
{
  //  if(!hasEvent(event))
  if (event == NULL)
    RuntimeException::selfThrow("EventsManager getIntTimeOfEvent,  Event == NULL (not present in the set?) ");
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
  if (simulation != NULL)
    cout << "- This manager belongs to the simulation named \" " << simulation->getName() << "\", of type " << simulation->getType() << "." << endl;
  else
    cout << "- No simulation linked to this manager." << endl;
  cout << " - Tick: " << tick << endl;
  EventsContainer::iterator it;
  cout << " - The number of already processed events is: " << pastEvents.size() << endl;
  for (it = pastEvents.begin(); it != pastEvents.end(); ++it)
    (*it)->display();
  cout << " - The number of unprocessed events (including current one) is: " << unProcessedEvents.size() << endl;
  for (it = unProcessedEvents.begin(); it != unProcessedEvents.end(); ++it)
    (*it)->display();
  cout << "===== End of EventsManager display =====" << endl;
}

const bool EventsManager::scheduleEvent(const string& type, double time)
{
  // Check that the new time is inside Model time-bounds.
  double t0 = simulation->getModelPtr()->getT0();
  double finalT = simulation->getModelPtr()->getFinalT();
  if (time < t0 || time > finalT)
    RuntimeException::selfThrow("EventsManager scheduleEvent(..., time), time out of bounds ([t0,T]).");

  // === Insert the event into the list ===
  bool isInsertOk = insertEvent(type, time);

  // Check if an event already exists at this time (=>isInsertOk == false), if so update its type
  // Usually this corresponds to a Non Smooth Event that occurs at a user time-step.
  if (!isInsertOk)
  {
    // convert time to the int format
    unsigned long int t_int = doubleToIntTime(time) ;
    // We first erase the existing element
    unProcessedEvents.erase(getEventPtr(t_int));
    // and then insert a new one at the same time but of a different type
    isInsertOk = insertEvent(type, time);
    if (currentEvent->getTimeOfEvent() == t_int)
      currentEvent = getEventPtr(t_int);
  }

  // update nextEvent value (may have change because of insertion).

  nextEvent = getFollowingEventPtr(currentEvent);
  return isInsertOk;
}

void EventsManager::removeEvent(Event* event)
{
  if (event == currentEvent)
    RuntimeException::selfThrow("EventsManager removeEvent(input), input = currentEvent, you can not remove it!");

  if (unProcessedEvents.size() < 3) // check that at least 2 events remains in the set.
    RuntimeException::selfThrow("EventsManager removeEvent(input), can not remove input, else only currentEvent will remain in the set!");

  EventsContainer::iterator it = unProcessedEvents.find(event);
  if (it != unProcessedEvents.end())
  {
    //if ((*it)!=NULL) delete (*it);  // causes exception. Is erase enough to free memory?
    unProcessedEvents.erase(event);
  }
  else
    RuntimeException::selfThrow("EventsManager removeEvent(input), input is not present in the set.");

  // update nextEvent value (may have change because of removal)
  //unsigned long int t = currentEvent->getTimeOfEvent();
  nextEvent = getFollowingEventPtr(currentEvent);
}

void EventsManager::processEvents()
{
  // Process all events simultaneous to nextEvent.
  process();

  // Shift current to next ...
  shiftEvents();

  // Set Model current time
  if (nextEvent != NULL)
    simulation->getModelPtr()->setCurrentT(getTimeOfEvent(nextEvent));
}

void EventsManager::process()
{
  // 4 - We get a range of all the Events at time tnext and process them.
  pair<EventsContainerIterator, EventsContainerIterator> rangeNew = unProcessedEvents.equal_range(nextEvent);
  EventsContainerIterator it;
  for (it = rangeNew.first; it != rangeNew.second ; ++it)
    (*it)->process(simulation);
}

void EventsManager::shiftEvents()
{
  // === update current and next event pointed values ===

  // new current event = old next event, old current move to past events.

  EventsContainerIterator check; // to check if insertion succeed or not.
  // 1 - Get time of current event
  unsigned long int told = currentEvent->getTimeOfEvent();

  // 2 - Get new currentEvent. currentEvent is shifted to the next event in time.
  currentEvent = getFollowingEventPtr(told);

  // 3 - Save Event(s) simultaneous to told into pastEvents
  // Warning: do not directly insert or remove currentEvent. Mind the pointer links!!
  // We get a range of all the Events at time told.
  pair<EventsContainerIterator, EventsContainerIterator> rangeOld = unProcessedEvents.equal_range(getEventPtr(told));
  EventsContainerIterator it;
  for (it = rangeOld.first; it != rangeOld.second ; ++it)
    pastEvents.insert(*it);

  // Remove events that occur at time told (ie old current event) from unProcessedEvents set
  unProcessedEvents.erase(rangeOld.first, rangeOld.second);

  // Get new nextEvent (return NULL if currentEvent is the last one)
  nextEvent = getFollowingEventPtr(currentEvent->getTimeOfEvent());
}


