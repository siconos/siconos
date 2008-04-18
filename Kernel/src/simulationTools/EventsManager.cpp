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
#include "EventsManager.h"
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
const bool EventsManager::insertEvent(int type, double time)
{
  EventsContainerIterator it; // to check if insertion succeed or not.
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  it = allEvents.insert(regEvent.instantiate(time, type));

  return (it != allEvents.end());
}

// PUBLIC METHODS

// Default/from data constructor
EventsManager::EventsManager(Simulation * newSimu):
  currentEvent(NULL), nextEvent(NULL), simulation(newSimu)
{}

EventsManager::~EventsManager()
{
  currentEvent = NULL;
  nextEvent = NULL;
  EventsContainer::iterator it;
  for (it = allEvents.begin(); it != allEvents.end(); ++it)
  {
    if (*it != NULL) delete(*it);
  }
  allEvents.clear();
}

void EventsManager::initialize()
{
  if (simulation == NULL)
    RuntimeException::selfThrow("EventsManager initialize, no simulation linked to the manager.");

  // === get original, user, time discretisation from Simulation object. ===
  TimeDiscretisation * td = simulation->getTimeDiscretisationPtr();

  // === insert TimeDiscretisation into the allEvents set. ===
  scheduleTimeDiscretisation(td, 1);

  // Todo: allow external events (provided by user -> xml reading or anything else) addition.
  // Mind to check that new time-events are superior to t0/currentEvent ...

  // === Set current and nextEvent ===
  currentEvent = *(allEvents.begin()); // First event in the set, normally the one at t0.
  nextEvent = getFollowingEventPtr(currentEvent);
  if (nextEvent == NULL)
    RuntimeException::selfThrow("EventsManager initialize, can not find next event since there is only one event in the set!");
}

void EventsManager::scheduleTimeDiscretisation(TimeDiscretisation* td, int type)
{
  // === get tk ===
  SiconosVector * tk = td->getTkPtr();
  unsigned int nSteps = tk->size(); // number of time steps

  // === insert tk into the allEvents set ===
  // Create unProcessed events list, by adding time discretisation events.
  bool isInsertOk;
  // isInsertOk is not really usefull: no test at the time.
  for (unsigned int i = 0 ; i < nSteps ; ++i)
    isInsertOk = insertEvent(type, (*tk)(i));
}

const bool EventsManager::insertEvents(const EventsContainer& e)
{
  EventsContainerIterator it;
  allEvents.insert(e.begin(), e.end());
  return (it != allEvents.end());
}

Event* EventsManager::getEventPtr(const mpz_t& inputTime) const
{
  EventsContainer::iterator current;
  Event * searchedEvent = NULL;
  // look for the event following the one which time is inputTime
  for (current = allEvents.begin(); current != allEvents.end(); ++current)
  {
    if (*(*current)->getTimeOfEvent() == inputTime)
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
  EventsContainer::iterator next, current = allEvents.find(inputEvent);
  // get iterator pointing to next event ...
  if (current != allEvents.end())
    next = allEvents.upper_bound(*current);
  else
    RuntimeException::selfThrow("EventsManager getFollowingEventPtr(inputEvent), Event input is not present in the set ");

  if (next == allEvents.end())
    return NULL; //RuntimeException::selfThrow("EventsManager getFollowingEventPtr(inputEvent), no next event, input is the last one in the list.");
  else
    return (*next);
}

Event* EventsManager::getFollowingEventPtr(const mpz_t& inputTime) const
{
  EventsContainer::iterator next = allEvents.upper_bound(getEventPtr(inputTime));

  if (next == allEvents.end())
    return NULL; //RuntimeException::selfThrow("EventsManager getFollowingEventPtr(inputTime), no next event, the one corresponding to inputTime is the last one in the list.");
  else
    return (*next);
}

const bool EventsManager::hasEvent(Event* event) const
{
  if (event == NULL) return false;
  EventsContainer::iterator it2 = allEvents.find(event);
  return ((it2 != allEvents.end()));
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
  return event->getDoubleTimeOfEvent();
}

const double EventsManager::getStartingTime() const
{
  if (currentEvent == NULL)
    RuntimeException::selfThrow("EventsManager getStartingTime, current event is NULL");
  return currentEvent->getDoubleTimeOfEvent();
}

const double EventsManager::getNextTime() const
{
  if (nextEvent == NULL)
    RuntimeException::selfThrow("EventsManager getNextTime, next event is NULL");
  return nextEvent->getDoubleTimeOfEvent();
}

void EventsManager::display() const
{
  cout << "=== EventsManager data display ===" << endl;
  if (simulation != NULL)
    cout << "- This manager belongs to the simulation named \" " << simulation->getName() << "\", of type " << simulation->getType() << "." << endl;
  else
    cout << "- No simulation linked to this manager." << endl;
  EventsContainer::iterator it;
  cout << " - The number of unprocessed events (including current one) is: " << allEvents.size() << endl;
  for (it = allEvents.begin(); it != allEvents.end(); ++it)
    (*it)->display();
  cout << "===== End of EventsManager display =====" << endl;
}

const bool EventsManager::scheduleEvent(int type, double time)
{
  // Check that the new time is inside Model time-bounds.
  double t0 = simulation->getModelPtr()->getT0();
  double finalT = simulation->getModelPtr()->getFinalT();
  if (time < t0 || time > finalT)
    RuntimeException::selfThrow("EventsManager scheduleEvent(..., time), time out of bounds ([t0,T]).");

  double currentTime = currentEvent->getDoubleTimeOfEvent();
  if (time < currentTime)
    RuntimeException::selfThrow("EventsManager scheduleEvent(..., time), time is lower than current event time while it is forbidden to step back.");

  // === Insert the event into the list ===
  bool isInsertOk = insertEvent(type, time);

  if (!isInsertOk)
    RuntimeException::selfThrow("EventsManager scheduleEvent(..., time): insertion of a new event failed.");

  // update nextEvent value (may have change because of insertion).
  nextEvent = getFollowingEventPtr(currentEvent);
  return isInsertOk;
}

void EventsManager::removeEvent(Event* event)
{
  if (event == currentEvent)
    RuntimeException::selfThrow("EventsManager removeEvent(input), input = currentEvent, you can not remove it!");

  if (allEvents.size() < 3) // check that at least 2 events remains in the set.
    RuntimeException::selfThrow("EventsManager removeEvent(input), can not remove input, else only currentEvent will remain in the set!");

  EventsContainer::iterator it = allEvents.find(event);
  if (it != allEvents.end())
  {
    //if ((*it)!=NULL) delete (*it);  // causes exception. Is erase enough to free memory?
    allEvents.erase(event);
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

  // Update index sets
  simulation->updateIndexSets();

  // Shift current to next ...
  shiftEvents();

  // Set Model current time
  if (nextEvent != NULL)
    simulation->getModelPtr()->setCurrentTime(getTimeOfEvent(nextEvent));
}

void EventsManager::process()
{
  // 4 - We get a range of all the Events at time tnext and process them.
  pair<EventsContainerIterator, EventsContainerIterator> rangeNew = allEvents.equal_range(nextEvent);
  EventsContainerIterator it;
  for (it = rangeNew.first; it != rangeNew.second ; ++it)
    (*it)->process(simulation);
}

void EventsManager::shiftEvents()
{
  // === update current and next event pointed values ===

  // new current event = old next event.

  EventsContainerIterator check; // to check if insertion succeed or not.
  // 1 - Get time of current event
  const mpz_t * told = currentEvent->getTimeOfEvent();

  // 2 - Get new currentEvent. currentEvent is shifted to the next event in time.
  currentEvent = getFollowingEventPtr(*told);

  // 3 - We get a range of all the Events at time told.
  pair<EventsContainerIterator, EventsContainerIterator> rangeOld = allEvents.equal_range(getEventPtr(*told));
  // Remove events that occur at time told (ie old current event) from allEvents set
  allEvents.erase(rangeOld.first, rangeOld.second);

  // Get new nextEvent (return NULL if currentEvent is the last one)
  nextEvent = getFollowingEventPtr(*currentEvent->getTimeOfEvent());
}



