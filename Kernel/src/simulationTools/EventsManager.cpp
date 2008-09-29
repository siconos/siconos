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


// Creation and insertion of a new event into allEvents set.
const bool EventsManager::createAndInsertEvent(int type, double time)
{
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  return ((allEvents.insert(regEvent.instantiate(time, type)))  != allEvents.end());
}

EventsManager::EventsManager(SP::Simulation newSimu):
  currentEvent(NULL), nextEvent(NULL), ETD(NULL), ENonSmooth(NULL), simulation(newSimu), hasNS(false), hasCM(false)
{
  //  === Checks connection with a simulation ===
  if (!simulation)
    RuntimeException::selfThrow("EventsManager initialize, no simulation linked to the manager.");

  //  === Creates and inserts two events corresponding
  // to times tk and tk+1 of the simulation time-discretisation  ===
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  currentEvent = regEvent.instantiate(simulation->getTk(), 1);
  allEvents.insert(currentEvent);
  ETD = regEvent.instantiate(simulation->getTkp1(), 1);
  allEvents.insert(ETD);
  // === Set nextEvent ===
  nextEvent = ETD;
}

EventsManager::~EventsManager()
{
  nextEvent = NULL;
  allEvents.erase(currentEvent);
  allEvents.erase(ETD);
  allEvents.erase(ENonSmooth);
  delete currentEvent;
  currentEvent = NULL;
  delete ETD;
  ETD = NULL;
  delete ENonSmooth;
  ENonSmooth = NULL;
  EventsContainer::iterator it;
  for (it = allEvents.begin(); it != allEvents.end(); ++it)
  {
    if (*it != NULL) delete(*it);
  }
  allEvents.clear();
}

void EventsManager::initialize()
{
  pair<EventsContainerIterator, EventsContainerIterator> rangeNew = allEvents.equal_range(currentEvent);
  // Note: a backup is required for rangeNew since any erase/insert in loop over equal range
  // may result in invalidation of iterator and hazardous results.
  // Moreover, any "setTime" of an event need an erase/insert of the event to force resorting of the set.

  EventsContainer bckUp(rangeNew.first, rangeNew.second);

  allEvents.erase(allEvents.lower_bound(currentEvent), allEvents.upper_bound(currentEvent));
  for (EventsContainerIterator it = bckUp.begin(); it != bckUp.end() ; ++it)
  {
    (*it)->process(simulation);
    // "synchronise" actuators/sensors events with nextEvent
    if ((*it)->getType() == 3 || (*it)->getType() == 4)
      (*it)->update();

    allEvents.insert(*it);
  }
  bckUp.clear();
}

void EventsManager::insertEvents(const EventsContainer& e)
{
  allEvents.insert(e.begin(), e.end());
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
  if (!event) return false;
  EventsContainer::iterator it2 = allEvents.find(event);
  return ((it2 != allEvents.end()));
}

const bool EventsManager::hasNextEvent() const
{
  return (nextEvent);
}

const double EventsManager::getTimeOfEvent(Event* event) const
{
  //  if(!hasEvent(event))
  if (!event)
    RuntimeException::selfThrow("EventsManager getTimeOfEvent, Event == NULL (not present in the set?) ");
  return event->getDoubleTimeOfEvent();
}

const double EventsManager::getStartingTime() const
{
  if (!currentEvent)
    RuntimeException::selfThrow("EventsManager getStartingTime, current event is NULL");
  return currentEvent->getDoubleTimeOfEvent();
}

const double EventsManager::getNextTime() const
{
  if (!nextEvent)
    RuntimeException::selfThrow("EventsManager getNextTime, next event is NULL");
  return nextEvent->getDoubleTimeOfEvent();
}

void EventsManager::display() const
{
  cout << "=== EventsManager data display ===" << endl;
  if (simulation)
    cout << "- This manager belongs to the simulation named \" " << simulation->getName() << "\", of type " << simulation->getType() << "." << endl;
  else
    cout << "- No simulation linked to this manager." << endl;
  EventsContainer::iterator it;
  cout << " - The number of unprocessed events (including current one) is: " << allEvents.size() << endl;
  for (it = allEvents.begin(); it != allEvents.end(); ++it)
    (*it)->display();
  cout << "===== End of EventsManager display =====" << endl;
}

// Insertion of a new EXISTING event.
void EventsManager::insertEvent(Event* newE)
{
  if (newE->getDoubleTimeOfEvent() < currentEvent->getDoubleTimeOfEvent())
    RuntimeException::selfThrow("EventsManager insertEvent(), time is lower than current event time while it is forbidden to step back.");
  allEvents.insert(newE);
  update();
  hasCM = true;
}

void EventsManager::update()
{
  // update nextEvent value (may have change because of an insertion).
  nextEvent = getFollowingEventPtr(currentEvent);
}

// Creates (if required) and update the non smooth event of the set
// Useful during simulation when a new event is detected.
void EventsManager::scheduleNonSmoothEvent(double time)
{
  if (!ENonSmooth)
  {
    EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
    ENonSmooth = regEvent.instantiate(time, 2);
  }
  else
  {
    //allEvents.erase(ENonSmooth);
    ENonSmooth->setTime(time);
  }
  allEvents.insert(ENonSmooth);

  update();
  hasNS = true;
}


void EventsManager::removeEvent(Event* event)
{
  if (event == currentEvent)
    RuntimeException::selfThrow("EventsManager removeEvent(input), input = currentEvent, you can not remove it!");

  if (allEvents.size() < 3) // check that at least 2 events remains in the set.
    RuntimeException::selfThrow("EventsManager removeEvent(input), can not remove input, else only currentEvent will remain in the set!");

  EventsContainer::iterator it = allEvents.find(event);
  if (it != allEvents.end())
    allEvents.erase(event);
  else
    RuntimeException::selfThrow("EventsManager removeEvent(input), input is not present in the set.");

  update();
}

void EventsManager::processEvents()
{
  // No non-smooth events and no control manager
  if (!hasNS && !hasCM)
    OptimizedProcessEvents();
  else
    GeneralProcessEvents();
}

void EventsManager::OptimizedProcessEvents()
{
  // ==== Valid only when no Non Smooth event occurs and without control manager ====

  nextEvent->process(simulation);
  simulation->updateIndexSets();
  currentEvent->setTime(nextEvent->getDoubleTimeOfEvent());
  simulation->getTimeDiscretisationPtr()->increment();
  ETD->setTime(simulation->getTkp1());
  update();

  // Set Model current time
  if (nextEvent)
    simulation->getModelPtr()->setCurrentTime(getTimeOfEvent(nextEvent));

}

void EventsManager::GeneralProcessEvents()
{
  // 1 - Process all events simultaneous to nextEvent.
  //  We get a range of all the Events at time tnext and process them.
  pair<EventsContainerIterator, EventsContainerIterator> rangeNew = allEvents.equal_range(nextEvent);
  for (EventsContainerIterator it = rangeNew.first; it != rangeNew.second ; ++it)
    (*it)->process(simulation);

  // 2 - Update index sets of the simulation
  simulation->updateIndexSets();

  allEvents.erase(currentEvent);
  currentEvent->setTime(nextEvent->getDoubleTimeOfEvent());

  // Note: a backup is required for rangeNew since any erase/insert in loop over equal range
  // may result in invalidation of iterator and hazardous results.
  // Moreover, any "setTime" of an event need an erase/insert of the event to force resorting of the set.
  EventsContainer bckUp(rangeNew.first, rangeNew.second);
  allEvents.erase(allEvents.lower_bound(nextEvent), allEvents.upper_bound(nextEvent));

  // 3 - update events at time tnext.
  for (EventsContainerIterator it = bckUp.begin(); it != bckUp.end() ; ++it)
  {
    // If the event is the end of the Simulation Time Discretisation time step
    if ((*it) == ETD)
    {
      simulation->getTimeDiscretisationPtr()->increment();
      ETD->setTime(simulation->getTkp1());
      allEvents.insert(*it);

    }
    // Non Smooth Event
    else if ((*it) == ENonSmooth)
    {
      hasNS = false; // false until next insertion
    }
    // Actuator or or Sensor event
    else
    {
      (*it)->update();
      allEvents.insert(*it);
    }
  }
  bckUp.clear();
  allEvents.insert(currentEvent);
  update();

  // Set Model current time
  if (nextEvent)
    simulation->getModelPtr()->setCurrentTime(getTimeOfEvent(nextEvent));
}

void EventsManager::SortEvent(Event* e)
{
  // Temporary function to deal with add/remove events with actuators and sensors.
  allEvents.erase(e);
  allEvents.insert(e);
}
