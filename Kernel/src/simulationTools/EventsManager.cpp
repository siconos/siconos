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
#include "EventsManager.hpp"
#include "EventFactory.hpp"
#include "TimeDiscretisation.hpp"
#include "Model.hpp"
#include "Simulation.hpp"
#include <math.h>
#include <limits> // for ULONG_MAX
using namespace std;


// Creation and insertion of a new event into allEvents set.
const bool EventsManager::createAndInsertEvent(int type, double time)
{
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  return ((_allEvents.insert(regEvent.instantiate(time, type)))  != _allEvents.end());
}

EventsManager::EventsManager(): _hasNS(false), _hasCM(false)
{}

void EventsManager::initialize(SP::Simulation sim)
{
  // Connection with the simulation
  assert(sim && "EventsManager::initialize(simu), input simulation object = null.");
  _simulation = sim;

  //  === Creates and inserts two events corresponding
  // to times tk and tk+1 of the simulation time-discretisation  ===
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  _currentEvent = regEvent.instantiate(simulation()->getTk(), 1);
  _allEvents.insert(_currentEvent);
  _ETD = regEvent.instantiate(simulation()->getTkp1(), 1);
  _allEvents.insert(_ETD);
  // === Set nextEvent ===
  _nextEvent = _ETD;
}

void EventsManager::preUpdate()
{
  pair<EventsContainerIterator, EventsContainerIterator> rangeNew = _allEvents.equal_range(_currentEvent);
  // Note: a backup is required for rangeNew since any erase/insert in loop over equal range
  // may result in invalidation of iterator and hazardous results.
  // Moreover, any "setTime" of an event need an erase/insert of the event to force resorting of the set.

  EventsContainer bckUp(rangeNew.first, rangeNew.second);

  _allEvents.erase(_allEvents.lower_bound(_currentEvent), _allEvents.upper_bound(_currentEvent));

  for (EventsContainerIterator it = bckUp.begin(); it != bckUp.end() ; ++it)
  {
    (*it)->process(simulation());
    // "synchronise" actuators/sensors events with nextEvent
    if ((*it)->getType() == 3 || (*it)->getType() == 4)
      (*it)->update();

    _allEvents.insert(*it);
  }
  bckUp.clear();
}

void EventsManager::insertEvents(const EventsContainer& e)
{
  _allEvents.insert(e.begin(), e.end());
}

SP::Event EventsManager::event(const mpz_t& inputTime) const
{
  EventsContainer::iterator current;
  SP::Event searchedEvent;
  // look for the event following the one which time is inputTime
  for (current = _allEvents.begin(); current != _allEvents.end(); ++current)
  {
    if (*(*current)->getTimeOfEvent() == inputTime)
    {
      searchedEvent = *current;
      break;
    }
  }
  return searchedEvent;
}

SP::Event EventsManager::followingEvent(SP::Event inputEvent) const
{
  // look for inputEvent in the unProcessed events list ...
  EventsContainer::iterator next, current = _allEvents.find(inputEvent);
  // get iterator pointing to next event ...
  if (current != _allEvents.end())
    next = _allEvents.upper_bound(*current);
  else
    RuntimeException::selfThrow("EventsManager followingEvent(inputEvent), Event input is not present in the set ");

  if (next == _allEvents.end())
    return SP::Event(); //RuntimeException::selfThrow("EventsManager followingEvent(inputEvent), no next event, input is the last one in the list.");
  else
    return (*next);
}

SP::Event EventsManager::followingEvent(const mpz_t& inputTime) const
{
  EventsContainer::iterator next = _allEvents.upper_bound(event(inputTime));

  if (next == _allEvents.end())
    return SP::Event(); //RuntimeException::selfThrow("EventsManager followingEvent(inputTime), no next event, the one corresponding to inputTime is the last one in the list.");
  else
    return (*next);
}

const bool EventsManager::hasEvent(SP::Event event) const
{
  if (!event) return false;
  EventsContainer::iterator it2 = _allEvents.find(event);
  return ((it2 != _allEvents.end()));
}

const bool EventsManager::hasNextEvent() const
{
  return (_nextEvent);
}

const double EventsManager::getTimeOfEvent(SP::Event event) const
{
  //  if(!hasEvent(event))
  if (!event)
    RuntimeException::selfThrow("EventsManager getTimeOfEvent, Event == NULL (not present in the set?) ");
  return event->getDoubleTimeOfEvent();
}

const double EventsManager::startingTime() const
{
  if (!_currentEvent)
    RuntimeException::selfThrow("EventsManager startingTime, current event is NULL");
  return _currentEvent->getDoubleTimeOfEvent();
}

const double EventsManager::nextTime() const
{
  if (!_nextEvent)
    RuntimeException::selfThrow("EventsManager nextTime, next event is NULL");
  return _nextEvent->getDoubleTimeOfEvent();
}

void EventsManager::display() const
{
  cout << "=== EventsManager data display ===" << endl;
  if (simulation())
    cout << "- This manager belongs to the simulation named \" " << simulation()->name() << "\", of type " << simulation()->getType() << "." << endl;
  else
    cout << "- No simulation linked to this manager." << endl;
  EventsContainer::iterator it;
  cout << " - The number of unprocessed events (including current one) is: " << _allEvents.size() << endl;
  for (it = _allEvents.begin(); it != _allEvents.end(); ++it)
    (*it)->display();
  cout << "===== End of EventsManager display =====" << endl;
}

// Insertion of a new EXISTING event.
void EventsManager::insertEvent(SP::Event newE)
{
  if (newE->getDoubleTimeOfEvent() < _currentEvent->getDoubleTimeOfEvent())
    RuntimeException::selfThrow("EventsManager insertEvent(), time is lower than current event time while it is forbidden to step back.");
  _allEvents.insert(newE);
  update();
  _hasCM = true;
}

void EventsManager::update()
{
  // update nextEvent value (may have change because of an insertion).
  _nextEvent = followingEvent(_currentEvent);
}

// Creates (if required) and update the non smooth event of the set
// Useful during simulation when a new event is detected.
void EventsManager::scheduleNonSmoothEvent(double time)
{
  if (!_ENonSmooth)
  {
    EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
    _ENonSmooth = regEvent.instantiate(time, 2);
  }
  else
  {
    //_allEvents.erase(ENonSmooth);
    _ENonSmooth->setTime(time);
  }
  _allEvents.insert(_ENonSmooth);

  update();
  _hasNS = true;
}


void EventsManager::removeEvent(SP::Event event)
{
  if (event == _currentEvent)
    RuntimeException::selfThrow("EventsManager removeEvent(input), input = _currentEvent, you can not remove it!");

  if (_allEvents.size() < 3) // check that at least 2 events remains in the set.
    RuntimeException::selfThrow("EventsManager removeEvent(input), can not remove input, else only _currentEvent will remain in the set!");

  EventsContainer::iterator it = _allEvents.find(event);
  if (it != _allEvents.end())
    _allEvents.erase(event);
  else
    RuntimeException::selfThrow("EventsManager removeEvent(input), input is not present in the set.");

  update();
}

void EventsManager::processEvents()
{
  // No non-smooth events and no control manager
  if (!_hasNS && !_hasCM)
    OptimizedProcessEvents();
  else
    GeneralProcessEvents();
}

void EventsManager::OptimizedProcessEvents()
{
  // ==== Valid only when no Non Smooth event occurs and without control manager ====

  _nextEvent->process(simulation());
  simulation()->updateIndexSets();
  _currentEvent->setTime(_nextEvent->getDoubleTimeOfEvent());
  simulation()->timeDiscretisation()->increment();
  _ETD->setTime(simulation()->getTkp1());
  update();

  // Set Model current time
  if (_nextEvent)
    simulation()->model()->setCurrentTime(getTimeOfEvent(_nextEvent));

}

void EventsManager::GeneralProcessEvents()
{
  // 1 - Process all events simultaneous to nextEvent.
  //  We get a range of all the Events at time tnext and process them.
  pair<EventsContainerIterator, EventsContainerIterator> rangeNew = _allEvents.equal_range(_nextEvent);
  for (EventsContainerIterator it = rangeNew.first; it != rangeNew.second ; ++it)
    (*it)->process(simulation());

  // 2 - Update index sets of the simulation
  simulation()->updateIndexSets();

  _allEvents.erase(_currentEvent);
  _currentEvent->setTime(_nextEvent->getDoubleTimeOfEvent());

  // Note: a backup is required for rangeNew since any erase/insert in loop over equal range
  // may result in invalidation of iterator and hazardous results.
  // Moreover, any "setTime" of an event need an erase/insert of the event to force resorting of the set.
  EventsContainer bckUp(rangeNew.first, rangeNew.second);
  _allEvents.erase(_allEvents.lower_bound(_nextEvent), _allEvents.upper_bound(_nextEvent));

  // 3 - update events at time tnext.
  for (EventsContainerIterator it = bckUp.begin(); it != bckUp.end() ; ++it)
  {
    // If the event is the end of the Simulation Time Discretisation time step
    if ((*it) == _ETD)
    {
      simulation()->timeDiscretisation()->increment();
      _ETD->setTime(simulation()->getTkp1());
      _allEvents.insert(*it);

    }
    // Non Smooth Event
    else if ((*it) == _ENonSmooth)
    {
      _hasNS = false; // false until next insertion
    }
    // Actuator or or Sensor event
    else
    {
      (*it)->update();
      _allEvents.insert(*it);
    }
  }
  bckUp.clear();
  _allEvents.insert(_currentEvent);
  update();

  // Set Model current time
  if (_nextEvent)
    simulation()->model()->setCurrentTime(getTimeOfEvent(_nextEvent));
}

void EventsManager::SortEvent(SP::Event e)
{
  // Temporary function to deal with add/remove events with actuators and sensors.
  _allEvents.erase(e);
  _allEvents.insert(e);
}
