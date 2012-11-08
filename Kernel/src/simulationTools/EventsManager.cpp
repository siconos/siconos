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
#include "EventsManager.hpp"
#include "EventFactory.hpp"
#include "TimeDiscretisation.hpp"
#include "Model.hpp"
#include "Simulation.hpp"
#include <cmath>
#include <limits> // for ULONG_MAX
using namespace std;


// Creation and insertion of a new event into allEvents set.
bool EventsManager::createAndInsertEvent(int type, double time)
{
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  return ((_allEvents.insert(regEvent.instantiate(time, type)))  != _allEvents.end());
}

EventsManager::EventsManager(): _hasNS(false), _hasCM(false)
{}


unsigned long int EventsManager::_GapLimit2Events = GAPLIMIT_DEFAULT;


void EventsManager::initialize(SP::Simulation sim)
{
  // Connection with the simulation
  assert(sim && "EventsManager::initialize(simu), input simulation object = null.");
  _simulation = sim;

  //  === Creates and inserts two events corresponding
  // to times tk and tk+1 of the simulation time-discretisation  ===
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  _currentEvent = regEvent.instantiate(simulation()->getTk(), TD_EVENT);
  _allEvents.insert(_currentEvent);
  _ETD = regEvent.instantiate(simulation()->getTkp1(), TD_EVENT);
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
    if ((*it)->getType() == SENSOR_EVENT || (*it)->getType() == ACTUATOR_EVENT)
    {
      (*it)->update();
      _allEvents.insert(*it);
    }
    // With non-smooth event
    if ((*it)->getType() == NS_EVENT)
    {
      _hasNS = false;
    }
    if ((*it) == _currentEvent)
    {
      _allEvents.insert(*it);
    }
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

bool EventsManager::hasEvent(SP::Event event) const
{
  if (!event) return false;
  EventsContainer::iterator it2 = _allEvents.find(event);
  return ((it2 != _allEvents.end()));
}

bool EventsManager::hasNextEvent() const
{
  return (bool(_nextEvent));
}

double EventsManager::getTimeOfEvent(SP::Event event) const
{
  //  if(!hasEvent(event))
  if (!event)
    RuntimeException::selfThrow("EventsManager getTimeOfEvent, Event == NULL (not present in the set?) ");
  return event->getDoubleTimeOfEvent();
}

double EventsManager::startingTime() const
{
  if (!_currentEvent)
    RuntimeException::selfThrow("EventsManager startingTime, current event is NULL");
  return _currentEvent->getDoubleTimeOfEvent();
}

double EventsManager::nextTime() const
{
  if (!_nextEvent)
    RuntimeException::selfThrow("EventsManager nextTime, next event is NULL");
  return _nextEvent->getDoubleTimeOfEvent();
}

void EventsManager::display() const
{
  cout << "=== EventsManager data display ===" << endl;
  if (simulation())
    cout << "- This manager belongs to the simulation named \"" <<
         simulation()->name() << "\", of type " <<
         Type::name(*simulation()) << "." << endl;
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
  // Check if the _ENonSmooth inserted is in the same range of the _currentEvent
  if (_hasNS) // if have one non-smooth event
  {
    const mpz_t * t1 = _currentEvent->getTimeOfEvent(); // time of the current event
    const mpz_t * t2 = _ENonSmooth->getTimeOfEvent();   // time of the non-smooth event
    if (mpz_cmp(*t1, *t2) > 0) // t1 > t2
    {
      RuntimeException::selfThrow("EventsManager::update, non-smooth event inserted must equal or greater than the current event!");
    }
    _nextEvent = _ENonSmooth; // in this case, the next event must be the non-smooth event in order to process the latter one
  }
  else
  {
    _nextEvent = followingEvent(_currentEvent);
  }
}

// Creates (if required) and update the non smooth event of the set
// Useful during simulation when a new event is detected.
void EventsManager::scheduleNonSmoothEvent(double time, bool yes_update)
{
  if (!_ENonSmooth)
  {
    EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
    _ENonSmooth = regEvent.instantiate(time, NS_EVENT);
  }
  else
  {
    //_allEvents.erase(ENonSmooth);
    _ENonSmooth->setTime(time);
  }
  _allEvents.insert(_ENonSmooth);
  _hasNS = true;
  if (yes_update) // if we need update the the next event
    update();
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
  // // For TimeStepping Scheme, need to update IndexSets, but not for EventDriven scheme,
  // if(Type::value(*simulation())== Type::TimeStepping)
  //   {
  //     simulation()->updateIndexSets();
  //   }
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
  // Note: a backup is required for rangeNew since any erase/insert in loop over equal range
  // may result in invalidation of iterator and hazardous results.
  // Moreover, any "setTime" of an event need an erase/insert of the event to force resorting of the set.
  EventsContainer bckUp(rangeNew.first, rangeNew.second);
  if (_hasNS)
  {
    // Ckeck if the _ENonSmooth is in the list of events to be processed
    EventsContainer::iterator it_nonSmooth = bckUp.find(_ENonSmooth);
    if (it_nonSmooth == bckUp.end())
    {
      RuntimeException::selfThrow("EventsManager::GeneralProcessEvents, _ENonSmooth is not in the list of events to be processed.");
    }
  }
  // Check if the current event is in the list of events to be processed
  EventsContainer::iterator it_current = bckUp.find(_currentEvent);

  if (it_current != bckUp.end()) // the _currentEvent is in the list (_currentEvent is in the same range of the _ENonSmooth)
  {
    // In this case, we suppose that all the events other than _ENonSmooth have been processed before.
    // We need process only the _ENonSmooth
    _ENonSmooth->process(simulation());
  }
  else  // _currentEvent is not in the list
  {
    // In this case, we process all the events in the list
    for (EventsContainerIterator it = rangeNew.first; it != rangeNew.second ; ++it)
    {
      (*it)->process(simulation());
    }
  }
  // 2 - Update index sets of the simulation
  simulation()->updateIndexSets();
  _allEvents.erase(_currentEvent);
  _currentEvent->setTime(_nextEvent->getDoubleTimeOfEvent());
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
    else if (((*it)->getType() == SENSOR_EVENT) || ((*it)->getType() == ACTUATOR_EVENT))
    {
      (*it)->update();
      _allEvents.insert(*it);
    }
  }
  bckUp.clear();
  //4. Check if the _currentEvent is so close to the event _ETD. If it is the case, we have to move _ETD to the end of the next time step
  const mpz_t *t1 = _currentEvent->getTimeOfEvent(); // Time of the current event
  const mpz_t *t2 = _ETD->getTimeOfEvent();          // Time of the _ETD event
  mpz_t delta_time;
  mpz_init(delta_time); // initialize delta_time
  mpz_sub(delta_time, *t2, *t1); // gap between the _currentEvent and _ETD event
  if (mpz_cmp_ui(delta_time, 0) < 0)
    RuntimeException::selfThrow("EventsManager::GeneralProcessEvents, the time of current event must be smaller than the time of ETD event.");

  if (mpz_cmp_ui(delta_time, _GapLimit2Events) <= 0) // _currentEvent is so close to _ETD
  {
    simulation()->timeDiscretisation()->increment();
    _ETD->setTime(simulation()->getTkp1());
    //    _allEvents.insert(_ETD);
  }
  // free memory
  mpz_clear(delta_time);
  //
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
