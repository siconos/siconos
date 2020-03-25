/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "EventsManager.hpp"
#include "EventFactory.hpp"
#include "TimeDiscretisationEvent.hpp"
#include "TimeDiscretisationEventNoSaveInMemory.hpp"
#include "Simulation.hpp"
#include <cmath>
#include <limits> // for ULONG_MAX
#include "CxxStd.hpp"
#include <gmp.h>
#include <iostream>
#include <set>

unsigned long int EventsManager::_GapLimit2Events = GAPLIMIT_DEFAULT;

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

EventsManager::EventsManager(SP::TimeDiscretisation td): _k(0), _td(td),
  _T(std::numeric_limits<double>::infinity()), _NSeventInsteadOfTD(false)
{
  //  === Creates and inserts two events corresponding
  // to times tk and tk+1 of the simulation time-discretisation  ===
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  _events.push_back(regEvent.instantiate(_td->getTk(0), TD_EVENT));
  _events[0]->setType(-1); // this is just a dumb event
  double tkp1 = _td->getTk(1);
  double tkp2 = _td->getTk(2);
  _events.push_back(regEvent.instantiate(tkp1, TD_EVENT));
  _events.push_back(regEvent.instantiate(tkp2, TD_EVENT));
  _events[1]->setTimeDiscretisation(_td);
  _events[2]->setTimeDiscretisation(_td);
  _events[1]->setK(_k+1);
  _events[2]->setK(_k+2);
}

void EventsManager::initialize(double T)
{
  _T = T;
}

// Creation and insertion of a new event into the event set.
Event& EventsManager::insertEvent(int type, double time)
{
  DEBUG_BEGIN("Event& EventsManager::insertEvent(int type, double time)\n");
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get());
  unsigned int pos = insertEv(regEvent.instantiate(time, type));
  DEBUG_END("Event& EventsManager::insertEvent(int type, double time)\n");
  return *_events[pos];
}

Event& EventsManager::insertEvent(const int type, SP::TimeDiscretisation td)
{
  Event& ev = insertEvent(type, td->getTk(_k));
  ev.setTimeDiscretisation(td);
  return ev;
}

void EventsManager::noSaveInMemory(const Simulation& sim)
{
  for(EventsContainer::iterator it = _events.begin();
      it != _events.end(); ++it)
  {
    Event& ev = **it;
    if(ev.getType() == TD_EVENT)
    {
      (*it).reset(new TimeDiscretisationEventNoSaveInMemory(ev.getDoubleTimeOfEvent(), 0));
      (*it)->setTimeDiscretisation(ev.getTimeDiscretisation());
    }
  }
}

void EventsManager::preUpdate(Simulation& sim)
{
  DEBUG_BEGIN("EventsManager::preUpdate(Simulation& sim)\n");
  DEBUG_EXPR(display(););
  const mpz_t *t1 = _events[0]->getTimeOfEvent();
  _events[0]->process(sim);
  for(unsigned int i = 1; i < _events.size() ; i++)
  {
    const  mpz_t *t2 =  _events[i]->getTimeOfEvent();
    int res = mpz_cmp(*t1, *t2);
    if(res == 0)
    {
      if(_events[i]->getType() == NS_EVENT)
      {
        _events[i]->process(sim);
        _events.erase(_events.begin()+i);
      }
    }
    else
      break;
  }
  DEBUG_END("EventsManager::preUpdate(Simulation& sim)\n");
}

double EventsManager::startingTime() const
{
  if(_events.size() == 0)
    RuntimeException::selfThrow("EventsManager::startingTime current event is nullptr");
  return _events[0]->getDoubleTimeOfEvent();
}

double EventsManager::nextTime() const
{
  if(_events.size() <= 1)
    RuntimeException::selfThrow("EventsManager nextTime, next event is nullptr");
  return _events[1]->getDoubleTimeOfEvent();
}

bool EventsManager::needsIntegration() const
{
  if(_events.size() <= 1)
    RuntimeException::selfThrow("EventsManager nextTime, next event is nullptr");
  return (mpz_cmp(*_events[0]->getTimeOfEvent(), *_events[1]->getTimeOfEvent()) < 0);
}

// Creates (if required) and update the non smooth event of the set
// Useful during simulation when a new event is detected.
void EventsManager::scheduleNonSmoothEvent(Simulation& sim, double time, bool yes_update)
{
  if(!_eNonSmooth)
  {
    EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
    _eNonSmooth = regEvent.instantiate(time, NS_EVENT);
  }
  else
  {
    _eNonSmooth->setTime(time);
  }

  // NonsmoothEvent is special, we need to take care of it.
  // If a NS event is scheduled too close to a TD event, LsodarOSI will refuse to
  // integrate from the NS event to the TD event. Thus we just delete the TD event.
  // In fact we just skip a t_k in this case
  //
  // First thing to do is to look for the next TD event
  const mpz_t *t1 = _eNonSmooth->getTimeOfEvent();
  unsigned int pos;
  pos = insertEv(_eNonSmooth);
  // looking for a TD event close to the NS one.
  mpz_t delta_time = {};
  mpz_init(delta_time); // initialize delta_time
  for(unsigned int j = 1; j < _events.size(); j++)
  {
    if(j == pos)
      continue;
    Event& ev = *_events[j];
    if(ev.getType() != TD_EVENT)  // current event is not of type TD
      continue;
    mpz_sub(delta_time, *ev.getTimeOfEvent(), *t1); // gap between the NS and TD events
    if(mpz_cmp_ui(delta_time, 0) < 0)  // ok
      continue;
    if(mpz_cmp_ui(delta_time, _GapLimit2Events) <= 0)  // the two are too close
    {
      // reschedule the TD event only if its time instant is less than T
      if(!isnan(getTkp3()))
      {
        _NSeventInsteadOfTD = true;
        static_cast<TimeDiscretisationEvent&>(ev).update(_k+3);
        insertEv(_events[j]);
      }
      // delete the TD event (that has to be done in all cases)
      _events.erase(_events.begin()+j);
      break;
    }
  }
  mpz_clear(delta_time);
}

void EventsManager::processEvents(Simulation& sim)
{
  //process next event
  _events[1]->process(sim);

  // update the event stack
  update(sim);
}

void EventsManager::update(Simulation& sim)
{
  // delete last event, since we have processed one
  int event0Type = _events[0]->getType();
  // reschedule an TD event if needed
  if(event0Type == TD_EVENT)
  {
    // this checks whether the next time instant is less than T or not
    // it is isn't then tkp1 is a NaN, in which case we don't reschedule the event
    // and the simulation will stop
    // TODO: create a TD at T if T ∈ (t_k, t_{k+1}), so the simulation effectively
    // run until T
    double tkp2 = getTkp2();
    std11::static_pointer_cast<TimeDiscretisationEvent>(_events[0])->update(_k+2);
    if(!isnan(tkp2))
    {
      insertEv(_events[0]);
    }
  }
  // reschedule if needed
  else if(_events[0]->reschedule())
  {
    _events[0]->update();
    if(_events[0]->getDoubleTimeOfEvent() < _T + 100.0*std::numeric_limits<double>::epsilon())
      insertEv(_events[0]);
  }
  // An NS_EVENT was schedule close to a TD_EVENT
  // the latter was removed, but we still need to increase
  // the current index
  else if(event0Type == NS_EVENT && _NSeventInsteadOfTD)
  {
    _NSeventInsteadOfTD = false;
    _k++;
  }

  // unconditionally remove previous processed event
  _events.erase(_events.begin());

  // Now we may update _k if we have processed a TD_EVENT
  if(_events[0]->getType() == TD_EVENT)
    _k++;
}

unsigned int EventsManager::insertEv(SP::Event e)
{
  mpz_t *t1 = const_cast<mpz_t*>(e->getTimeOfEvent());
  const unsigned int eType = e->getType();
  bool inserted = false;
  unsigned int pos = 0;
  mpz_t delta_time;
  mpz_init(delta_time); // initialize delta_time
  mpz_t abs_delta_time;
  mpz_init(abs_delta_time); // initialize delta_time
  // Find a place for the event in the vector
  for(EventsContainer::iterator it = _events.begin();
      it != _events.end(); ++it)
  {
    Event& ev = **it;
    mpz_sub(delta_time, *ev.getTimeOfEvent(), *t1); // delta = t_existing_event - t_event_to _insert
    int res = mpz_cmp_ui(delta_time, _GapLimit2Events);
    if(res > 0)  // insert
    {
      _events.insert(it, e);
      inserted = true;
      break;
    }
    else
    {
      mpz_abs(abs_delta_time, delta_time);
      if(mpz_cmp_ui(abs_delta_time, _GapLimit2Events) <= 0)  // the two are too close
      {
        // reschedule the TD event only if its time instant is less than T
        mpz_set(*t1, *ev.getTimeOfEvent());
        res = eType - ev.getType();
        if(res < 0)
        {
          _events.insert(it, e);
          inserted = true;
          break;
        }
      }
    }
    pos++;
  }

  if(!inserted)
    _events.push_back(e);

  mpz_clear(delta_time);
  mpz_clear(abs_delta_time);
  return pos;
}

void EventsManager::display() const
{
  std::cout << "=== EventsManager data display ===" <<std::endl;
  std::cout << " - The number of unprocessed events (including current one) is: " << _events.size() <<std::endl;
  for(EventsContainer::const_iterator it = _events.begin(); it != _events.end(); ++it)
    (*it)->display();
  std::cout << "===== End of EventsManager display =====" <<std::endl;
}
