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
#include "TimeDiscretisationEvent.hpp"
#include "Model.hpp"
#include "Simulation.hpp"
#include <cmath>
#include <limits> // for ULONG_MAX
#include "TimeDiscretisationEventNoSaveInMemory.hpp"
#include "CxxStd.hpp"

unsigned long int EventsManager::_GapLimit2Events = GAPLIMIT_DEFAULT;

EventsManager::EventsManager(SP::TimeDiscretisation td): _k(0), _td(td),
   _NSeventInsteadOfTD(false)
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
Event& EventsManager::insertEvent(const int type, const double& time)
{
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get());
  unsigned int pos = insertEv(regEvent.instantiate(time, type));
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
  for (EventsContainer::iterator it = _events.begin();
       it != _events.end(); ++it)
  {
    Event& ev = **it;
    if (ev.getType() == TD_EVENT)
    {
      (*it).reset(new TimeDiscretisationEventNoSaveInMemory(ev.getDoubleTimeOfEvent(), 0));
      (*it)->setTimeDiscretisation(ev.getTimeDiscretisation());
    }
  }
}

void EventsManager::preUpdate(Simulation& sim)
{
  const mpz_t *t1 = _events[0]->getTimeOfEvent();
  for (EventsContainer::iterator it = _events.begin(); it != _events.end() ; ++it)
  {
    const  mpz_t *t2 = (*it)->getTimeOfEvent();
    int res = mpz_cmp(*t1, *t2);
    if (res == 0)
    {
      if ((*it)->getType() != SENSOR_EVENT && (*it)->getType() != ACTUATOR_EVENT && (*it)->getType() != OBSERVER_EVENT)
      {
        (*it)->process(sim);
      // "synchronise" actuators/sensors events
      // XXX needed ???
//      if ((*it)->getType() == SENSOR_EVENT || (*it)->getType() == ACTUATOR_EVENT || (*it)->getType() == OBSERVER_EVENT)
//      {
//        (*it)->update();
//        insertEv(*it);
//        _events.erase(it);
//      }
        if (it == _events.begin())
          continue;
        else _events.erase(it);
      }
    }
  }
}

double EventsManager::startingTime() const
{
  if (_events.size() == 0)
    RuntimeException::selfThrow("EventsManager::startingTime current event is NULL");
  return _events[0]->getDoubleTimeOfEvent();
}

double EventsManager::nextTime() const
{
  if (_events.size() <= 1)
    RuntimeException::selfThrow("EventsManager nextTime, next event is NULL");
  return _events[1]->getDoubleTimeOfEvent();
}


// Creates (if required) and update the non smooth event of the set
// Useful during simulation when a new event is detected.
void EventsManager::scheduleNonSmoothEvent(Simulation& sim, double time, bool yes_update)
{
  if (!_eNonSmooth)
  {
    EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
    _eNonSmooth = regEvent.instantiate(time, NS_EVENT);
  }
  else
  {
    _eNonSmooth->setTime(time);
  }

  // NonsmoothEvent is special, we need to take care of it.
  // If a NS event is scheduled too close to a TD event, Lsodar will refuse to
  // integrate from the NS event to the TD event. Thus we just delete the TD event.
  // In fact we just skip a t_k in this case
  //
  // First thing to do is to look for the next TD event
  const mpz_t *t1 = _eNonSmooth->getTimeOfEvent();
  unsigned int pos;
  pos = insertEv(_eNonSmooth);
  // looking for a TD event close to the NS one.
  mpz_t delta_time;
  mpz_init(delta_time); // initialize delta_time
  for (unsigned int j = 1; j < _events.size(); j++)
  {
    if (j == pos)
      continue;
    Event& ev = *_events[j];
    if (ev.getType() != TD_EVENT) // current event is not of type TD
      continue;
    mpz_sub(delta_time, *ev.getTimeOfEvent(), *t1); // gap between the NS and TD events
    if (mpz_cmp_ui(delta_time, 0) < 0) // ok
      continue;
    if (mpz_cmp_ui(delta_time, _GapLimit2Events) <= 0) // the two are too close
    {
      // reschedule the TD event only if its time instant is less than T
      if (!isnan(getTkp3()))
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
  if (event0Type == TD_EVENT)
  {
    // this checks whether the next time instant is less than T or not
    // it is isn't then tkp1 is a NaN, in which case we don't reschedule the event
    // and the simulation will stop
    // TODO: create a TD at T if T ∈ (t_k, t_{k+1}), so the simulation effectively
    // run until T
    double tkp2 = getTkp2();
    std11::static_pointer_cast<TimeDiscretisationEvent>(_events[0])->update(_k+2);
    if (!isnan(tkp2))
    {
      insertEv(_events[0]);
    }
  }
  // reschedule if needed
  else if (_events[0]->reschedule())
  {
    _events[0]->update();
    if (_events[0]->getDoubleTimeOfEvent() < _T + 100.0*std::numeric_limits<double>::epsilon())
      insertEv(_events[0]);
  }
  // An NS_EVENT was schedule close to a TD_EVENT
  // the latter was removed, but we still need to increase
  // the current index
  else if (event0Type == NS_EVENT && _NSeventInsteadOfTD)
  {
    _NSeventInsteadOfTD = false;
    _k++;
  }

  // unconditionally remove previous processed event
  _events.erase(_events.begin());

  // Now we may update _k if we have processed a TD_EVENT
  if (_events[0]->getType() == TD_EVENT)
    _k++;
}

unsigned int EventsManager::insertEv(SP::Event e)
{
  const mpz_t *t1 = e->getTimeOfEvent();
  const unsigned int eType = e->getType();
  int res;
  bool inserted = false;
  unsigned int pos = 0;
  // Find a place for the event in the vector
  for (EventsContainer::iterator it = _events.begin();
       it != _events.end(); ++it)
  {
    Event& ev = **it;
    const mpz_t *t2 = ev.getTimeOfEvent();
    res = mpz_cmp(*t1, *t2);
    if (res == 0)
      res = eType - ev.getType();
    if (res < 0)
    {
      _events.insert(it, e);
      inserted = true;
      break;
    }
    pos++;
  }

  if (!inserted)
    _events.push_back(e);

  return pos;
}

void EventsManager::display() const
{
  std::cout << "=== EventsManager data display ===" <<std::endl;
  std::cout << " - The number of unprocessed events (including current one) is: " << _events.size() <<std::endl;
  for (EventsContainer::const_iterator it = _events.begin(); it != _events.end(); ++it)
    (*it)->display();
  std::cout << "===== End of EventsManager display =====" <<std::endl;
}
