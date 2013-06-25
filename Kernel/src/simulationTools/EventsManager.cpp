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


unsigned long int EventsManager::_GapLimit2Events = GAPLIMIT_DEFAULT;

EventsManager::EventsManager()
{}

void EventsManager::initialize(const Simulation& sim)
{
  //  === Creates and inserts two events corresponding
  // to times tk and tk+1 of the simulation time-discretisation  ===
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  _events.push_back(regEvent.instantiate(sim.getTk(), TD_EVENT));
  _events[0]->setType(-1); // this is just a dumb event
#if __cplusplus >= 201103L
  assert(!::isnan(sim.getTkp1()));
#else
  assert(!std::isnan(sim.getTkp1()));
#endif
#if __cplusplus >= 201103L
  assert(!::isnan(sim.getTkp2()));
#else
  assert(!std::isnan(sim.getTkp2()));
#endif
  _events.push_back(regEvent.instantiate(sim.getTkp1(), TD_EVENT));
  _events.push_back(regEvent.instantiate(sim.getTkp2(), TD_EVENT));
}

// Creation and insertion of a new event into the event set.
Event& EventsManager::insertEvent(const int type, const double& time)
{
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get());
  unsigned int pos = insertEv(regEvent.instantiate(time, type));
  return *_events[pos];
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
      (*it)->process(sim);
      // "synchronise" actuators/sensors events
      // XXX needed ???
      if ((*it)->getType() == SENSOR_EVENT || (*it)->getType() == ACTUATOR_EVENT || (*it)->getType() == OBSERVER_EVENT)
      {
        (*it)->update();
        insertEv(*it);
        _events.erase(it);
      }
      else if (it == _events.begin())
        continue;
      else _events.erase(it);
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
      sim.timeDiscretisation()->increment();
      // reschedule the TD event only if its time instant is less than T
#if __cplusplus >= 201103L
      if (!::isnan(sim.getTkp1()))
#else
      if (!std::isnan(sim.getTkp1()))
#endif
      {
        ev.setTime(sim.getTkp1());
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

  // If last processed event is a TD event, increment TD in the simulation
  // We have a problem at the start of the simulation, since the Event at t=t0
  // is just here to fill the first stop
  if (_events[0]->getType() == TD_EVENT)
    sim.timeDiscretisation()->increment();

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
    // TODO: create a TD at T if T ∈ (t_k, t_{k+1}), so the simulation effectivly
    // run until T
    double tkp2 = sim.getTkp2();
#if __cplusplus >= 201103L
    if (!::isnan(tkp2))
#else
    if (!std::isnan(tkp2))
#endif
    {
      _events[0]->setTime(tkp2);
      insertEv(_events[0]);
    }
  }
  // reschedule a Actuator or Sensor event if needed
  else if ((event0Type == SENSOR_EVENT) || (event0Type == ACTUATOR_EVENT) || (event0Type == OBSERVER_EVENT))
  {
    _events[0]->update();
    if (_events[0]->getDoubleTimeOfEvent() < sim.model()->finalT() + 100.0*std::numeric_limits<double>::epsilon())
      insertEv(_events[0]);
  }

  // remove previous processed event
  _events.erase(_events.begin());
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
