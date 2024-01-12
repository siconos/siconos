/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "Event.hpp"
#include "SiconosConst.hpp"  // siconos::internal::GAPLIMIT_DEFAULT
#include "SiconosException.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeDiscretisationEvent.hpp"

unsigned long int siconos::simulation::EventsManager::_GapLimit2Events =
    siconos::internal::GAPLIMIT_DEFAULT;

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::simulation::EventsManager::EventsManager(std::shared_ptr<TimeDiscretisation> td)
    : _td(td) {
  //  === Creates and inserts two events corresponding
  // to times tk and tk+1 of the simulation time-discretisation  ===
  _events.push_back(EventFactory::instance()->create(_td->getTk(0), EventType::TD));
  _events.push_back(EventFactory::instance()->create(_td->getTk(1), EventType::TD));
  _events.push_back(EventFactory::instance()->create(_td->getTk(2), EventType::TD));
  _events[1]->setTimeDiscretisation(_td);
  _events[2]->setTimeDiscretisation(_td);
  _events[1]->setK(_k + 1);
  _events[2]->setK(_k + 2);
}

void siconos::simulation::EventsManager::initialize(double T) { _T = T; }

// Creation and insertion of a new event into the event set.
siconos::simulation::Event& siconos::simulation::EventsManager::insertEvent(
    siconos::simulation::EventType type, double time) {
  DEBUG_BEGIN(
      "Event& siconos::simulation::EventsManager::insertEvent(int type, double time)\n");
  // Uses the events factory to insert the new event.
  auto pos = insertEv(EventFactory::instance()->create(time, type));
  DEBUG_END("Event& siconos::simulation::EventsManager::insertEvent(int type, double time)\n");
  return *_events[pos];
}

siconos::simulation::Event& siconos::simulation::EventsManager::insertEvent(
    siconos::simulation::EventType type, std::shared_ptr<TimeDiscretisation> td) {
  auto& ev = insertEvent(type, td->getTk(_k));
  ev.setTimeDiscretisation(td);
  return ev;
}

double siconos::simulation::EventsManager::getTk() { return _td->getTk(_k); }

double siconos::simulation::EventsManager::getTkp1() const {
  double tkp1 = _td->getTk(_k + 1);
  if (tkp1 <= _T + 100.0 * std::numeric_limits<double>::epsilon())
    return tkp1;
  else
    return std::numeric_limits<double>::quiet_NaN();
};

double siconos::simulation::EventsManager::getTkp2() const {
  double tkp2 = _td->getTk(_k + 2);
  if (tkp2 <= _T + 100.0 * std::numeric_limits<double>::epsilon())
    return tkp2;
  else
    return std::numeric_limits<double>::quiet_NaN();
};

double siconos::simulation::EventsManager::getTkp3() const {
  double tkp3 = _td->getTk(_k + 3);
  if (tkp3 <= _T + 100.0 * std::numeric_limits<double>::epsilon())
    return tkp3;
  else
    return std::numeric_limits<double>::quiet_NaN();
};

void siconos::simulation::EventsManager::noSaveInMemory(const Simulation& sim) {
  for (auto& it : _events) {
    if (it->getType() == EventType::TD) {
      static_pointer_cast<TimeDiscretisationEvent>(it)->noSaveInMemory();
    }
  }
  //   Event& ev = *it;
  //   if (ev.getType() == EventType::TD) {
  //     it =
  //     std::make_shared<TimeDiscretisationEventNoSaveInMemory>(ev.getDoubleTimeOfEvent(),
  //                                                                  0);
  //     it->setTimeDiscretisation(ev.timeDiscretisation());
  //   }
  // }
}

void siconos::simulation::EventsManager::preUpdate(Simulation& sim) {
  // Note FP: seems to be unused/obsolete. To be reviewed.
  DEBUG_BEGIN("siconos::simulation::EventsManager::preUpdate(Simulation& sim)\n");
  DEBUG_EXPR(display(););
  const mpz_t* t1 = _events[0]->getTimeOfEvent();
  _events[0]->process(sim);
  for (unsigned int i = 1; i < _events.size(); i++) {
    const mpz_t* t2 = _events[i]->getTimeOfEvent();
    int res = mpz_cmp(*t1, *t2);
    if (res == 0) {
      if (_events[i]->getType() == EventType::NS) {
        _events[i]->process(sim);
        _events.erase(_events.begin() + i);
      }
    } else
      break;
  }
  DEBUG_END("siconos::simulation::EventsManager::preUpdate(Simulation& sim)\n");
}

double siconos::simulation::EventsManager::startingTime() const {
  if (_events.size() == 0)
    THROW_EXCEPTION(
        "siconos::simulation::EventsManager::startingTime current event is nullptr");
  return _events[0]->getDoubleTimeOfEvent();
}

double siconos::simulation::EventsManager::nextTime() const {
  if (_events.size() <= 1) THROW_EXCEPTION("EventsManager nextTime, next event is nullptr");
  return _events[1]->getDoubleTimeOfEvent();
}

bool siconos::simulation::EventsManager::needsIntegration() const {
  if (_events.size() <= 1) THROW_EXCEPTION("EventsManager nextTime, next event is nullptr");
  return (mpz_cmp(*_events[0]->getTimeOfEvent(), *_events[1]->getTimeOfEvent()) < 0);
}

// Creates (if required) and update the non smooth event of the set
// Useful during simulation when a new event is detected.
void siconos::simulation::EventsManager::scheduleNonSmoothEvent(Simulation& sim, double time,
                                                                bool yes_update) {
  if (!_eNonSmooth) {
    _eNonSmooth = EventFactory::instance()->create(time, EventType::NS);
  } else {
    _eNonSmooth->setTime(time);
  }

  // NonsmoothEvent is special, we need to take care of it.
  // If a NS event is scheduled too close to a TD event, LsodarOSI will refuse to
  // integrate from the NS event to the TD event. Thus we just delete the TD event.
  // In fact we just skip a t_k in this case
  //
  // First thing to do is to look for the next TD event
  const mpz_t* t1 = _eNonSmooth->getTimeOfEvent();
  auto pos = insertEv(_eNonSmooth);
  // looking for a TD event close to the NS one.
  mpz_t delta_time = {};
  mpz_init(delta_time);  // initialize delta_time
  for (unsigned int j = 1; j < _events.size(); j++) {
    if (j == pos) continue;
    Event& ev = *_events[j];
    if (ev.getType() != EventType::TD)  // current event is not of type TD
      continue;
    mpz_sub(delta_time, *ev.getTimeOfEvent(), *t1);  // gap between the NS and TD events
    if (mpz_cmp_ui(delta_time, 0) < 0)               // ok
      continue;
    if (mpz_cmp_ui(delta_time, _GapLimit2Events) <= 0)  // the two are too close
    {
      // reschedule the TD event only if its time instant is less than T
      if (!std::isnan(getTkp3())) {
        _NSeventInsteadOfTD = true;
        static_cast<TimeDiscretisationEvent&>(ev).update(_k + 3);
        insertEv(_events[j]);
      }
      // delete the TD event (that has to be done in all cases)
      _events.erase(_events.begin() + j);
      break;
    }
  }
  mpz_clear(delta_time);
}

void siconos::simulation::EventsManager::processEvents(Simulation& sim) {
  // process next event
  _events[1]->process(sim);

  // update the event stack
  update(sim);
}

void siconos::simulation::EventsManager::update(Simulation& sim) {
  // delete last event, since we have processed one
  auto event0Type = _events[0]->getType();
  // reschedule a TD event if needed
  if (event0Type == EventType::TD && _k != 0) {
    // this checks whether the next time instant is less than T or not
    // if it isn't then tkp1 is a NaN, in which case we don't reschedule the event
    // and the simulation will stop
    // TODO: create a TD at T if T ∈ (t_k, t_{k+1}), so the simulation effectively
    // run until T
    double tkp2 = getTkp2();
    std::static_pointer_cast<TimeDiscretisationEvent>(_events[0])->update(_k + 2);
    if (!std::isnan(tkp2)) {
      insertEv(_events[0]);
    }
  }
  // reschedule if needed
  else if (_events[0]->reschedule()) {
    _events[0]->update();
    if (_events[0]->getDoubleTimeOfEvent() <
        _T + 100.0 * std::numeric_limits<double>::epsilon())
      insertEv(_events[0]);
  }
  // A NS_EVENT was schedule close to a TD_EVENT
  // the latter was removed, but we still need to increase
  // the current index
  else if (event0Type == EventType::NS && _NSeventInsteadOfTD) {
    _NSeventInsteadOfTD = false;
    _k++;
  }

  // unconditionally remove the previous processed event
  _events.erase(_events.begin());

  // Now we can update _k if we have processed a TD_EVENT
  if (_events[0]->getType() == EventType::TD) _k++;
}

unsigned int siconos::simulation::EventsManager::insertEv(std::shared_ptr<Event> new_event) {
  mpz_t* t1 = const_cast<mpz_t*>(new_event->getTimeOfEvent());
  // const auto eType = new_event->getType();
  bool inserted = false;
  unsigned int pos = 0;
  mpz_t delta_time;
  mpz_init(delta_time);  // initialize delta_time
  mpz_t abs_delta_time;
  mpz_init(abs_delta_time);  // initialize delta_time

  // Find a place for the event in the vector
  for (auto it = _events.begin(); it != _events.end(); ++it) {
    auto ev = *it;
    // // delta = t(ev) - t(new_event)
    mpz_sub(delta_time, *(ev->getTimeOfEvent()), *t1);
    auto compare = mpz_cmp_ui(delta_time, _GapLimit2Events);

    // if (mpz_cmp_ui(delta_time, _GapLimit2Events) > 0)
    if (compare > 0)  // current event time > new event time > ==> insert just before
    {
      _events.insert(it, new_event);
      inserted = true;
      break;
    } else {
      mpz_abs(abs_delta_time, delta_time);
      //Let us check if the new time is not too close to the time of the current event
      if (mpz_cmp_ui(abs_delta_time, _GapLimit2Events) <= 0)  // the two are too close
      {
        // reschedule new_event at the same time as the one of the current event
        // and only if its type is different from the type of the current event
        mpz_set(*t1, *(ev->getTimeOfEvent()));
      }
      //   // Warning FP: we use the type (enum) of the event
      //   // to decide if it must be inserted before or after the current event.
      //   // So the priority is defined in Event.hpp, with the enum EventType.
      //   if (eType != ev->getType()) {
      //     if (ev->getType() == EventType::TD) {
      // 	    _events.insert(std::next(it), new_event);
      // 	  pos++;
      //     //} else
      //     //  _events.insert(it, new_event);
      //     inserted = true;
      //     break;
      //   }
      // }
    }
    pos++;
  }

  if (!inserted) _events.push_back(new_event);

  mpz_clear(delta_time);
  mpz_clear(abs_delta_time);
  return pos;
}

double siconos::simulation::EventsManager::currentTimeStep() { return _td->timeStep(_k); }

void siconos::simulation::EventsManager::display() const {
  std::cout << "=== EventsManager data display ===\n";
  std::cout << " - The number of unprocessed events (including current one) is: "
            << _events.size() << "\n";
  for (auto it : _events) it->display();
  std::cout << "===== End of EventsManager display =====\n";
}

bool siconos::simulation::EventsManager::compareEvents::operator()(std::shared_ptr<Event> e1,
                                                                   std::shared_ptr<Event> e2) {
  mpz_t* t1 = const_cast<mpz_t*>(e1->getTimeOfEvent());
  mpz_t d_time;
  mpz_init(d_time);  // initialize delta_time
  mpz_sub(d_time, *(e2->getTimeOfEvent()), *t1);
  auto compare = mpz_cmp_ui(d_time, _GapLimit2Events);
  e1->display();
  e2->display();
  std::cout << "Compare ?" << compare << "\n";
  if (compare > 0)  // e1 > e2
    return false;
  else {
    // mpz_t abs_delta_time;
    // mpz_init(abs_delta_time);  // initialize delta_time

    // mpz_abs(abs_delta_time, d_time);
    // // Let us check if the new time is not too close to the time of the current event
    // if (mpz_cmp_ui(abs_delta_time, _GapLimit2Events) <= 0)  // the two are too close
    // {  // e1.time == e2.time, let's compare types
    //   auto e1Type = e1->getType();
    //   auto e2Type = e2->getType();
    //   if (e1Type <= e2Type)
    //     return true;
    //   else if (e1Type > e2Type)
    //     return false;
    // }
    // else
    return true;
  }
};
