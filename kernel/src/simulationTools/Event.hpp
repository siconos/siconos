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
/*! \file Event.hpp
  General interface for Events
*/

#ifndef Event_H
#define Event_H

// #include <cmath>
#include <gmp.h>

#include <cassert>
#include <cmath>  // rint
#include <functional>
#include <map>
#include <memory>

#include "SiconosSerialization.hpp"

// As always, MSVC miss C99
// #if defined(_MSC_VER) && _MSC_VER < 1800
// extern "C" double rint(double x);
// #endif

namespace siconos::simulation {

class TimeDiscretisation;
class Simulation;

/** To set the 'id' (type) of events */
enum class EventType {
  // Warning FP: when two events have the same time, their type (EventType)
  // is used to decide which one will be the first in the events vector.
  // The smallest the type, the earlier the event.
  /** Time discretisation event */
  TD = 1,
  /** Nonsmooth event */
  NS = 2,
  /** Sensor (control toolbox) */
  Sensor,
  /** Observer (control toolbox) */
  Observer,
  /** Actuator (control toolbox) */
  Actuator,
  /** User defined: extra ids to allow users to define their own events */
  UserDefined1,
  UserDefined2
};

/**
  Abstract class that represents generic time events.

  This base class simply records the time at which the event will take place. A pure virtual
 function named process will be invoked to execute the event. The time is represented with a
 mpz_t, from gmp library. See http://gmplib.org.

  Derived classes:
  - TimeDiscretisationEvent: events that corresponds to user-defined time-discretisation
 points
  - NonSmoothEvent: specific events, detected during simulation, when constraints are
 violated (thanks to roots-finding algorithm)
  - SensorEvent: event dedicated to data capture through user-defined sensors.

*/

class Event {
 protected:
  ACCEPT_SERIALIZATION(Event);

  /** Date of the present event,
       represented with a mpz_t */
  mpz_t _timeOfEvent;

  /** Number of ticks corresponding to a timestep */
  mpz_t _tickIncrement;

  /** Id or type of the Event */
  EventType _type{EventType::TD};

  /** Date of the present event,
   *  represented with a double */
  double _dTime{0.};

  /** confidence interval used to convert double time value to mpz_t */
  static double _tick;

  /** True if at least one Event object has already been instanciated.
      Used to detect potentially dangerous cases in setTick
  */
  static bool _eventCreated;

  /** index for the current Event*/
  unsigned int _k{0};

  /** TimeDiscretisation for the Event (unused only in the NonSmoothEvent) */
  std::shared_ptr<TimeDiscretisation> _td{nullptr};

  /** For automatic rescheduling */
  bool _reschedule{false};

  // /** Default constructor */
  // Event()
  // {
  //   // what is the default values of timeOfEvent and tickIncrement (mpz_t default) ?
  //   mpz_init(_timeOfEvent);
  //   mpz_init(_tickIncrement);
  // };

  // Rule of five
  Event() = delete;
  Event(const Event&) = delete;
  Event(Event&&) = delete;
  Event& operator=(const Event&) = delete;
  Event& operator=(Event&&) = delete;

 public:
  /** constructor with time value and type as input
   *
   *  \param time the starting type (a double)
   *  \param newType the Event type (an int)
   *  \param reschedule set this to true if the event has to be rescheduled
   */
  Event(double time, EventType newType, bool reschedule = false);

  /** destructor
   */
  virtual ~Event() noexcept;

  /** get tick value
   *
   *  \return a double
   */
  inline double getTick() const { return _tick; };

  /** set tick value
   *
   *  \param newTick the new tick value
   */
  static void setTick(double newTick);

  /** get the time of the present event (mpz_t format)
   *
   *  \return a mpz_t
   */
  inline const mpz_t* getTimeOfEvent() const { return &_timeOfEvent; };

  /** get the time of the present event (double format)
   *
   *  \return a double
   */
  inline double getDoubleTimeOfEvent() const { return _dTime; }

  inline void incrementTime(unsigned int step = 1) {
    for (unsigned int i = 0; i < step; i++)
      mpz_add(_timeOfEvent, _timeOfEvent, _tickIncrement);
    _dTime = mpz_get_d(_timeOfEvent) * _tick;
  }

  /** set the time of the present event (double format)
   *
   *  \param time the new time
   */
  inline void setTime(double time) {
    _dTime = time;
    mpz_set_d(_timeOfEvent, rint(_dTime / _tick));
  };

  /** \return the type (enum, either TD or NS) of the event
   */
  inline auto getType() const { return _type; };

  /** set the type (enum, either TD or NS) of the event
   *
   *  \param newType the new Event type
   */
  inline void setType(EventType newType) { _type = newType; };

  /** Set the current step k
   *
   *  \param newK the new value of _k
   */
  inline void setK(unsigned int newK) { _k = newK; };

  /** Set the TimeDiscretisation
   *
   *  \param td a TimeDiscretisation for this Event
   */
  void setTimeDiscretisation(std::shared_ptr<TimeDiscretisation> td);

  /** Get the TimeDiscretisation
   *
   *  \return the TimeDiscretisation used in this Event
   */
  inline std::shared_ptr<TimeDiscretisation> timeDiscretisation() const { return _td; };

  /** display Event data
   */
  void display() const;

  /** virtual function which actions depends on event type
   *
   *  \param sim the simulation that owns this Event (through the EventsManager)
   */
  virtual void process(Simulation& sim) = 0;

  /** virtual function which actions depends on event type.
   *  The generic implementation present in this object is to increment the
   *  TimeDiscretisation and to chamge the time of the current Event
   *
   *  \param k meaning depends on the type of event. See derived class.
   */
  virtual void update(unsigned int k = 0);

  inline bool reschedule() const { return _reschedule; };
};

/** A class to handle events creation

  Requirements:
  - the Event type must be known and register

  For a XXXXEvent, add in the file describing the XXXXEvent class:

  static siconos::simulation::EventRegistration<siconos::simulation::XXXEvent>
 reg(siconos::simulation::EventType::XXXX);

  See EventType enum for the available names.

  Usage:

  auto event = EventFactory::instance()->create(time, EventType::XXXX)

*/
class EventFactory {
  // Signature of Event constructor
  using EventCreator = std::function<std::shared_ptr<Event>(double)>;

  /** map to connect event type and the function used to create them */
  std::map<EventType, EventCreator> m_factories;

 public:
  /** Factory function which creates and returns an event

      \param time instant time of the event
      \param type type of the event (must be a EventType (enum))
      \return a pointer to event
  */
  std::shared_ptr<Event> create(double time, EventType type) {
    assert(m_factories.contains(type) && "unknown Event type");
    return m_factories[type](time);
  }

  /** access to the (singleton) factory instance */
  static EventFactory* instance() {
    static EventFactory factory;
    return &factory;
  }

  void registerCreator(EventType newtype, EventCreator caller) {
    m_factories[newtype] = caller;
  }
};

template <class T>
class EventRegistration {
 public:
  EventRegistration(EventType newtype) {
    EventFactory::instance()->registerCreator(newtype,
                                              [](double a) { return std::make_shared<T>(a); });
  }
};

}  // namespace siconos::simulation
#endif  // Event_H
