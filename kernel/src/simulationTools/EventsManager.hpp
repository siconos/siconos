/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
/*! \file EventsManager.hpp
  \brief Management of a Set (STL) of Events for the Simulation process
*/

#ifndef EventsManager_H
#define EventsManager_H

#include <vector>
#include <limits>
#include "SiconosFwd.hpp"
#include "TimeDiscretisation.hpp"

const unsigned long int GAPLIMIT_DEFAULT = 100;

/** set of events, with an ordering based on Event time value (mpz_t) to compare Events */
typedef std::vector<SP::Event> EventsContainer; // Event are already sorted

/** Tools to handle a set of Events for the Simulation

   \author SICONOS Development Team - copyright INRIA
   \version 3.6.0.
   \date (Creation) December 22, 2012

   The EventsManager handles a set of events (from user time-discretisation, sensors, non-smooth ...),
   and is supposed to provide to the simulation the values of "current" and "next" events to define
   the time-integration interval.

   Events:
   - currentEvent: starting time for integration. Initialized with t0 of the simulation time-discretisation.
   - ETD: corresponds to the instant t[k+1] of the Simulation (user) TimeDiscretisation.
   - ENonSmooth: for EventDriven simulation only. Present only if one or more non-smooth events have been
   detected between currentEvent and the next event.
   - Sensor or Actuators Events. To each Sensor or Actuator declared in the ControlManager, corresponds
   an Event in the manager. When this event is processed, its time value is increased to next instant
   in the time-discretisation of the sensor/actuator.

   Examples:
   - for a TimeStepping, with one Sensor, the EventsManager looks like:\n
   {currentEvent(tk), ESensor(tsensor), ETD(tk+1)}
   - for an EventDriven, with one actuator, an non-smooth event detected at time tns: \n
   {currentEvent(tk), EActuator(tact), ENonSmooth(tns), ETD(tk+1)}.

   After each process, the time values of each event are updated and nextEvent points to the first event after currentEvent.

   \section EMMfunc Main functions
   - initialize(): process all events which have the same time as currentEvent
   - processEvents(): process all events simultaneous to nextEvent, increment them to next step, update index sets,
   increment currentEvent.


*/
class EventsManager
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(EventsManager);

  /** list of events
   * The first element is the last processed event. All the others
   * are unprocessed events.
   * This list is not fixed and can be updated at any time
   * depending on the simulation, user add ...
   */
  EventsContainer _events;

  /** Placeholder for the non smooth event */
  SP::Event _eNonSmooth;

  /** Current index (for time instants) */
  unsigned int _k;

  /** TimeDiscretisation for the time integration*/
  SP::TimeDiscretisation _td;

  /** End of the Simulation */
  double _T;

  /** unsigned long int variable to check if two events are too close */
  static unsigned long int _GapLimit2Events;

  /** boolean to remember that a TD_EVENT has been deleted since a
   * NS_EVENT was too close */
  bool _NSeventInsteadOfTD;

  /** Insert an event in the event stack
   * \param e the event to insert
   * \return the position of the inserted event in the stack
   */
  unsigned int insertEv(SP::Event e);

  /** Update the set of events
   * \param sim the Simulation using this EventsManager
   */
  void update(Simulation& sim);

  /** default constructor */
  EventsManager() {};

public:

  /**  default constructor
   * \param td the TimeDiscretisation used in the Simulation
   */
  EventsManager(SP::TimeDiscretisation td);

  /** destructor
   */
  virtual ~EventsManager() {};

  /** Initialize: just set the final time
   * \param T the final time of the Simulation
   */
  void initialize(double T);

  /** Set the gap limit between two events
   * \param var the new _GapLimit2Events
   */
  inline void setGapLimitEvents(unsigned long int var)
  {
    _GapLimit2Events = var;
  };

  /** Get the gap limit between two events
   * \return the gap limit
   */
  inline unsigned long int getGapLimitEvents() const
  {
    return _GapLimit2Events;
  };

  /** Change TimeDiscretisationEvent to TimeDiscretisationEventNoSaveInMemory
   * \warning use this at your own risk, many integrators needs previous values
   * to integrate properly
   * \param sim the Simulation that owns this EventsManager
   */
  void noSaveInMemory(const Simulation& sim);

  /** get the current event
   *  \return a pointer to Event
   */
  inline SP::Event currentEvent() const
  {
    return _events[0];
  };

  /** get the next event to be processed.
   *  \return a pointer to Event
   */
  inline SP::Event nextEvent() const
  {
    return _events[1];
  };

  /** return all the events
   * \return a reference to the events set
   */
  inline EventsContainer& events()
  {
    return _events;
  };

  /** check if there are some unprocessed events
   *  \return true if there are unprocessed events
   */
  inline bool hasNextEvent() const
  {
    return _events.size() > 1;
  };

  /** get the time of current event, in double format
   *  \return the time of the last processed events
   */
  double startingTime() const ;

  /** get the time of next event, in double format
   *  \return the time of the next events
   */
  double nextTime() const ;

  /** is an integration step required ? The current event and the next one may have
   * the same time instant in which case no integration as to be performed
   * \return true if the simulation needs to be integrate, no otherwise
   */
  bool needsIntegration() const;

  /** display EventsManager data
   */
  void display() const ;

  /** add a new Event in the allEvents list and update nextEvent value
   * \param sim the simulation that owns this EventsManager
   * \param time the time (double format) of occurence of the event
   * \param yes_update indicator to update or not the next event (default value is true)
   */
  void scheduleNonSmoothEvent(Simulation& sim, double time, bool yes_update = true);

  /** Process the next event, update the indexSets if necessary
   * \param sim the simulation that owns this EventsManager
   */
  void processEvents(Simulation& sim);

  /** Function to be called once after initialization.
   * It is used to process NonSmoothEvents at the beginning of
   * the Simulation, if there is any.
   * \param sim the simulation that owns this EventsManager
   */
  void preUpdate(Simulation& sim);

  /** insert an event of a certain type. The event is created on the fly.
   * \param type the type of the event
   * \param time the time of the event
   * \return a reference to the Event
   */
  Event& insertEvent(int type, double time);

  /** insert an event of a certain type. The event is created on the fly,
   * and the SP::TimeDiscretisation given in argument is stored inside
   * \param type the type of the event
   * \param td a TimeDiscretisation for the Event
   * \return a reference to the Event
   */
  Event& insertEvent(int type, SP::TimeDiscretisation td);

  double getTk()
  {
    return _td->getTk(_k);
  }

  /** get time instant k+1 of the time discretisation
   * \return a double. If the simulation is near the end (t_{k+1} >= T), it returns NaN.
   */
  inline double getTkp1() const
  {
    double tkp1 = _td->getTk(_k+1);
    if (tkp1 <= _T + 100.0*std::numeric_limits<double>::epsilon())
      return tkp1;
    else
      return std::numeric_limits<double>::quiet_NaN();
  };

  /** get time instant k+2 of the time discretisation
      \return a double. If the simulation is near the end (t_{k+2} >= T), it returns NaN.
  */
  inline double getTkp2() const
  {
    double tkp2 = _td->getTk(_k+2);
    if (tkp2 <= _T + 100.0*std::numeric_limits<double>::epsilon())
      return tkp2;
    else
      return std::numeric_limits<double>::quiet_NaN();
  };

  /** get time instant k+3 of the time discretisation.
   *  It is used when we have to reschedule a TD Event in scheduleNonSmoothEvent
   *  \return a double. If the simulation is near the end (t_{k+3} >= T), it returns NaN.
  */
  inline double getTkp3() const
  {
    double tkp3 = _td->getTk(_k+3);
    if (tkp3 <= _T + 100.0*std::numeric_limits<double>::epsilon())
      return tkp3;
    else
      return std::numeric_limits<double>::quiet_NaN();
  };

  /** Get current timestep
   * \return the current timestep
   */
  inline double currentTimeStep()
  {
    return _td->currentTimeStep(_k);
  }

  /** get TimeDiscretisation
   * \return the TimeDiscretisation in use for the time integration
   */
  inline const TimeDiscretisation& timeDiscretisation() const { return *_td;};

  /** update final time
   * \param T the new final time
   * */
  inline void updateT(double T) { _T = T; };
};

#endif // EventsManager_H
