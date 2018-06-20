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
/*! \file ObserverEvent.hpp
  Observer Events
*/
#ifndef ObserverEvent_H
#define ObserverEvent_H

#include "Event.hpp"
#include "SiconosControlFwd.hpp"
#include "ControlTypeDef.hpp"

/** Events when the observer updates the state estimate
 *
 *
 */
class ObserverEvent : public Event
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(ObserverEvent);


  /** The observer linked to the present event */
  SP::Observer _observer;

  /** Default constructor */
  ObserverEvent(): Event(0.0, OBSERVER_EVENT, true) {};

public:

  /** constructor with time value as a parameter
   *  \param time the starting time of the Event
   *  \param name the type of the Event
   */
  ObserverEvent(double time, int name): Event(time, name, true) {};

  /** destructor
   */
  ~ObserverEvent() {};

  /** get the Observer linked to this Event
   *  \return a SP::Observer to the Observer
   */
  inline SP::Observer observer() const
  {
    return _observer;
  };

  /** set the Observer linked to this Event
   *  \param newObserver the SP::Observer
   */
  void setObserverPtr(SP::Observer newObserver)
  {
    _observer = newObserver;
  };

  /** Call the capture method of the linked Observer
   *  \param sim a SP::Simulation (ignored).
   */
  void process(Simulation& sim);

};

#endif // ObserverEvent_H
