/* Siconos-Kernel, Copyright INRIA 2005-2013.
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
/*! \file ObserverEvent.hpp
  Observer Events
*/
#ifndef ObserverEvent_H
#define ObserverEvent_H

#include "Event.hpp"
#include "Observer.hpp"

/** Events when the observer updates the state estimate
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.6.0
 *  \date (Creation) May 21, 2013
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
  ObserverEvent(): Event(0.0, OBSERVER_EVENT) {};

public:

  /** constructor with time value as a parameter
   *  \param time the starting time of the Event
   *  \param name the type of the Event
   */
  ObserverEvent(double time, int name): Event(time, name) {};

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

  /** Increment time of the present event according to
      the time discretisation of the linked Actuator
  */
  void update();
};

#endif // ObserverEvent_H
