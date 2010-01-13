/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
/*! \file
Time Discretisation Events
*/
#ifndef TIMEDISCREVENT_H
#define TIMEDISCREVENT_H

#include "Event.hpp"

/** Event that corresponds to user-defined time discretisation points
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 21, 2006
 *
 */
class TimeDiscretisationEvent : public Event
{

private:

  /** Default constructor */
  TimeDiscretisationEvent();

public:

  /** constructor with time value as a parameter
  *  \param a double
  *  \param an int, type of Event
  */
  TimeDiscretisationEvent(double, int);

  /** destructor
  */
  ~TimeDiscretisationEvent();

  /**
  *  \param SP::Simulation, the simulation that owns this Event (through the EventsManager)
  */
  void process(SP::Simulation);
};

#endif // TimeDiscretisationEvent_H
