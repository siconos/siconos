/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
/*! \file
 Non-Smooth Events
*/
#ifndef NONSMOOTHEVENT_H
#define NONSMOOTHEVENT_H

/** Events due to non smooth behavior (contact occurence...)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 21, 2006
 *
 * Those events are detected during Simulation process (integration of the smooth part with a roots-finding algorithm)
 * and scheduled into the EventsManager.
 *
 */

#include "Event.h"

class NonSmoothEvent : public Event
{

private:

  /** Default constructor */
  NonSmoothEvent();

public:

  /** constructor with time value as a parameter
  *  \param a double
  *  \param an int, type of Event
  */
  NonSmoothEvent(double, int);

  /** destructor
  */
  ~NonSmoothEvent();

  /** OSNS solving and IndexSets updating
  *  \param SP::Simulation, the simulation that owns this Event (through the EventsManager)
  */
  void process(SP::Simulation);
};

#endif // NonSmoothEvent_H
