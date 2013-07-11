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
/*! \file
Time Discretisation Events
*/
#ifndef TIMEDISCRETISATIONEVENTNOSAVEINMEMORY_H
#define TIMEDISCRETISATIONEVENTNOSAVEINMEMORY_H

#include "Event.hpp"

/** Event that corresponds to user-defined time discretisation points
 *  This Event does not automatically save in memory some variables.
 *  Use it at your own risk
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.6.0.
 *  \date (Creation) July 10, 2013
 *
 */
class TimeDiscretisationEventNoSaveInMemory : public Event
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(TimeDiscretisationEventNoSaveInMemory);


  /** Default constructor */
  TimeDiscretisationEventNoSaveInMemory();

public:

  /** constructor with time value as a parameter
  *  \param time starting time (a double)
  *  \param notUsed unused int
  */
  TimeDiscretisationEventNoSaveInMemory(double time, int notUsed);

  /** destructor
  */
  ~TimeDiscretisationEventNoSaveInMemory();

  /**
  *  \param simulation the simulation that owns this Event (through the EventsManager)
  */
  void process(Simulation& simulation);
};

#endif // TimeDiscretisationEventNoSaveInMemory_H
