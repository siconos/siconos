/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
#ifndef TIMEDISCREVENT_H
#define TIMEDISCREVENT_H

/** \class TimeDiscrEvent
 *  \brief class derived from Event one: events that corresponds to original time discretisation steps.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.1.
 *  \date (Creation) February 21, 2006
 *
 */

#include "Event.h"

class TimeDiscrEvent : public Event
{

private:

  /** Default constructor */
  TimeDiscrEvent();

public:

  /** \fn TimeDiscrEvent(const TimeDiscrEvent&)
   *  \brief copy constructor
   *  \param the timeDiscrEvent to be copied
   */
  TimeDiscrEvent(const TimeDiscrEvent&);

  /** \fn TimeDiscrEvent(const unsigned long int&)
   *  \brief constructor with time value as a parameter
   *  \param an unsigned long int
   */
  TimeDiscrEvent(const unsigned long int&);

  /** \fn ~TimeDiscrEvent()
   *  \brief destructor
   */
  ~TimeDiscrEvent();

  /** \fn void process()
   *  \brief
   */
  void process();
};

#endif // TimeDiscrEvent_H
