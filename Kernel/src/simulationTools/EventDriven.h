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
#ifndef EVENTDRIVEN_H
#define EVENTDRIVEN_H

#include "Strategy.h"

/** \class EventDriven
 *  \brief Strategy based on event driven method, ie events detection (see theoretical manual for more details).
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 */

class EventDriven : public Strategy
{
public:

  /** \fn EventDriven()
   *  \brief defaut constructor
   */
  EventDriven();

  /** \fn EventDriven(StrategyXML*, Model*)
   *  \brief constructor with XML object of the EventDriven
   *  \param StrategyXML* : the XML object corresponding
   *  \param Model* : the Model which contains the Strategy
   */
  EventDriven(StrategyXML*, Model*);

  ~EventDriven();

  /** \fn EventDriven* convert (Strategy* str)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Strategy* : the Strategy which must be converted
   * \return a pointer on the Strategy if it is of the right type, NULL otherwise
   */
  static EventDriven* convert(Strategy* str);

};

#endif // EVENTDRIVEN_H
