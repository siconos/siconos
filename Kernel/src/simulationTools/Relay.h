/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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

/** \class Relay
 *  \brief Non smooth problem written as a relay.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.3.
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
#ifndef RELAY_H
#define RELAY_H

#include "OneStepNSProblem.h"

class Relay : public OneStepNSProblem
{
public:
  /** \fn Relay()
   *  \brief default constructor
   */
  Relay();

  /** \fn Relay(OneStepNSProblemXML*, Strategy*=NULL)
   *  \brief xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Strategy *: the strategy that owns the problem (optional)
   */
  Relay(OneStepNSProblemXML*, Strategy * = NULL);

  ~Relay();

  /** \fn void compute(const double&)
   *  \brief make the computation so solve the NS problem
   *  param double : current time
   */
  void compute(const double&);

  /** \fn void saveNSProblemToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSProblemToXML();

  /** \fn Relay* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static Relay* convert(OneStepNSProblem* osnsp);
};

#endif // RELAY_H
