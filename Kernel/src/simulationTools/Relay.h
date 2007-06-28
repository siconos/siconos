/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
*/
#ifndef RELAY_H
#define RELAY_H

#include "OneStepNSProblem.h"

/** Non smooth problem written as a relay.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class Relay : public OneStepNSProblem
{
public:
  /** default constructor
  */
  Relay();

  /** xml constructor
  *  \param OneStepNSProblemXML* : the XML linked-object
  *  \param Simulation *: the simulation that owns the problem (optional)
  */
  Relay(OneStepNSProblemXML*, Simulation * = NULL);

  ~Relay();

  /** make the computation so solve the NS problem
  *  param double : current time
  */
  void compute(const double&);

  /** copy the data of the OneStepNSProblem to the XML tree
  *  \exception RuntimeException
  */
  void saveNSProblemToXML();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param OneStepNSProblem* : the one step problem which must be converted
  * \return a pointer on the problem if it is of the right type, NULL otherwise
  */
  static Relay* convert(OneStepNSProblem* osnsp);
};

#endif // RELAY_H
