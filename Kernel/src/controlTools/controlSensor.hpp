/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file controlSensor.hpp
 * A generic control sensor
*/

#ifndef controlSensor_H
#define controlSensor_H

#include "SiconosKernel.hpp"

class SiconosMatrix;
class SimpleMatrix;
/** \class controlSensor
 *  \brief Generic control Sensor to get the output of the system
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.3.0.
 *  \date (Creation) november 09, 2011
 *
 * A generic control sensor
 *
 */
class controlSensor : public Sensor
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(controlSensor);

protected:
  /** Dimension of the output */
  unsigned int _YDim;
  /** A vector for the current value of the output */
  SP::SimpleVector _storedY;

public:

  controlSensor(int name, SP::TimeDiscretisation t, SP::Model m): Sensor(name, t, m) {}
  /** Get the dimension of the output
   * \return an unsigned int
   */
  inline unsigned int getYDim() const
  {
    return _YDim;
  };

  /** Get a pointer to the output
   * \return pointer
   */
  inline SP::SimpleVector y() const
  {
    return _storedY;
  };

  /** Encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Sensor*
   * \return a pointer on the derived type
   */
  static controlSensor* convert(Sensor* s);
};
DEFINE_SPTR(controlSensor)
#endif
