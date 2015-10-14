/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2011.
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

/*! \file SensorX.h
  A specific sensor, to capture q vector ("position") in Lagrangian systems.
  Used as an example on how to implement new user-sensors.
*/

#ifndef SensorX_H
#define SensorX_H

#include "Sensor.hpp"

class SiconosVector;

/** \class SensorX
 *  \brief Specific Sensor to get position
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) february 01, 2007
 *
 * A specific sensor, to capture q vector ("position") in Lagrangian systems.
 * Used as an example on how to implement new user-sensors.
 *
 */
class SensorX : public Sensor
{
private:

  /** Default constructor
   */
  SensorX();

  /** A copy vector of X at the last capture time
   */
  SP::SiconosVector storedX;

public:

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Sensor, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   */
  SensorX(int, SP::TimeDiscretisation);

  /** Destructor
   */
  ~SensorX();

  /** initialize sensor data.
   */
  void initialize();

  /** capture data when the SensorEvent is processed ( for example set data[SensorEvent]=... )
   */
  void capture();

  /** Encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Sensor*
   * \return a pointer on the derived type
   */
  static SensorX* convert(Sensor* s);
};

#endif
