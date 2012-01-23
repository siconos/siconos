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

/*! \file linearSensor.hpp
 * A generic linear sensor, to capture the output y defined as y = Cx + Du
*/

#ifndef linearSensor_H
#define linearSensor_H

#include "SiconosKernel.hpp"

class SiconosMatrix;
class SimpleMatrix;
/** \class linearSensor
 *  \brief Common linear Sensor to get output of the system
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.3.0.
 *  \date (Creation) november 08, 2011
 *
 * A generic linear sensor, to capture the output y defined as y = Cx + Du
 *
 */
class linearSensor : public controlSensor
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(linearSensor);


  /** A matrix for output */
  SP::SiconosMatrix _data;
  /** A matrix for saving all values */
  SP::SimpleMatrix _dataPlot;
  /** counter */
  unsigned int _k;

  /** Canonical matrices */
  SP::SimpleMatrix _matC;
  SP::SimpleMatrix _matD;

  /** Number of time steps*/
  unsigned int _nSteps;

  /** Default constructor
   */
  linearSensor() {};

public:

  /** Constructor with a TimeDiscretisation and a Model.
   * \param an int, the type of the Sensor, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   */
  linearSensor(int, SP::TimeDiscretisation, SP::Model);

  /** Constructor with a TimeDiscretisation, a Model and two matrices.
   * \param an int, the type of the Sensor, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   * \param a SP::SiconosMatrix C
   * \param a SP::SiconosMatrix D (optional)
   */
  linearSensor(int, SP::TimeDiscretisation, SP::Model, SP::SimpleMatrix, SP::SimpleMatrix);


  /** Destructor
   */
  ~linearSensor();

  /** initialize sensor data.
   */
  void initialize();

  /** capture data when the SensorEvent is processed ( for example set data[SensorEvent]=... )
   */
  void capture();

  /** Set the C matrix.
   * \param a SimpleMatrix
   */
  void setC(const SimpleMatrix& C)
  {
    *_matC = C;
  };

  /** Set the C matrix
   * \param a SP::SimpleMatrix
   */
  void setCPtr(SP::SimpleMatrix C)
  {
    _matC = C;
  };

  /** Set the D matrix
   * \param a SimpleMatrix
   */
  void setD(const SimpleMatrix& D)
  {
    *_matD = D;
  };

  /** Set the D matrix
   * \param a SP::SimpleMatrix
   */
  void setDPtr(SP::SimpleMatrix D)
  {
    _matD = D;
  };
};
TYPEDEF_SPTR(linearSensor)
#endif
