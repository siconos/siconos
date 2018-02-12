/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file LinearSensor.hpp
 * A generic linear sensor, to capture the output y defined as y = Cx + Du
*/

#ifndef LinearSensor_H
#define LinearSensor_H

#include "ControlSensor.hpp"

/** \class LinearSensor
 *  \brief Common linear Sensor to get output of the system
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.3.0.
 *  \date (Creation) november 08, 2011
 *
 * A generic linear sensor, to capture the output y defined as y = Cx + Du
 *
 */
class LinearSensor : public ControlSensor
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LinearSensor);

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
  LinearSensor() {};

public:

  /** Constructor for the SensorFactory
   * \param ds the SP::DynamicalSystem it observes
   */
  LinearSensor(SP::DynamicalSystem ds);

  /** Constructor with the full set of data
   * \param ds the SP::DynamicalSystem it observes.
   * \param C a SP::SiconosMatrix.
   * \param D a SP::SiconosMatrix (optional).
   */
  LinearSensor(SP::DynamicalSystem ds,
      SP::SimpleMatrix C, SP::SimpleMatrix D = SP::SimpleMatrix());

  /** Destructor
   */
  virtual ~LinearSensor();

  /** initialize sensor data
   * \param m the Model
   */
  virtual void initialize(const NonSmoothDynamicalSystem& nsds);

  /** capture data when the SensorEvent is processed ( for example set data[SensorEvent]=... )
   */
  void capture();

  /** Set the C matrix.
   * \param C a SimpleMatrix
   */
  void setC(const SimpleMatrix& C);

  /** Set the C matrix
   * \param C a SP::SimpleMatrix
   */
  inline void setCPtr(SP::SimpleMatrix C)
  {
    _matC = C;
  };

  /** Set the D matrix
   * \param D a SimpleMatrix
   */
  void setD(const SimpleMatrix& D);

  /** Set the D matrix
   * \param D a SP::SimpleMatrix
   */
  inline void setDPtr(SP::SimpleMatrix D)
  {
    _matD = D;
  };
};

#endif
