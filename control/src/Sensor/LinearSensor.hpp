/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

namespace siconos::control {
/**
   A generic linear sensor, to capture the output y defined as y = Cx + Du
*/
class LinearSensor : public ControlSensor {
 private:
  ACCEPT_SERIALIZATION(LinearSensor);

  /** A matrix for output */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _data{nullptr};
  /** A matrix for saving all values */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _dataPlot{nullptr};
  /** counter */
  unsigned int _k{0};

  /** Canonical matrices */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _matC{nullptr};
  std::shared_ptr<siconos::algebra::SimpleMatrix> _matD{nullptr};

  /** Number of time steps*/
  unsigned int _nSteps{0};

  // /** Default constructor
  //  */
  // LinearSensor(){};

 public:
  /** Constructor for the SensorFactory
   *
   *  \param ds the std::shared_ptr<siconos::modeling::DynamicalSystem> it observes
   */
  LinearSensor(std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** Constructor with the full set of data
   *
   *  \param ds the std::shared_ptr<siconos::modeling::DynamicalSystem> it observes.
   *  \param C a std::shared_ptr<siconos::algebra::SiconosMatrix>.
   *  \param D a std::shared_ptr<siconos::algebra::SiconosMatrix> (optional).
   */
  LinearSensor(std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
               std::shared_ptr<siconos::algebra::SimpleMatrix> C,
               std::shared_ptr<siconos::algebra::SimpleMatrix> D = nullptr);

  /** Destructor
   */
  virtual ~LinearSensor() noexcept = default;

  /** initialize sensor data
   *  \param nsds current nonsmooth dynamical system
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** capture data when the SensorEvent is processed ( for example set data[SensorEvent]=... )
   */
  void capture();

  /** Set the C matrix.
   *
   *  \param C a SimpleMatrix
   */
  void setC(const siconos::algebra::SimpleMatrix& C);

  /** Set the C matrix
   *
   *  \param C a std::shared_ptr<siconos::algebra::SimpleMatrix>
   */
  inline void setCPtr(std::shared_ptr<siconos::algebra::SimpleMatrix> C) { _matC = C; };

  /** Set the D matrix
   *
   *  \param D a SimpleMatrix
   */
  void setD(const siconos::algebra::SimpleMatrix& D);

  /** Set the D matrix
   *
   *  \param D a std::shared_ptr<siconos::algebra::SimpleMatrix>
   */
  inline void setDPtr(std::shared_ptr<siconos::algebra::SimpleMatrix> D) { _matD = D; };
};

// Register the sensor into the factory
static SensorRegistration<LinearSensor> reg_SL(SensorType::Linear);

}  // namespace siconos::control
#endif
