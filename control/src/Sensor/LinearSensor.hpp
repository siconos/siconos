// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
  std::shared_ptr<siconos::algebra::SiconosMatrix> _dataPlot{nullptr};
  /** counter */
  size_t _k{0};

  /** Canonical matrices */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _matC{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _matD{nullptr};

  /** Number of time steps*/
  size_t _nSteps{0};

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
               std::shared_ptr<siconos::algebra::SiconosMatrix> C,
               std::shared_ptr<siconos::algebra::SiconosMatrix> D = nullptr);

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
   *  \param C a SiconosMatrix
   */
  void setC(const siconos::algebra::SiconosMatrix& C);

  /** Set the C matrix
   *
   *  \param C a std::shared_ptr<siconos::algebra::SiconosMatrix>
   */
  inline void setCPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> C) { _matC = C; };

  /** Set the D matrix
   *
   *  \param D a SiconosMatrix
   */
  void setD(const siconos::algebra::SiconosMatrix& D);

  /** Set the D matrix
   *
   *  \param D a std::shared_ptr<siconos::algebra::SiconosMatrix>
   */
  inline void setDPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> D) { _matD = D; };
};

// Register the sensor into the factory
static SensorRegistration<LinearSensor> reg_SL(SensorType::Linear);

}  // namespace siconos::control
#endif
