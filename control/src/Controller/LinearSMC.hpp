// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file LinearSMC.hpp
 *
 *  \brief Linear sliding mode controller
 */

#ifndef LinearSMC_H
#define LinearSMC_H

#include "CommonSMC.hpp"
namespace siconos::control {

/**
   Linear sliding mode controller

   This controller implements the following sliding mode control strategy:
   \f[
   \begin{equation} u = u_{\mathrm{eq}} + u_s\qquad\text{where}\begin{cases} u_{\mathrm{eq}}
   \text{is the equivalent control input}\\u_s = -sgn(\sigma)\end{cases},\end{equation} \f]

   where \f$ \sigma = Cx \f$ is the user-defined sliding variable.

*/
class LinearSMC : public CommonSMC {
 private:
  ACCEPT_SERIALIZATION(LinearSMC);

 public:
  /** Constructor for the ActuatorFactory
   *
   *  \param sensor the ControlSensor feeding the Actuator
   *  \param type do not set this yourself ! this is used in derived classes
   */
  LinearSMC(std::shared_ptr<ControlSensor> sensor,
            ActuatorType type = ActuatorType::LinearSMC);

  /** Constructor with all the data
   *
   *  \param sensor the ControlSensor feeding the Actuator
   *  \param B the matrix
   *  \param D the D matrix in the FirstOrderLinearR
   *  \param type do not set this yourself ! this is used in derived classes
   */
  LinearSMC(std::shared_ptr<ControlSensor> sensor,
            std::shared_ptr<siconos::algebra::SiconosMatrix> B,
            std::shared_ptr<siconos::algebra::SiconosMatrix> D = nullptr,
            ActuatorType type = ActuatorType::LinearSMC);

  /** destructor
   */
  virtual ~LinearSMC() noexcept = default;

  /** Compute the new control law at each event
   *  Here we are using the following formula:
   *  \f$ u = u_{\mathrm{eq}} + u_s \f$
   */
  virtual void actuate();

  /** Set the D matrix
   *
   *  \param D the new D matrix
   */
  inline void setDPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> D) { _D = D; };
};
// Register the observer into the factory
static ActuatorRegistration<LinearSMC> reg_ALSMC(ActuatorType::LinearSMC);

}  // namespace siconos::control
#endif
