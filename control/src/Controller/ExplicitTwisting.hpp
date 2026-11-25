// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file ExplicitTwisting.hpp
  \brief twisting algorithm with an explicit discretization
*/

#ifndef ExplicitTwisting_H
#define ExplicitTwisting_H

#include "CommonSMC.hpp"

namespace siconos::control {

/** Twisting Controller with an explicit discretization */
class ExplicitTwisting : public CommonSMC {
 private:
  ACCEPT_SERIALIZATION(ExplicitTwisting);

 protected:
 public:
  /** Constructor for a nonlinear system or the ActuatorFactory
   * \param sensor the ControlSensor feeding the Actuator
   */
  ExplicitTwisting(std::shared_ptr<ControlSensor> sensor);

  /** Constructor for the linear case
   * \param sensor the ControlSensor feeding the Actuator
   * \param gain control magnitude
   * \param beta twisting parameter
   */
  ExplicitTwisting(std::shared_ptr<ControlSensor> sensor, double gain, double beta);

  /** destructor
   */
  virtual ~ExplicitTwisting() noexcept = default;

  /** Compute the new control law at each event
   *  Here we are using the following formula:
   */
  virtual void actuate();

  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);
};
// Register the observer into the factory
static ActuatorRegistration<ExplicitTwisting> reg_AET(ActuatorType::ExplicitTwisting);

}  // namespace siconos::control
#endif
