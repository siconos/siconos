// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file Twisting.hpp
  \brief twisting algorithm with a modification
*/

#ifndef Twisting_H
#define Twisting_H

#include "CommonSMC.hpp"
namespace siconos::control {

/** Modified Twisting Controller with an implicit discretization
 *
 *  This is the twisting controller that gives the best result
 */
class Twisting : public CommonSMC {
 private:
  ACCEPT_SERIALIZATION(Twisting);

 public:
  /** Constructor for the ActuatorFactory
   *
   *  \param sensor the ControlSensor feeding the Actuator
   */
  Twisting(std::shared_ptr<ControlSensor> sensor)
      : CommonSMC(ActuatorType::Twisting, sensor) {};

  /** Constructor for a nonlinear system.
   *
   *  \param sensor the ControlSensor feeding the Actuator
   *  \param hControl sampling period
   */
  Twisting(std::shared_ptr<ControlSensor> sensor, double hControl);

  /** Constructor for the linear case
   *
   *  \param sensor the ControlSensor feeding the Actuator
   *  \param gain control magnitude
   *  \param beta twisting parameter
   *  \param hControl sampling period
   */
  Twisting(std::shared_ptr<ControlSensor> sensor, double gain, double beta, double hControl);

  /** destructor
   */
  virtual ~Twisting() noexcept = default;

  /** Compute the new control law at each event
   *  Here we are using the following formula:
   */
  virtual void actuate();

  /** set nonsmooth data: NormalConeNSL and AVI osnsp
   *
   *  \param hControl sampling period
   */
  virtual void setNSdata(double hControl);

  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);
};
// Register the observer into the factory
static ActuatorRegistration<Twisting> reg_AT(ActuatorType::Twisting);

}  // namespace siconos::control
#endif
