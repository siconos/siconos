// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file RegularTwisting.hpp
  \brief twisting algorithm with a direct discretization
*/

#ifndef RegularTwisting_H
#define RegularTwisting_H

#include "Twisting.hpp"

namespace siconos::control {

/** Twisting Controller with a straightforward discretization */
class RegularTwisting : public Twisting {
 private:
  ACCEPT_SERIALIZATION(RegularTwisting);

 public:
  /** Constructor for a nonlinear system or the ActuatorFactory
   * \param sensor the ControlSensor feeding the Actuator
   * \param hControl sampling period
   */
  RegularTwisting(std::shared_ptr<ControlSensor> sensor);

  /** Constructor for the linear case
   * \param sensor the ControlSensor feeding the Actuator
   * \param gain control magnitude
   * \param beta twisting parameter
   * \param hControl sampling period
   */
  RegularTwisting(std::shared_ptr<ControlSensor> sensor, double gain, double beta);

  /** destructor
   */
  virtual ~RegularTwisting() noexcept = default;
};
// Register the observer into the factory
static ActuatorRegistration<RegularTwisting> reg_ART(ActuatorType::RegularTwisting);

}  // namespace siconos::control
#endif
