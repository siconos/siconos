// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "RegularTwisting.hpp"

siconos::control::RegularTwisting::RegularTwisting(std::shared_ptr<ControlSensor> sensor,
                                                   double gain, double beta)
    : Twisting(sensor, gain, beta, 0.) {
  _type = ActuatorType::RegularTwisting;
}

siconos::control::RegularTwisting::RegularTwisting(std::shared_ptr<ControlSensor> sensor)
    : Twisting(sensor, 0.) {
  _type = ActuatorType::RegularTwisting;
}
