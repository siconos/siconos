// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "ExplicitTwisting.hpp"

#include <math.h>  // for copysign

#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>
#include <iostream>

#include "ControlSensor.hpp"

siconos::control::ExplicitTwisting::ExplicitTwisting(std::shared_ptr<ControlSensor> sensor,
                                                     double gain, double beta)
    : CommonSMC(ActuatorType::ExplicitTwisting, sensor) {
  _u = std::make_shared<siconos::algebra::SiconosVector>(2);
  _u->setZero();

  if (beta <= 0.0 || beta >= 1.0) {
    std::cout << "ExplicitTwisting constructor: beta is not in (0, 1)" << std::endl;
  }

  _B = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  _B->setZero();
  (*_B)(1, 0) = gain;
  (*_B)(1, 1) = gain * beta;
}

siconos::control::ExplicitTwisting::ExplicitTwisting(std::shared_ptr<ControlSensor> sensor)
    : CommonSMC(ActuatorType::ExplicitTwisting, sensor) {
  _u = std::make_shared<siconos::algebra::SiconosVector>(2);
  _u->setZero();
}

void siconos::control::ExplicitTwisting::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s) {
  // \TODO(Xhub) this is quite unnecessary
  CommonSMC::initialize(nsds, s);
}

void siconos::control::ExplicitTwisting::actuate() {
  const auto& sigma = _sensor->y();

  // discontinous part
  (*_u)(0) = std::copysign(1., -sigma(0));
  (*_u)(1) = std::copysign(1., -sigma(1));
  *_us = *_u;
  _indx++;
}
