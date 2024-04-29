
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "ExplicitTwisting.hpp"

#include <math.h>  // for copysign

#include <SiconosVector.hpp>
#include <SiconosMatrix.hpp>
#include <iostream>

#include "ControlSensor.hpp"

siconos::control::ExplicitTwisting::ExplicitTwisting(std::shared_ptr<ControlSensor> sensor,
                                                     double gain, double beta)
    : CommonSMC(ActuatorType::ExplicitTwisting, sensor)
{
  _u = std::make_shared<siconos::algebra::SiconosVector>(2);
  if (beta <= 0.0 || beta >= 1.0) {
    std::cout << "ExplicitTwisting constructor: beta is not in (0, 1)" << std::endl;
  }

  _B = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  (*_B)(1, 0) = gain;
  (*_B)(1, 1) = gain * beta;
}

siconos::control::ExplicitTwisting::ExplicitTwisting(std::shared_ptr<ControlSensor> sensor)
    : CommonSMC(ActuatorType::ExplicitTwisting, sensor)
{
  _u = std::make_shared<siconos::algebra::SiconosVector>(2);
}

void siconos::control::ExplicitTwisting::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s)
{
  // \TODO(Xhub) this is quite unnecessary
  CommonSMC::initialize(nsds, s);
}

void siconos::control::ExplicitTwisting::actuate()
{
  const auto& sigma = _sensor->y();

  // discontinous part
  _u->setValue(0, std::copysign(1., -sigma(0)));
  _u->setValue(1, std::copysign(1., -sigma(1)));
  *_us = *_u;
  _indx++;
}
