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

#include "Twisting.hpp"

#include <iostream>

#include "AVI.hpp"
#include "ControlSensor.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "NormalConeNSL.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "TimeStepping.hpp"

siconos::control::Twisting::Twisting(std::shared_ptr<ControlSensor> sensor, double gain,
                                     double beta, double hControl)
    : CommonSMC(ActuatorType::Twisting, sensor)
{
  _u = std::make_shared<siconos::algebra::SiconosVector>(2);
  if (beta <= 0.0 || beta >= 1.0) {
    std::cout << "Twisting constructor: beta is not in (0, 1)" << std::endl;
  }

  _B = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2);
  (*_B)(1, 0) = gain;
  (*_B)(1, 1) = gain * beta;

  setNSdata(hControl);
}

siconos::control::Twisting::Twisting(std::shared_ptr<ControlSensor> sensor, double hControl)
    : CommonSMC(ActuatorType::Twisting, sensor)
{
  _u = std::make_shared<siconos::algebra::SiconosVector>(2);
  setNSdata(hControl);
}

void siconos::control::Twisting::setNSdata(double hControl)
{
  std::shared_ptr<siconos::algebra::SimpleMatrix> H =
      std::make_shared<siconos::algebra::SimpleMatrix>(4, 2);
  (*H)(0, 0) = 1.0;
  (*H)(1, 0) = -hControl / 2.0;
  (*H)(2, 0) = -1.0;
  (*H)(3, 0) = hControl / 2.0;
  (*H)(1, 1) = 1.0;
  (*H)(3, 1) = -1.0;

  auto K = std::make_shared<siconos::algebra::SiconosVector>(4);
  (*K)(0) = -1.0;
  (*K)(1) = -1.0;
  (*K)(2) = -1.0;
  (*K)(3) = -1.0;

  _nsLawSMC = std::make_shared<siconos::modeling::NormalConeNSL>(2, H, K);
  _OSNSPB_SMC = std::make_shared<siconos::nonsmooth_formulations::AVI>();
  _numericsSolverId = SICONOS_AVI_CAOFERRIS;
}

void siconos::control::Twisting::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s)
{
  // basic check
  if (!_nsLawSMC || !_OSNSPB_SMC) {
    THROW_EXCEPTION(
        "siconos::control::Twisting::initialize - nslaw or osnsp not set. If you used the "
        "constructor with only the ControlSensor as argument, you need to manually call "
        "setNSdata");
  }

  CommonSMC::initialize(nsds, s);
}

void siconos::control::Twisting::actuate()
{
  *(_DS_SMC->x()) = _sensor->y();

  _simulationSMC->computeOneStep();
  _simulationSMC->nextStep();

  // discontinous part
  *_us = *_lambda;
  *_u = *_us;
  _indx++;
}
