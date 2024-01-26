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

#include "ExplicitLinearSMC.hpp"

#include "ControlSensor.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

siconos::control::ExplicitLinearSMC::ExplicitLinearSMC(std::shared_ptr<ControlSensor> sensor)
    : CommonSMC(ActuatorType::ExplicitLinearSMC, sensor) {}

siconos::control::ExplicitLinearSMC::ExplicitLinearSMC(
    std::shared_ptr<ControlSensor> sensor, std::shared_ptr<siconos::algebra::SimpleMatrix> B)
    : CommonSMC(ActuatorType::ExplicitLinearSMC, sensor, B) {}

void siconos::control::ExplicitLinearSMC::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s) {
  CommonSMC::initialize(nsds, s);

  _sigma = std::make_shared<siconos::algebra::SiconosVector>(_u->size());
}

void siconos::control::ExplicitLinearSMC::actuate() {
  if (!_noUeq) {
    computeUeq();
  }

  siconos::algebra::prod(*_Csurface, _sensor->y(), *_sigma);

  auto sDim = _u->size();

  if (_D)  // we are using a saturation
  {
    for (decltype(sDim) i = 0; i < sDim; i++) {
      if ((*_sigma)(i) > (*_D)(i, i))
        (*_us)(i) = -_alpha;
      else if ((*_sigma)(i) < -(*_D)(i, i))
        (*_us)(i) = _alpha;
      else {
        if ((*_D)(i, i) != 0)
          (*_us)(i) = -(*_sigma)(i) / (*_D)(i, i);
        else
          (*_us)(i) = 0;
      }
    }
  } else {
    for (decltype(sDim) i = 0; i < sDim; i++) {
      if ((*_sigma)(i) > 0)
        (*_us)(i) = -_alpha;
      else if ((*_sigma)(i) < 0)
        (*_us)(i) = _alpha;
      else
        (*_us)(i) = 0;
    }
  }

  *_lambda = *_us;
  *_u = *_us;
  *_u += *_ueq;
}
