/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "ExplicitLinearSMC.hpp"

#include "ControlSensor.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::control::ExplicitLinearSMC::ExplicitLinearSMC(std::shared_ptr<ControlSensor> sensor)
    : CommonSMC(ActuatorType::ExplicitLinearSMC, sensor) {}

siconos::control::ExplicitLinearSMC::ExplicitLinearSMC(
    std::shared_ptr<ControlSensor> sensor, std::shared_ptr<siconos::algebra::SiconosMatrix> B)
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
  *_sigma = *_Csurface * _sensor->y();
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
