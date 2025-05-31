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

#include "PID.hpp"

#include <iostream>

#include "ControlSensor.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::control::PID::PID(std::shared_ptr<ControlSensor> sensor,
                           std::shared_ptr<siconos::algebra::SiconosMatrix> B)
    : Actuator(ActuatorType::PID, sensor, B)
{
  _u = std::make_shared<siconos::algebra::SiconosVector>(1);
  _u->setZero();
}

void siconos::control::PID::initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                                       const siconos::simulation::Simulation& s)
{
  Actuator::initialize(nsds, s);

  _curDeltaT = s.currentTimeStep();

  // initialize _err
  _err = std::make_shared<boost::circular_buffer<double>>(3);
  for (auto i = 0; i < 3; ++i) (*_err).push_front(0.0);
}

void siconos::control::PID::actuate()
{
  DEBUG_BEGIN("void siconos::control::PID::actuate()\n");
  /** \todo We have to distinguish two cases : linear or nonlinear
   *  support the nonlinear case
   */

  // Compute the new error

  (*_err).push_front(_ref - _sensor->y()(0));
  DEBUG_PRINTF("_curDeltaT = %g\n", _curDeltaT);
  DEBUG_PRINTF("_ref = %g\n", _ref);
  DEBUG_EXPR(siconos::algebra::print(_sensor->y()););
  DEBUG_EXPR(siconos::algebra::print(*_u));;
  DEBUG_PRINTF("added term  = %g\n",
               ((*_K)(0) + (*_K)(2) / _curDeltaT + (*_K)(1) * _curDeltaT) * (*_err)[0] +
                   (-(*_K)(0) - 2 * (*_K)(2) / _curDeltaT) * (*_err)[1] +
                   (*_K)(2) / _curDeltaT * (*_err)[2]);
  // compute the new control and update it
  (*_u)(0) += ((*_K)(0) + (*_K)(2) / _curDeltaT + (*_K)(1) * _curDeltaT) * (*_err)[0] +
              (-(*_K)(0) - 2 * (*_K)(2) / _curDeltaT) * (*_err)[1] +
              (*_K)(2) / _curDeltaT * (*_err)[2];
  DEBUG_EXPR(siconos::algebra::print(*_u));;
  DEBUG_END("void siconos::control::PID::actuate()\n");
}

void siconos::control::PID::setK(std::shared_ptr<siconos::algebra::SiconosVector> K)
{
  // check dimensions ...
  if (K->size() != 3) {
    THROW_EXCEPTION("siconos::control::PID::setK - the size of K should be 3");
  }
  else {
    _K = K;
  }
}

void siconos::control::PID::setTimeDiscretisation(
    const siconos::simulation::TimeDiscretisation& td)
{
}

void siconos::control::PID::display() const
{
  Actuator::display();
  std::cout << "current error vector: ";
  std::cout << (*_err)[0] << " " << (*_err)[1] << " " << (*_err)[2] << std::endl;
}
