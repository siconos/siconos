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

#include "LinearSMC.hpp"

#include "ControlSensor.hpp"
#include "FirstOrderLinearDS.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "TimeStepping.hpp"
// #define DEBUG_WHERE_MESSAGES
//  #define DEBUG_NOCOLOR
//  #define DEBUG_STDOUT
//  #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::control::LinearSMC::LinearSMC(std::shared_ptr<ControlSensor> sensor,
                                       ActuatorType type)
    : CommonSMC(type, sensor) {}

siconos::control::LinearSMC::LinearSMC(std::shared_ptr<ControlSensor> sensor,
                                       std::shared_ptr<siconos::algebra::SimpleMatrix> B,
                                       std::shared_ptr<siconos::algebra::SimpleMatrix> D,
                                       ActuatorType type)
    : CommonSMC(type, sensor, B, D) {}

void siconos::control::LinearSMC::actuate() {
  DEBUG_BEGIN("void siconos::control::LinearSMC::actuate()\n")

  if (!_noUeq) {
    computeUeq();
    auto& LinearDS_SMC =
        *std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC);
    siconos::algebra::prod(*_B, *_ueq, *(LinearDS_SMC.b()));
  }

  DEBUG_EXPR(_DS_SMC->xMemory().display(););

  *(_DS_SMC->x()) = _sensor->y();

  // SS: Really need to modify stored xMemory?
  _DS_SMC->xMemory().getSiconosVectorMutable(0) = _sensor->y();

  if (not std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC)) {
    _DS_SMC->computef(_simulationSMC->startingTime(), _DS_SMC->x());
    _DS_SMC->swapInMemory();
    //    _DS_SMC->computef(_simulationSMC->startingTime());
    //    *_DS_SMC->fold() = *_DS_SMC->f();
  }

  _simulationSMC->computeOneStep();
  //  if (_indx > 0)
  { _simulationSMC->nextStep(); }

  // discontinous part
  *_us = *_lambda;
  *_u = *_us;
  *_u += *_ueq;
  DEBUG_EXPR(_u->display(););
  _indx++;
  DEBUG_END("void siconos::control::LinearSMC::actuate()\n")
}
