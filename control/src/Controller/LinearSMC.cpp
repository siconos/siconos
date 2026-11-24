// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LinearSMC.hpp"

#include "ControlSensor.hpp"
#include "FirstOrderLinearDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "TimeStepping.hpp"  // IWYU pragma: keep
//  #define DEBUG_WHERE_MESSAGES
//   #define DEBUG_NOCOLOR
//   #define DEBUG_STDOUT
//   #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::control::LinearSMC::LinearSMC(std::shared_ptr<ControlSensor> sensor,
                                       ActuatorType type)
    : CommonSMC(type, sensor) {}

siconos::control::LinearSMC::LinearSMC(std::shared_ptr<ControlSensor> sensor,
                                       std::shared_ptr<siconos::algebra::SiconosMatrix> B,
                                       std::shared_ptr<siconos::algebra::SiconosMatrix> D,
                                       ActuatorType type)
    : CommonSMC(type, sensor, B, D) {}

void siconos::control::LinearSMC::actuate() {
  DEBUG_BEGIN("void siconos::control::LinearSMC::actuate()\n")

  if (!_noUeq) {
    computeUeq();
    auto& LinearDS_SMC =
        *std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC);
    bSMC_.resize(_ueq->size());
    bSMC_ = *_B * *_ueq;
    LinearDS_SMC.setConstantbVector(bSMC_, siconos::algebra::alias_t);  // Shared memory view
  }

  DEBUG_EXPR(siconos::algebra::print(_DS_SMC->xMemory()););

  *(_DS_SMC->x()) = _sensor->y();

  // SS: Really need to modify stored xMemory?
  _DS_SMC->xMemory().getSiconosVectorMutable(0) = _sensor->y();

  if (not std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC)) {
    _DS_SMC->computefVector(*_DS_SMC->x(), _simulationSMC->startingTime());
    _DS_SMC->swapInMemory();
  }

  _simulationSMC->computeOneStep();
  //  if (_indx > 0)
  {
    _simulationSMC->nextStep();
  }

  // discontinous part
  *_us = *_lambda;
  *_u = *_us;
  *_u += *_ueq;
  DEBUG_EXPR(siconos::algebra::print(*_u));
  ;
  _indx++;
  DEBUG_END("void siconos::control::LinearSMC::actuate()\n")
}
