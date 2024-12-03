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

#include "LinearSMCOT2.hpp"

#include "ControlSensor.hpp"
#include "EventDriven.hpp"
#include "FirstOrderLinearDS.hpp"
#include "LsodarOSI.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"
#include "TimeDiscretisation.hpp"

siconos::control::LinearSMCOT2::LinearSMCOT2(std::shared_ptr<ControlSensor> sensor)
    : CommonSMC(ActuatorType::LinearSMCOT2, sensor) {}

void siconos::control::LinearSMCOT2::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s) {
  if (!_Csurface) {
    THROW_EXCEPTION(
        "CommonSMC::initialize - you have to set either _Csurface or h(.) before initializing "
        "the Actuator");
  } else {
    if (_Csurface && !_u)
      _u = std::make_shared<siconos::algebra::SiconosVector>(_Csurface->size(0), 0);
  }

  Actuator::initialize(nsds, s);

  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // We can use the Visitor mighty power to check if we have the right type
  auto DS = _sensor->getDS();

  if (auto folds = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearTIDS>(DS)) {
    _DSPhi = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(*folds);
    _DSPred = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(*folds);
  } else if (auto fods =
                 std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(DS)) {
    auto x0 = std::make_shared<siconos::algebra::SiconosVector>(*fods->x0());
    _DSPhi = std::make_shared<siconos::modeling::FirstOrderLinearDS>(x0);
    _DSPred = std::make_shared<siconos::modeling::FirstOrderLinearDS>(x0);
    if (fods->A()) {
      _DSPhi->setA(*fods->A());
      _DSPred->setA(*fods->A());
    }
    if (fods->b()) {
      _DSPhi->setb(*fods->b());
      _DSPred->setb(*fods->b());
    }
  } else
    THROW_EXCEPTION("LinearSMCOT2 implemented only for first order systems");

  // We have to reset _pluginb
  // _DSPhi->setComputebFunction(NULL);
  // _DSPred->setComputebFunction(NULL);
  // XXX What if there is more than one sensor ...

  _indx = 0;
  //  _Phi= std::make_shared<SiconosMatrix(_nDim, _nDim));
  //  _Phi->setIdentity();
  //  _Xold= std::make_shared<siconos::algebra::SiconosVector>(_nDim));
  //  *_Xold = *(_sensor->y());
  auto _t0 = nsds.t0();
  auto _T = nsds.finalT() + _tdPhi->timeStep(0);

  _XPhi = _DSPhi->x();

  _Xhat = _DSPred->x();
  auto dummyb = std::make_shared<siconos::algebra::SiconosVector>(_B->size(0), 0);
  _DSPred->setbPtr(dummyb);
  siconos::algebra::prod(*_B, *_u, *_DSPred->b());

  //  _Xhat= std::make_shared<siconos::algebra::SiconosVector>(_nDim, 0));
  //  _DSPred->setXPtr(_Xhat);

  _nsdsPhi = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  _PhiOSI = std::make_shared<siconos::integrators::LsodarOSI>();
  _nsdsPhi->insertDynamicalSystem(_DSPhi);
  _simulPhi = std::make_shared<siconos::simulation::EventDriven>(_nsdsPhi, _tdPhi, 0);
  _simulPhi->associate(_PhiOSI, _DSPhi);

  // Integration for Gamma
  _nsdsPred = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  _PredOSI = std::make_shared<siconos::integrators::LsodarOSI>();
  _nsdsPred->insertDynamicalSystem(_DSPred);
  _simulPred = std::make_shared<siconos::simulation::EventDriven>(_nsdsPred, _tdPred, 0);
  _simulPred->associate(_PredOSI, _DSPred);

  _X = _sensor->yTk();
}

void siconos::control::LinearSMCOT2::actuate() {
  auto hCurrent = _tdPhi->timeStep(_indx);
  // Get current value of the state
  // Update it
  *_XPhi = *_X;

  // We change the values of the state each time, so we need to change istate to 1
  // See LsodarOSI.cpp for the meaning of istate
  if (_indx > 0) {
    _simulPhi->setIstate(1);
    _simulPred->setIstate(1);
  }
  // Compute _XPhi = \Phi*X
  _simulPhi->advanceToEvent();
  _simulPhi->processEvents();
  // XXX small hack here
  auto CS = std::make_shared<siconos::algebra::SiconosVector>(_B->size(0));
  *CS = _Csurface->row(0);
  _coeff = -1 / (CS->vector_sum() * hCurrent);

  double uEq = CS->dot(_coeff * (*_XPhi + *_X - *_Xhat));
  double uEqP;
  // We need to project
  // TODO this should work in more than 1D
  uEqP = std::min(uEq, 2.0);
  uEqP = std::max(uEqP, -2.0);
  (*_u)(_u->size() - 1) = uEqP;
  siconos::algebra::prod(*_B, *_u, *_DSPred->b());
  _indx++;
  *_Xhat = *_X;
  // Compute \hat{x}_k
  _simulPred->advanceToEvent();
  _simulPred->processEvents();
}

void siconos::control::LinearSMCOT2::setTimeDiscretisation(
    const siconos::simulation::TimeDiscretisation& td) {
  _tdPhi = std::make_shared<siconos::simulation::TimeDiscretisation>(td);
  _tdPred = std::make_shared<siconos::simulation::TimeDiscretisation>(td);
}
