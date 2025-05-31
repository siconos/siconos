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

#include "LinearSMCimproved.hpp"

#include <boost/circular_buffer.hpp>

#include "ControlSensor.hpp"
#include "FirstOrderLinearDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"
#include "TimeStepping.hpp"
#include "ZeroOrderHoldOSI.hpp"

siconos::control::LinearSMCimproved::LinearSMCimproved(std::shared_ptr<ControlSensor> sensor)
    : LinearSMC(sensor, ActuatorType::LinearSMCimproved) {}

siconos::control::LinearSMCimproved::LinearSMCimproved(
    std::shared_ptr<ControlSensor> sensor, std::shared_ptr<siconos::algebra::SiconosMatrix> B,
    std::shared_ptr<siconos::algebra::SiconosMatrix> D)
    : LinearSMC(sensor, B, D, ActuatorType::LinearSMCimproved) {}

void siconos::control::LinearSMCimproved::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s) {
  LinearSMC::initialize(nsds, s);
  _up = std::make_shared<siconos::algebra::SiconosVector>(_us->size());
  _measuredPert = std::make_shared<
      boost::circular_buffer<std::shared_ptr<siconos::algebra::SiconosVector>>>(0);
  _predictedPert = std::make_shared<
      boost::circular_buffer<std::shared_ptr<siconos::algebra::SiconosVector>>>(0);
}

void siconos::control::LinearSMCimproved::predictionPerturbation(
    const siconos::algebra::SiconosVector& xTk,
    const Eigen::FullPivLU<siconos::algebra::SiconosMatrix>& LUCBstar) {
  if (siconos::algebra::normInf(*_us) < _alpha) {
    if (_inDisceteTimeSlidingPhase) {
      if (_measuredPert->full()) {
        if (_measuredPert->size() > 1) {
          _measuredPert->rotate(_measuredPert->end() - 1);
          _predictedPert->rotate(_predictedPert->end() - 1);
        }
      } else {
        // inject new vector in the case where the measurement vector is not full.
        auto sp1 = std::make_shared<siconos::algebra::SiconosVector>(_us->size(), 0);
        auto sp2 = std::make_shared<siconos::algebra::SiconosVector>(_us->size(), 0);
        _measuredPert->push_front(sp1);
        _predictedPert->push_front(sp2);
      }

      // inject new measured value and also perturbation prediction
      auto& predictedPertC = *(*_predictedPert)[0];
      auto& measuredPertC = *(*_measuredPert)[0];

      // Cp_k = s_k + Cp_k-tilde
      siconos::algebra::prod(*_Csurface, xTk, measuredPertC);
      measuredPertC += *(*_predictedPert)[std::min((unsigned int)1,
                                                   (unsigned int)_predictedPert->size() - 1)];

      // compute prediction
      switch (_measuredPert->size() - 1) {
        case 0:
          predictedPertC = measuredPertC;
          break;
        case 1:
          predictedPertC = 2 * measuredPertC - *(*_measuredPert)[1];
          break;
        case 2:
          predictedPertC =
              3 * measuredPertC - 3 * (*(*_measuredPert)[1]) + *(*_measuredPert)[2];
          break;
        default:
          THROW_EXCEPTION(
              "siconos::control::LinearSMCimproved::predictionPerturbation: unknown order " +
              std::to_string(_measuredPert->size()));
      }

      // Compute the control to counteract the perturbation
      *_up = predictedPertC;
      *_up *= -1;
      *_up = LUCBstar.solve(*_up);

      // project onto feasible set
      double norm = _up->norm();
      if (norm > _ubPerturbation) {
        *_up *= _ubPerturbation / norm;
        predictedPertC *= _ubPerturbation / norm;
      }
    } else
      _inDisceteTimeSlidingPhase = true;
  } else if (_inDisceteTimeSlidingPhase) {
    _inDisceteTimeSlidingPhase = false;
    _up->setZero();
  }
}

void siconos::control::LinearSMCimproved::actuate() {
  auto xTk = std::make_shared<siconos::algebra::SiconosVector>(_sensor->y());

  auto& zoh =
      *std::static_pointer_cast<siconos::integrators::ZeroOrderHoldOSI>(_integratorSMC);

  // equivalent part
  auto tmpM1 = -*_Csurface * zoh.Ad(_DS_SMC) + *_Csurface;
  auto CBstar = *_Csurface * zoh.Bd(_DS_SMC);

  // compute C(I-e^{Ah})x_k
  *_ueq = tmpM1 * *xTk;
  // compute the solution u^eq of the system CB^{*}u^eq = C(I-e^{Ah})x_k
  Eigen::FullPivLU<siconos::algebra::SiconosMatrix> luCBstar(CBstar);
  *_ueq = luCBstar.solve(*_ueq);
  *(_DS_SMC->x()) = *xTk;

  bSMC_ = *_B * *_ueq;  // update DS_SMC b vector because of shared memory view
  _simulationSMC->computeOneStep();
  _simulationSMC->nextStep();

  *_us = *_lambda;

  // inject those in the system
  *_u = *_us;
  *_u += *_ueq;

  // prediction of the perturbation
  if (_predictionPerturbation) {
    predictionPerturbation(*xTk, luCBstar);
    *_u += *_up;
  }

  _indx++;
}

void siconos::control::LinearSMCimproved::setPredictionOrder(unsigned int order) {
  _measuredPert->set_capacity(order + 1);
  _predictedPert->set_capacity(order + 1);
}
