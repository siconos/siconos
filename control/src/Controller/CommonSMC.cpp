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

#include "CommonSMC.hpp"

#include <Eigen/LU>

#include "ControlSensor.hpp"
#include "EulerMoreauOSI.hpp"
#include "FirstOrderLinearDS.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderNonLinearR.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relay.hpp"
#include "RelayNSL.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "SolverOptions.h"
#include "StorageTools.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "ZeroOrderHoldOSI.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::control::CommonSMC::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s) {
  DEBUG_BEGIN(
      "siconos::control::CommonSMC::initialize(const "
      "siconos::modeling::NonSmoothDynamicalSystem & nsds, const "
      "Simulation & s)\n");
  if (!_Csurface && !computeh_) {
    THROW_EXCEPTION(
        "siconos::control::CommonSMC::initialize - you have to set either _Csurface or h(.) "
        "before initializing "
        "the Actuator");
  } else {
    if (_Csurface && !_u) {
      _u = std::make_shared<siconos::algebra::SiconosVector>(_Csurface->rows());
      _u->setZero();
    }

    Actuator::initialize(nsds, s);
  }
  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS
  // We can use the Visitor mighty power to check if we have the right type
  auto DS = _sensor->getDS();
  // create the DS for the controller
  // if the DS we use is different from the DS we are controlling
  // when we want for instant to see how well the controller behaves
  // if the plant model is not exact, we can use the setSimulatedDS
  // method
  if (auto folds = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(DS)) {
    // Copy but use a different b vector
    _DS_SMC = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*folds);
    bSMC_.resize(_DS_SMC->dimension());
    bSMC_.setZero();
    std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC)
        ->setConstantbVector(bSMC_,
                             siconos::algebra::alias_t);  // Shared memory view, bSMC_ is DS.b
  } else if (auto fonlds =
                 std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(DS)) {
    _DS_SMC = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(*fonlds);
  } else {
    THROW_EXCEPTION("LinearSMC is only  implemented for FirstOrder DS");
  }
  _DS_SMC->setNumber(999999);
  _DS_SMC->initMemory(1);
  _DS_SMC->swapInMemory();

  // Get the dimension of the output
  // XXX What if there is more than one sensor ...
  double t0 = nsds.t0();
  double T = nsds.finalT() + _td->timeStep(0);
  // create the SMC Model
  _nsdsSMC = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
  // Set up the simulation
  _simulationSMC = std::make_shared<siconos::simulation::TimeStepping>(_nsdsSMC, _td);

  auto sDim = _u->size();
  // create the interaction

  // Note FP: it seems that the type of relation is deduced from the coeff or
  // plugins that are set ...

  if (isRelationLinear_)  // deduced from the setCompute... functions call
  {
    if (computee_) {
      DEBUG_PRINT("A FirstOrderLinearR is created for the _relationSMC\n");
      _relationSMC = std::make_shared<siconos::modeling::FirstOrderLinearR>();
      auto folr = std::static_pointer_cast<siconos::modeling::FirstOrderLinearR>(_relationSMC);
      folr->setConstantC(*_Csurface);
      folr->setConstantB(*_B);
      folr->setComputeeVectorFunction(computee_);
      folr->setConstantD(*_D);
    } else {
      DEBUG_PRINT("A FirstOrderLinearTIR is created for the _relationSMC\n");
      _relationSMC = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*_Csurface, *_B);
      auto foltir =
          std::static_pointer_cast<siconos::modeling::FirstOrderLinearTIR>(_relationSMC);
      if (_D) foltir->setConstantD(*_D);
    }
  } else {
    // Non linear relations. NonLinearR, Type1R, Type2R.
    // To simplify, let's consider only NonlinearR.
    DEBUG_PRINT("A FirstOrderNonLinearR is created for the _relationSMC\n");
    _relationSMC = std::make_shared<siconos::modeling::FirstOrderNonLinearR>();

    auto fonlr =
        std::static_pointer_cast<siconos::modeling::FirstOrderNonLinearR>(_relationSMC);
    if (_Csurface) fonlr->setConstantJacobianhOver_state(*_Csurface);
    if (computejacobianhOver_state_)
      fonlr->setComputeJacobianhOver_stateFunction(computejacobianhOver_state_);

    if (_D) fonlr->setConstantJacobianhOver_lambda(*_D);
    if (computejacobianhOver_lambda_)
      fonlr->setComputeJacobianhOver_lambdaFunction(computejacobianhOver_lambda_);

    if (computejacobiangOver_state_)
      fonlr->setComputeJacobiangOver_stateFunction(computejacobiangOver_state_);

    if (_B) fonlr->setConstantJacobiangOver_lambda(*_B);
    if (computejacobiangOver_lambda_)
      fonlr->setComputeJacobiangOver_lambdaFunction(computejacobiangOver_lambda_);

    if (computeh_) fonlr->setComputehFunction(computeh_);
    if (computeg_) fonlr->setComputegFunction(computeg_);

    if (_computeResidus) {
      _simulationSMC->setComputeResiduY(true);
      _simulationSMC->setComputeResiduR(true);
      _simulationSMC->setUseRelativeConvergenceCriteron(false);
    }
  }

  // _nsLawSMC and the OSNSP can be defined in derived classes, like twisting
  if (!_nsLawSMC)
    _nsLawSMC = std::make_shared<siconos::modeling::RelayNSL>(sDim, -_alpha, _alpha);
  if (!_OSNSPB_SMC)
    _OSNSPB_SMC = std::make_shared<siconos::nonsmooth_formulations::Relay>(_numericsSolverId);

  _interactionSMC = std::make_shared<siconos::modeling::Interaction>(_nsLawSMC, _relationSMC);

  if (std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(DS))
    _integratorSMC = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();
  else if (std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(DS))
    _integratorSMC = std::make_shared<siconos::integrators::EulerMoreauOSI>(_thetaSMC);
  else
    THROW_EXCEPTION("LinearSMC is only  implemented for FirstOrderNonLinearDS");

  _nsdsSMC->insertDynamicalSystem(_DS_SMC);
  _nsdsSMC->setName(_DS_SMC, "plant_SMC");
  _nsdsSMC->link(_interactionSMC, _DS_SMC);
  _nsdsSMC->setControlProperty(_interactionSMC, true);
  _nsdsSMC->topology()->setName(_interactionSMC, "Sgn_SMC");

  _simulationSMC->setName("linear sliding mode controller simulation");
  _simulationSMC->insertIntegrator(_integratorSMC);
  _simulationSMC->associate(_integratorSMC, _DS_SMC);

  // OneStepNsProblem
  _OSNSPB_SMC->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = _precision;
  _simulationSMC->insertNonSmoothProblem(_OSNSPB_SMC);
  // Finally we can initialize everything ...
  _simulationSMC->associate(_integratorSMC, _DS_SMC);

  // _SMC->setSimulation(_simulationSMC);
  _simulationSMC->initialize();

  // Handy
  _eventsManager = _simulationSMC->eventsManager();
  _lambda = std::make_shared<siconos::algebra::SiconosVector>(sDim);
  _lambda = _interactionSMC->lambda(0);
  _us = std::make_shared<siconos::algebra::SiconosVector>(sDim);
  _us->setZero();
  _ueq = std::make_shared<siconos::algebra::SiconosVector>(sDim);
  _ueq->setZero();

  if (_Csurface) {
    _invCB = std::make_shared<siconos::algebra::SiconosMatrix>(_Csurface->rows(), _B->cols());
    *_invCB = (*_Csurface * *_B).inverse();
  }
  DEBUG_END(
      "siconos::control::CommonSMC::initialize(const "
      "siconos::modeling::NonSmoothDynamicalSystem & nsds, const "
      "Simulation & s)\n");
}

void siconos::control::CommonSMC::computeUeq() {
  DEBUG_BEGIN("void siconos::control::CommonSMC::computeUeq()\n");
  assert(std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC) &&
         "siconos::control::CommonSMC::computeUeq the DS should be linear");
  assert(_Csurface &&
         "siconos::control::CommonSMC::computeUeq the sliding variable should be linear "
         "subpsace of the state");

  auto LinearDS_SMC = std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC);

  if (LinearDS_SMC->hasA()) {
    auto n = LinearDS_SMC->dimension();
    auto xTk = std::make_shared<siconos::algebra::SiconosVector>(_sensor->y());

    auto zoh =
        std::static_pointer_cast<siconos::integrators::ZeroOrderHoldOSI>(_integratorSMC);

    // compute (CB)^{-1}CA
    siconos::algebra::SiconosDenseMatrix quasiProjB_A =
        *_invCB * *_Csurface * LinearDS_SMC->A();
    *_ueq = (_thetaSMC - 1) * quasiProjB_A * *xTk;

    // equivalent part, explicit contribution
    // tmpN = B^{*}(CB)^{-1}CA
    siconos::algebra::SiconosDenseMatrix tmpN = zoh->Bd(_DS_SMC) * quasiProjB_A;

    // W = I + \theta B^{*})CB)^{-1}CA

    siconos::algebra::SiconosMatrix tmpW =
        siconos::algebra::SiconosMatrix::Identity(n, n) + _thetaSMC * tmpN;
    siconos::algebra::SiconosDenseLUMatrix luW(tmpW);

    // compute e^{Ah}x_k
    *xTk = zoh->Ad(_DS_SMC) * *xTk;

    // xTk = (e^{Ah}-(1-\theta)\Psi_k\Pi_B A)x_k
    *xTk += (_thetaSMC - 1) * tmpN * _sensor->y();
    // compute the solution x_{k+1} of the system W*x_{k+1} = x_k
    *xTk = luW.solve(*xTk);

    // add the contribution from the implicit part to ueq
    *_ueq -= _thetaSMC * quasiProjB_A * *xTk;

  } else {
    _ueq->setZero();
  }
  DEBUG_END("void siconos::control::CommonSMC::computeUeq()\n");
}

void siconos::control::CommonSMC::setCsurface(
    std::shared_ptr<siconos::algebra::SiconosMatrix> newC) {
  // check dimensions ...
  _Csurface = newC;
}

void siconos::control::CommonSMC::setSaturationMatrix(
    std::shared_ptr<siconos::algebra::SiconosMatrix> newSat) {
  // check dimensions ...
  if (newSat->cols() != _B->cols()) {
    THROW_EXCEPTION(
        "siconos::control::CommonSMC::setSaturationMatrixPtr - inconstency between the "
        "dimension of the state "
        "space and D");
  } else {
    _D = newSat;
  }
}

void siconos::control::CommonSMC::setTimeDiscretisation(
    const siconos::simulation::TimeDiscretisation& td) {
  _td = std::make_shared<siconos::simulation::TimeDiscretisation>(td);
}

void siconos::control::CommonSMC::sete(
    const siconos::modeling::func_prototypes::FunctionS_V& fct) {
  computee_ = fct;
  isRelationLinear_ = true;
}

void siconos::control::CommonSMC::setComputehFunction(

    const siconos::modeling::func_prototypes::FunctionBVSV_V& fct) {
  computeh_ = fct;
  isRelationLinear_ = false;
}

void siconos::control::CommonSMC::setComputeJacobianhOver_stateFunction(

    const siconos::modeling::func_prototypes::FunctionBVSV_M& fct) {
  computejacobianhOver_state_ = fct;
  isRelationLinear_ = false;
}

void siconos::control::CommonSMC::setComputeJacobianhOver_lambdaFunction(

    const siconos::modeling::func_prototypes::FunctionBVSV_M& fct) {
  computejacobianhOver_lambda_ = fct;
  isRelationLinear_ = false;
}

void siconos::control::CommonSMC::setComputeJacobiangOver_lambdaFunction(

    const siconos::modeling::func_prototypes::FunctionBVSV_M& fct) {
  computejacobiangOver_lambda_ = fct;
  isRelationLinear_ = false;
}
