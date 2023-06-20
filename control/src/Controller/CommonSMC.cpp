/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "ControlSensor.hpp"
#include "EulerMoreauOSI.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderNonLinearR.hpp"
#include "FirstOrderR_helpers.hpp"
#include "FirstOrderType1R.hpp"
#include "FirstOrderType2R.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NumericsSolversNamespace.h"
#include "Relay.hpp"
#include "RelayNSL.hpp"
#include "SiconosAlgebraTools.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "Topology.hpp"
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
  if (!_Csurface && _pluginhName.empty()) {
    THROW_EXCEPTION(
        "siconos::control::CommonSMC::initialize - you have to set either _Csurface or h(.) "
        "before initializing "
        "the Actuator");
  } else {
    if (_Csurface && !_u)
      _u = std::make_shared<siconos::algebra::SiconosVector>(_Csurface->size(0), 0);

    Actuator::initialize(nsds, s);
  }
  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // We can use the Visitor mighty power to check if we have the right type
  auto DS = _sensor->getDS();
  // create the DS for the controller
  // if the DS we use is different from the DS we are controlling
  // when we want for instant to see how well the controller behaves
  // if the plant model is not exact, we can use the setSimulatedDS
  // method
  if (auto folds = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearTIDS>(DS)) {
    _DS_SMC = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(*folds);
    // We have to reset the _pluginb
    auto dummyb = std::make_shared<siconos::algebra::SiconosVector>(_DS_SMC->n(), 0);
    std::static_pointer_cast<siconos::modeling::FirstOrderLinearTIDS>(_DS_SMC)->setbPtr(
        dummyb);
  } else if (auto folds =
                 std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(DS)) {
    _DS_SMC = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*folds);
    std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC)
        ->setComputebFunction(nullptr);
    // We have to reset the _pluginb
    auto dummyb = std::make_shared<siconos::algebra::SiconosVector>(_DS_SMC->n(), 0);
    std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC)->setbPtr(dummyb);
  } else if (auto folds =
                 std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(DS)) {
    _DS_SMC = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(*folds);
  } else {
    THROW_EXCEPTION("LinearSMC is only  implemented for FirstOrderNonLinearDS");
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
  if (!_plugingName.empty()) {
    if (_pluginhName.empty())
      THROW_EXCEPTION(
          "LinearSMC::initialize - the Controller has a function g set but _pluginhName is "
          "not set\n You must supply a function to compute y=h(x,...)");
    if (!_pluginJacgxName.empty())  // Is the relation the most generic NL one ?
    {
      DEBUG_PRINT("A FirstOrderNonLinearR is created for the _relationSMC\n");
      _relationSMC = std::make_shared<siconos::modeling::FirstOrderNonLinearR>();
      _relationSMC->setComputeJacgxFunction(
          siconos::plugins::getPluginName(_pluginJacgxName),
          siconos::plugins::getPluginFunctionName(_pluginJacgxName));

      siconos::modeling::FirstOrderRHelpers::JacglambdaSetter(*_relationSMC, _B,
                                                              _pluginJacglambdaName);
      siconos::modeling::FirstOrderRHelpers::JachxSetter(*_relationSMC, _Csurface,
                                                         _pluginJachxName);
      if (_pluginJachlambdaName.empty() && !_D)
        _D = std::make_shared<siconos::algebra::SimpleMatrix>(sDim, sDim, 0);
      siconos::modeling::FirstOrderRHelpers::JachlambdaSetter(*_relationSMC, _D,
                                                              _pluginJachlambdaName);
    } else if (!_pluginJachlambdaName.empty() || _D)  // Type2R ?
    {
      DEBUG_PRINT("A FirstOrderType2R is created for the _relationSMC\n");
      _relationSMC = std::make_shared<siconos::modeling::FirstOrderType2R>();
      siconos::modeling::FirstOrderRHelpers::JacglambdaSetter(*_relationSMC, _B,
                                                              _pluginJacglambdaName);
      siconos::modeling::FirstOrderRHelpers::JachxSetter(*_relationSMC, _Csurface,
                                                         _pluginJachxName);
      siconos::modeling::FirstOrderRHelpers::JachlambdaSetter(*_relationSMC, _D,
                                                              _pluginJachlambdaName);
    } else  // Type1R
    {
      DEBUG_PRINT("A FirstOrderType1R is created for the _relationSMC\n");
      _relationSMC = std::make_shared<siconos::modeling::FirstOrderType1R>();
      siconos::modeling::FirstOrderRHelpers::JacglambdaSetter(*_relationSMC, _B,
                                                              _pluginJacglambdaName);
      siconos::modeling::FirstOrderRHelpers::JachxSetter(*_relationSMC, _Csurface,
                                                         _pluginJachxName);
    }

    _relationSMC->setComputehFunction(siconos::plugins::getPluginName(_pluginhName),
                                      siconos::plugins::getPluginFunctionName(_pluginhName));
    _relationSMC->setComputegFunction(siconos::plugins::getPluginName(_plugingName),
                                      siconos::plugins::getPluginFunctionName(_plugingName));

    if (_computeResidus) {
      _simulationSMC->setComputeResiduY(true);
      _simulationSMC->setComputeResiduR(true);
      _simulationSMC->setUseRelativeConvergenceCriteron(false);
    }
  } else {
    if (!_plugineName.empty()) {
      DEBUG_PRINT("A FirstOrderLinearR is created for the _relationSMC\n");
      _relationSMC = std::make_shared<siconos::modeling::FirstOrderLinearR>(_Csurface, _B);
      _relationSMC->setComputeEFunction(siconos::plugins::getPluginName(_plugineName),
                                        siconos::plugins::getPluginFunctionName(_plugineName));
    } else {
      DEBUG_PRINT("A FirstOrderLinearTIR is created for the _relationSMC\n");
      _relationSMC = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(_Csurface, _B);
    }
    std::static_pointer_cast<siconos::modeling::FirstOrderLinearTIR>(_relationSMC)
        ->setDPtr(_D);
  }

  // _nsLawSMC and the OSNSP can be defined in derived classes, like twisting
  if (!_nsLawSMC)
    _nsLawSMC = std::make_shared<siconos::modeling::RelayNSL>(sDim, -_alpha, _alpha);
  if (!_OSNSPB_SMC)
    _OSNSPB_SMC = std::make_shared<siconos::nonsmooth_formulations::Relay>(_numericsSolverId);

  _interactionSMC = std::make_shared<siconos::modeling::Interaction>(_nsLawSMC, _relationSMC);

  if (std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearTIDS>(DS) || std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(DS))
    _integratorSMC = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();
  else if(std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(DS))
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
  _ueq = std::make_shared<siconos::algebra::SiconosVector>(sDim);

  if (_Csurface) {
    auto tmpM =
        std::make_shared<siconos::algebra::SimpleMatrix>(_Csurface->size(0), _B->size(1));
    _invCB = std::make_shared<siconos::algebra::SimpleMatrix>(*tmpM);
    siconos::algebra::prod(*_Csurface, *_B, *tmpM);
    siconos::algebra::InvertMatrix(*tmpM->dense(), *_invCB->dense());
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
  auto& LinearDS_SMC =
      *std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS_SMC);
  auto n = LinearDS_SMC.A()->size(1);
  // equivalent part, explicit contribution
  auto tmpM1 = std::make_shared<siconos::algebra::SimpleMatrix>(_Csurface->size(0), n);
  auto tmpN = std::make_shared<siconos::algebra::SimpleMatrix>(n, n);
  auto quasiProjB_A = std::make_shared<siconos::algebra::SimpleMatrix>(_invCB->size(0), n);
  auto tmpW = std::make_shared<siconos::algebra::SimpleMatrix>(n, n, 0);
  auto xTk = std::make_shared<siconos::algebra::SiconosVector>(_sensor->y());
  tmpW->eye();
  siconos::algebra::prod(*_Csurface, *LinearDS_SMC.A(), *tmpM1);
  // compute (CB)^{-1}CA
  siconos::algebra::prod(*_invCB, *tmpM1, *quasiProjB_A);
  siconos::algebra::prod(_thetaSMC - 1, *quasiProjB_A, *xTk, *_ueq);

  // equivalent part, implicit contribution
  // XXX when to call this ?
  auto& zoh =
      *std::static_pointer_cast<siconos::integrators::ZeroOrderHoldOSI>(_integratorSMC);
  zoh.updateMatrices(_DS_SMC);

  // tmpN = B^{*}(CB)^{-1}CA
  siconos::algebra::prod(zoh.Bd(_DS_SMC), *quasiProjB_A, *tmpN, true);
  // W = I + \theta B^{*})CB)^{-1}CA
  siconos::algebra::scal(_thetaSMC, *tmpN, *tmpW, false);
  // compute e^{Ah}x_k
  siconos::algebra::prod(zoh.Ad(_DS_SMC), *xTk, *xTk);
  // xTk = (e^{Ah}-(1-\theta)\Psi_k\Pi_B A)x_k
  siconos::algebra::prod(_thetaSMC - 1, *tmpN, _sensor->y(), *xTk, false);
  // compute the solution x_{k+1} of the system W*x_{k+1} = x_k
  tmpW->Solve(*xTk);
  // add the contribution from the implicit part to ueq
  siconos::algebra::prod(-_thetaSMC, *quasiProjB_A, *xTk, *_ueq, false);
  DEBUG_END("void siconos::control::CommonSMC::computeUeq()\n");
}

void siconos::control::CommonSMC::setCsurface(
    std::shared_ptr<siconos::algebra::SimpleMatrix> newC) {
  // check dimensions ...
  _Csurface = newC;
}

void siconos::control::CommonSMC::setSaturationMatrix(
    std::shared_ptr<siconos::algebra::SimpleMatrix> newSat) {
  // check dimensions ...
  if (newSat->size(1) != _B->size(1)) {
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

void siconos::control::CommonSMC::sete(const std::string& plugin) { _plugineName = plugin; }

void siconos::control::CommonSMC::seth(const std::string& plugin) { _pluginhName = plugin; }
void siconos::control::CommonSMC::setJachx(const std::string& plugin) {
  _pluginJachxName = plugin;
}
void siconos::control::CommonSMC::setJachlambda(const std::string& plugin) {
  _pluginJachlambdaName = plugin;
}
void siconos::control::CommonSMC::setg(const std::string& plugin) { _plugingName = plugin; }
void siconos::control::CommonSMC::setJacgx(const std::string& plugin) {
  _pluginJacgxName = plugin;
}
void siconos::control::CommonSMC::setJacglambda(const std::string& plugin) {
  _pluginJacglambdaName = plugin;
}
