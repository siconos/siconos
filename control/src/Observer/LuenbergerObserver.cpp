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

#include "LuenbergerObserver.hpp"

#include "ControlSensor.hpp"
#include "ControlZOHAdditionalTerms.hpp"
#include "FirstOrderLinearDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "TimeStepping.hpp"
#include "Topology.hpp"
#include "ZeroOrderHoldOSI.hpp"
// #define DEBUG_BEGIN_END_ONLY
//  #define DEBUG_NOCOLOR
//  #define DEBUG_STDOUT
//  #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::control::LuenbergerObserver::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s) {
  DEBUG_BEGIN(
      "void siconos::control::LuenbergerObserver::initialize(const "
      "siconos::modeling::NonSmoothDynamicalSystem& nsds, const Simulation &s)\n");
  if (!_C) {
    THROW_EXCEPTION(
        "siconos::control::LuenbergerObserver::initialize - you have to set C before "
        "initializing the Observer");
  } else {
    Observer::initialize(nsds, s);
  }
  bool isDSinDSG0 = true;
  auto& originalDSG0 = *nsds.topology()->dSG(0);
  siconos::graphs::DynamicalSystemsGraph::VDescriptor originaldsgVD;
  if (!_DS)  // No DynamicalSystem was given
  {
    // We can only work with FirstOrderNonLinearDS and FirstOrderLinearDS
    // We can use the Visitor mighty power to check if we have the right type
    auto observedDS = _sensor->getDS();
    // create the DS for the controller
    // if the DS we use is different from the DS we are controlling
    // when we want for instant to see how well the controller behaves
    // if the plant model is not exact, we can use the setSimulatedDS
    // method
    if (auto folds =
            std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(observedDS)) {
      _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*folds);
    } else
      THROW_EXCEPTION("LuenbergerObserver is only implemented for FirstOrderLinearDS");
    // is it controlled ?
    originaldsgVD = originalDSG0.descriptor(_sensor->getDS());
  } else {
    // is it controlled ?
    if (originalDSG0.is_vertex(_DS))
      originaldsgVD = originalDSG0.descriptor(_DS);
    else
      isDSinDSG0 = false;
  }

  // Initialize with the guessed state
  _DS->setX0(*_xHat);  // Shared memory view!
  _DS->resetToInitialState();
  DEBUG_EXPR(_DS->display(););
  _e = std::make_shared<siconos::algebra::SiconosVector>(_C->rows());
  _y = std::make_shared<siconos::algebra::SiconosVector>(_C->rows());

  auto t0 = nsds.t0();
  auto h = s.currentTimeStep();
  auto T = nsds.finalT() + h;
  _nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
  _integrator = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();

  std::static_pointer_cast<siconos::integrators::ZeroOrderHoldOSI>(_integrator)
      ->setExtraAdditionalTerms(std::make_shared<ControlZOHAdditionalTerms>());

  _nsds->insertDynamicalSystem(_DS);

  // Add the necessary properties
  auto& DSG0 = *_nsds->topology()->dSG(0);
  auto dsgVD = DSG0.descriptor(_DS);
  // Observer part
  DSG0.L[dsgVD] = _L;
  DSG0.e[dsgVD] = _e;

  // Was the original DynamicalSystem controlled ?
  if (isDSinDSG0 && originalDSG0.B.hasKey(originaldsgVD)) {
    DSG0.B[dsgVD] = originalDSG0.B[originaldsgVD];
    assert(
        originalDSG0.u[originaldsgVD] &&
        "A DynamicalSystem is controlled but its control input has not been initialized yet");
    DSG0.u[dsgVD] = originalDSG0.u[originaldsgVD];
  }

  // all necessary things for simulation
  _simulation = std::make_shared<siconos::simulation::TimeStepping>(_nsds, _td, 0);
  _simulation->associate(_integrator, _DS);

  // initialize error
  *_y = _sensor->y();
  DEBUG_END(
      "void siconos::control::LuenbergerObserver::initialize(const "
      "siconos::modeling::NonSmoothDynamicalSystem& nsds, const Simulation &s)\n");
}

void siconos::control::LuenbergerObserver::process() {
  if (!_pass)
    _pass = true;
  else {
    *_e = *_C * *_xHat - *_y;

    // get measurement from sensor
    const auto& y = _sensor->y();

    // TODO theta method on the error
    _simulation->computeOneStep();
    //    _simulation->nextStep();

    // update the current measured value
    *_y = y;
    *_xHat = _DS->x_read();  // Copy
  }
}
