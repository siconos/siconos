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

#include "MatrixIntegrator.hpp"

#include "EventDriven.hpp"
#include "FirstOrderLinearDS.hpp"
#include "LsodarOSI.hpp"
#include "SiconosConst.hpp"  // MACHINE_PREC
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "SubPluggedObject.hpp"
#include "TimeDiscretisation.hpp"
// #define DEBUG_WHERE_MESSAGES

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

siconos::simulation::MatrixIntegrator::MatrixIntegrator(
    const siconos::modeling::DynamicalSystem& ds,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    std::shared_ptr<TimeDiscretisation> td,
    const Eigen::Ref<const siconos::algebra::SiconosMatrix>& E) {
  // Copy td
  auto tmp = std::make_shared<siconos::simulation::TimeDiscretisation>(*td);
  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(*td);

  DEBUG_EXPR(ds.display(););

  if (auto folds = dynamic_cast<const siconos::modeling::FirstOrderLinearDS*>(&ds)) {
    _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*folds);
    _isConst = (_TD->hConst()) && !(folds->isTimeInvariant()) ? true : false;
  }

  _DS->setNumber(9999999);
  DEBUG_EXPR(_DS->display(););
  // integration stuff
  _nsds =
      std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(nsds.t0(), nsds.finalT());

  _OSI = std::make_shared<siconos::integrators::LsodarOSI>();
  _nsds->insertDynamicalSystem(_DS);
  _sim = std::make_shared<siconos::simulation::EventDriven>(_nsds, _TD, 0);
  _sim->associate(_OSI, _DS);
  _sim->setName("Matrix integrator simulation");
  // change tolerance
  _OSI->setTol(1, 10 * siconos::internal::MACHINE_PREC, 5 * siconos::internal::MACHINE_PREC);

  _E = std::make_shared<siconos::algebra::ConstMapType>(E.data(), E.rows(), E.cols());
  _mat = std::make_shared<siconos::algebra::SiconosMatrix>(E);
  _mat->setZero();
}

// --> FP: TO BE REVIEWED.
// siconos::simulation::MatrixIntegrator::MatrixIntegrator(
//     const siconos::modeling::DynamicalSystem& ds,
//     const siconos::modeling::NonSmoothDynamicalSystem& nsds,
//     std::shared_ptr<TimeDiscretisation> td,
//     std::shared_ptr<siconos::plugins::PluggedObject> plugin, const unsigned int p)
//     : MatrixIntegrator{ds, nsds, td, nullptr} {
//   _plugin = plugin;
//   _isConst = false;
//   auto n = ds.dimension();
//   _mat = std::make_shared<siconos::algebra::SiconosMatrix>(n, p);
//   _mat->setZero();
//   _spo = std::make_shared<siconos::plugins::SubPluggedObject>(*_plugin, n, p);
//   std::static_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS)->setPluginB(_spo);
// }

siconos::simulation::MatrixIntegrator::MatrixIntegrator(
    const siconos::modeling::DynamicalSystem& ds,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    std::shared_ptr<TimeDiscretisation> td) {
  // Copy td
  auto tmp = std::make_shared<siconos::simulation::TimeDiscretisation>(*td);
  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(*td);

  DEBUG_EXPR(ds.display(););

  if (auto folds = dynamic_cast<const siconos::modeling::FirstOrderLinearDS*>(&ds)) {
    _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*folds);
    _isConst = (_TD->hConst()) && !(folds->isTimeInvariant()) ? true : false;
  }

  _DS->setNumber(9999999);
  DEBUG_EXPR(_DS->display(););
  // integration stuff
  _nsds =
      std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(nsds.t0(), nsds.finalT());

  _OSI = std::make_shared<siconos::integrators::LsodarOSI>();
  _nsds->insertDynamicalSystem(_DS);
  _sim = std::make_shared<siconos::simulation::EventDriven>(_nsds, _TD, 0);
  _sim->associate(_OSI, _DS);
  _sim->setName("Matrix integrator simulation");
  // change tolerance
  _OSI->setTol(1, 10 * siconos::internal::MACHINE_PREC, 5 * siconos::internal::MACHINE_PREC);

  auto n = ds.dimension();
  _mat = std::make_shared<siconos::algebra::SiconosMatrix>(n, n);
  _mat->setZero();
}

void siconos::simulation::MatrixIntegrator::integrate() {
  DEBUG_BEGIN("siconos::simulation::MatrixIntegrator::integrate()\n");

  auto x = *_DS->x();

  auto fods = std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(_DS);
  assert(fods && "MatrixIntegrator only available for first order DS.");

  auto p = _mat->cols();
  auto x0 = fods->x0_ptr();  // Shared memory !!

  std::shared_ptr<siconos::algebra::SiconosVector> Ecol;
  if (auto linear_ds =
          std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(fods)) {
    for (unsigned int i = 0; i < p; i++) {
      x0->setZero();
      if (_E) {
        siconos::algebra::SiconosVector b = _E->col(i);
        linear_ds->setConstantbVector(b);
      } else if (_plugin)   // ?????
        _spo->setIndex(i);  /// ?????
      else
        (*x0)(i) = 1.;

      // Reset LsodarOSI
      //_OSI->setIsInitialized(false);
      _DS->resetToInitialState();
      _sim->setIstate(1);
      _sim->advanceToEvent();
      _mat->col(i) = x;
    }

  } else  // Non linear case
  {
    for (unsigned int i = 0; i < p; i++) {
      x0->setZero();
      if (_plugin)
        _spo->setIndex(i);
      else
        (*x0)(i) = 1.;

      // Reset LsodarOSI
      //_OSI->setIsInitialized(false);
      _DS->resetToInitialState();
      _sim->setIstate(1);
      _sim->advanceToEvent();
      _mat->col(i) = x;
    }
  }

  _sim->processEvents();
  //_DS->resetToInitialState();

  DEBUG_EXPR(_mat->display(););
  DEBUG_EXPR(_DS->display(););
  DEBUG_END("siconos::simulation::MatrixIntegrator::integrate()\n");
}
