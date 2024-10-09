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

#include "FirstOrderNonLinearDS.hpp"

#include "PluggedObject.hpp"  // for getPluginName ...
#include "PluginTypes.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include <iostream>
#include <memory>

#include "siconos_debug.h"

// From a minimum set of data
siconos::modeling::FirstOrderNonLinearDS::FirstOrderNonLinearDS(
    std::shared_ptr<siconos::algebra::SiconosVector> initial_state)
    : DynamicalSystem(initial_state->size()) {
  DEBUG_PRINT("internal _init from FirstOrderNonLinearDS\n");

  // Memory allocation only for required parts of the DS:
  // state (initial and current). All other operators are optional and
  // allocated with 'set'-like methods.
  assert(_n > 0 && "dynamical system dimension should be greater than 0.");
  // Set initial conditions
  _x0 = initial_state;

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.
  // _x.resize(2); done in base class constructor.
  _x[0] = std::make_shared<siconos::algebra::SiconosVector>(*_x0);
  _x[1] = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _x[1]->setZero();
  _r = std::make_shared<siconos::algebra::SiconosVector>(
      _n);  // FP: move this to initializeNonSmoothInput?
  _r->setZero();
  _zeroPlugin();
  // dot x = r
}

// From a minimum set of data
siconos::modeling::FirstOrderNonLinearDS::FirstOrderNonLinearDS(
    std::shared_ptr<siconos::algebra::SiconosVector> initial_state, const std::string &fPlugin,
    const std::string &jacobianfxPlugin)
    : FirstOrderNonLinearDS(initial_state) {
  // == f and its jacobian ==
  // Allocation and link with the plug-in
  setComputeFFunction(siconos::plugins::getPluginName(fPlugin),
                      siconos::plugins::getPluginFunctionName(fPlugin));
  setComputeJacobianfxFunction(siconos::plugins::getPluginName(jacobianfxPlugin),
                               siconos::plugins::getPluginFunctionName(jacobianfxPlugin));
  // dot x  = f(x, z , t) + r
}

// Copy constructor
siconos::modeling::FirstOrderNonLinearDS::FirstOrderNonLinearDS(
    const FirstOrderNonLinearDS &FONLDS)
    : DynamicalSystem(FONLDS) {
  if (FONLDS.M()) _M = std::make_shared<siconos::algebra::SiconosMatrix>(*(FONLDS.M()));
  if (FONLDS.f()) _f = std::make_shared<siconos::algebra::SiconosVector>(*(FONLDS.f()));
  if (FONLDS.jacobianfx())
    _jacobianfx = std::make_shared<siconos::algebra::SiconosMatrix>(*(FONLDS.jacobianfx()));
  if (FONLDS.getPluginF())
    _pluginf = std::make_shared<siconos::plugins::PluggedObject>(*(FONLDS.getPluginF()));
  if (FONLDS.getPluginJacxf())
    _pluginJacxf =
        std::make_shared<siconos::plugins::PluggedObject>(*(FONLDS.getPluginJacxf()));
  if (FONLDS.getPluginM())
    _pluginM = std::make_shared<siconos::plugins::PluggedObject>(*(FONLDS.getPluginM()));
  if (FONLDS.LU_M()) {
    LU_M_ = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*_M);
    hasLU_M_ = true;
  }
  // Memory stuff to me moved to graph/osi
  if (FONLDS.fold())
    _fold = std::make_shared<siconos::algebra::SiconosVector>(*(FONLDS.fold()));
  _rMemory = FONLDS.rMemory();
}

void siconos::modeling::FirstOrderNonLinearDS::_zeroPlugin() {
  _pluginf = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJacxf = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginM = std::make_shared<siconos::plugins::PluggedObject>();
}

void siconos::modeling::FirstOrderNonLinearDS::initRhs(double time) {
  computeRhs(time);

  // !! jacxRhs must always be allocated (we must check this?)!!
  if (!_jacxRhs)  // if not allocated with a set or anything else
  {
    if (_jacobianfx && !_M)  // if M is not defined, then jacobianfx = jacobianRhsx, no memory
                             // allocation for that one.
      _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(_jacobianfx);
    else {  //  if (_jacobianfx && _M) or if(!jacobianRhsx)
      auto tmp = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
      tmp->setZero();
      _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(*tmp);
    }
    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianRhsx(time);
}

void siconos::modeling::FirstOrderNonLinearDS::updatePlugins(double time) {
  if (_M) computeM(time);
  if (_f) {
    computef(time, _x[0]);
    computeJacobianfx(time, _x[0]);
  }
}

void siconos::modeling::FirstOrderNonLinearDS::initializeNonSmoothInput(unsigned int level) {
  /**\warning V.A. _r should be initialized here and not in  the constructor
   * The level should also be used if we need more that one _r
   */
  if (!_r) _r = std::make_shared<siconos::algebra::SiconosVector>(_n);
}

void siconos::modeling::FirstOrderNonLinearDS::resetToInitialState() { *(_x[0]) = *_x0; }

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void siconos::modeling::FirstOrderNonLinearDS::initMemory(unsigned int steps) {
  DynamicalSystem::initMemory(steps);

  if (_f && !_fold) _fold = std::make_shared<siconos::algebra::SiconosVector>(_n);

  if (steps == 0)
    std::cout << "Warning : siconos::modeling::FirstOrderNonLinearDS::initMemory with size "
                 "equal to zero"
              << std::endl;
  else
    _rMemory.setMemorySize(steps, _n);
}

void siconos::modeling::FirstOrderNonLinearDS::swapInMemory() {
  DEBUG_BEGIN("void siconos::modeling::FirstOrderNonLinearDS::swapInMemory()\n");
  _xMemory.swap(*_x[0]);
  _rMemory.swap(*_r);
  if (_f) {
    assert(_fold);
    *_fold = *_f;
  }
  DEBUG_EXPR(_xMemory.display());
  DEBUG_END("void siconos::modeling::FirstOrderNonLinearDS::swapInMemory()\n");
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void siconos::modeling::FirstOrderNonLinearDS::setComputeMFunction(
    const std::string &pluginPath, const std::string &functionName) {
  if (!_M) _M = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);

  _pluginM->setComputeFunction(pluginPath, functionName);
}

void siconos::modeling::FirstOrderNonLinearDS::setComputeMFunction(
    siconos::plugins::FPtr1 fct) {
  if (!_M) _M = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);

  _pluginM->setComputeFunction((void *)fct);
}

void siconos::modeling::FirstOrderNonLinearDS::setComputeFFunction(
    const std::string &pluginPath, const std::string &functionName) {
  if (!_f) _f = std::make_shared<siconos::algebra::SiconosVector>(_n);

  _pluginf->setComputeFunction(pluginPath, functionName);
}

void siconos::modeling::FirstOrderNonLinearDS::setComputeFFunction(
    siconos::plugins::FPtr1 fct) {
  if (!_f) _f = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _pluginf->setComputeFunction((void *)fct);
}

void siconos::modeling::FirstOrderNonLinearDS::setComputeJacobianfxFunction(
    const std::string &pluginPath, const std::string &functionName) {
  if (!_jacobianfx) _jacobianfx = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _pluginJacxf->setComputeFunction(pluginPath, functionName);
}

void siconos::modeling::FirstOrderNonLinearDS::setComputeJacobianfxFunction(
    siconos::plugins::FPtr1 fct) {
  if (!_jacobianfx) _jacobianfx = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _pluginJacxf->setComputeFunction((void *)fct);
}

void siconos::modeling::FirstOrderNonLinearDS::computeM(double time) {
  if (_pluginM->fPtr && _M) {
    ((FNLDSPtrfct)_pluginM->fPtr)(time, _n, &((*(_x[0]))(0)), &(*_M)(0, 0), _z->size(),
                                  &(*_z)(0));
  }
}

void siconos::modeling::FirstOrderNonLinearDS::computef(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> state) {
  if (_f && _pluginf->fPtr)
    ((FNLDSPtrfct)_pluginf->fPtr)(time, _n, &(*state)(0), &(*_f)(0), _z->size(), &(*_z)(0));
}

void siconos::modeling::FirstOrderNonLinearDS::computeJacobianfx(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> state) {
  if (_jacobianfx && _pluginJacxf->fPtr)
    ((FNLDSPtrfct)_pluginJacxf->fPtr)(time, _n, state->data(), &(*_jacobianfx)(0, 0),
                                      _z->size(), _z->data());
}

void siconos::modeling::FirstOrderNonLinearDS::computeRhs(double time) {
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *_x[1] = *_r;  // Warning: p update is done in Interactions/Relations

  if (_f) {
    computef(time, _x[0]);
    *(_x[1]) += *_f;
  }

  if (_M) {
    computeM(time);
    // allocate invM at the first call of the present function
    LU_M_ = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*_M);
    *(_x[1]) = LU_M_->solve(*(_x[1]));
  }
}

void siconos::modeling::FirstOrderNonLinearDS::computeJacobianRhsx(double time) {
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianfx + jacobianX(T.u))
  // At the time, second term is set to zero.
  // assert(!_pluginJacxf->fPtr &&
  // "siconos::modeling::FirstOrderNonLinearDS::computeJacobianRhsx: there is no plugin to
  // compute the jacobian of f");

  computeJacobianfx(time, _x[0]);
  // solve M*jacobianXRhS = jacobianfx
  if (_M && _jacobianfx) {
    (_jacxRhs->copyBlock(0, 0, _jacobianfx));
    // copy _M into LU_M_ for LU-factorisation, at the first call of this function.

    computeM(time);
    LU_M_ = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*_M);
    *(_jacxRhs->block(0, 0)) = LU_M_->solve(*(_jacxRhs->block(0, 0)));
  }
  // else jacobianRhsx = jacobianfx, pointers equality set in initRhs
}

// ===== MISCELLANEOUS ====

void siconos::modeling::FirstOrderNonLinearDS::display(bool brief) const {
  std::cout << " =====> First Order Non Linear DS (number: " << _number << ")." << std::endl;
  std::cout << "- n (size) : " << _n << std::endl;
  std::cout << "- x " << std::endl;
  if (_x[0])
    _x[0]->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- x0 " << std::endl;
  if (_x0)
    _x0->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- M: " << std::endl;
  if (_M)
    _M->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << " ============================================" << std::endl;
}

void siconos::modeling::FirstOrderNonLinearDS::resetAllNonSmoothParts() { _r->setZero(); }

void siconos::modeling::FirstOrderNonLinearDS::resetNonSmoothPart(unsigned int level) {
  // V.A. 28/05/2012:  for the moment various level are not used for First Order systems
  // assert(0);
  _r->setZero();
}

void siconos::modeling::FirstOrderNonLinearDS::acceptSP(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
