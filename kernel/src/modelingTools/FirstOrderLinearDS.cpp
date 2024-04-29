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
#include "FirstOrderLinearDS.hpp"

#include "SiconosException.hpp"
#include "SiconosMatrixVectorOp.hpp"  // For prod matrix-vectors
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include <iostream>

#include "siconos_debug.h"

// --- Constructors ---
// From a minimum set of data, A and b connected to a plug-in
siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    std::shared_ptr<siconos::algebra::SiconosVector> newX0,
    const std::string &APlugin, const std::string &bPlugin)
    : FirstOrderNonLinearDS(newX0) {
  _zeroPlugin();
  setComputeAFunction(siconos::plugins::getPluginName(APlugin),
                      siconos::plugins::getPluginFunctionName(APlugin));
  setComputebFunction(siconos::plugins::getPluginName(bPlugin),
                      siconos::plugins::getPluginFunctionName(bPlugin));
  // dot x = A(t)x + b(t)
}

// From a minimum set of data, A from a given matrix
siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    std::shared_ptr<siconos::algebra::SiconosVector> newX0,
    std::shared_ptr<siconos::algebra::SiconosMatrix> newA)
    : FirstOrderNonLinearDS(newX0), _hasConstantA(true) {
  _zeroPlugin();
  if ((newA->size(0) != _n) || (newA->size(1) != _n))
    THROW_EXCEPTION(
        "FirstOrderLinearDS - constructor(number,x0,A): inconsistent "
        "dimensions "
        "with problem size for input matrix A");
  _A = newA;
}

siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    std::shared_ptr<siconos::algebra::SiconosVector> newX0)
    : FirstOrderNonLinearDS(newX0) {
  _zeroPlugin();
}

// From a minimum set of data, A from a given matrix
siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    std::shared_ptr<siconos::algebra::SiconosVector> newX0,
    std::shared_ptr<siconos::algebra::SiconosMatrix> newA,
    std::shared_ptr<siconos::algebra::SiconosVector> newB)
    : FirstOrderNonLinearDS(newX0), _hasConstantA(true), _hasConstantB(true) {
  _zeroPlugin();

  if ((newA->size(0) != _n) || (newA->size(1) != _n))
    THROW_EXCEPTION(
        "FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions "
        "with "
        "problem size for input matrix A");

  if (newB->size() != _n)
    THROW_EXCEPTION(
        "FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions "
        "with "
        "problem size for input vector b ");

  _A = newA;
  _b = newB;
}

// Copy constructor
siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    const FirstOrderLinearDS &FOLDS)
    : FirstOrderNonLinearDS(FOLDS) {
  _zeroPlugin();
  if (FOLDS.A())
    _A = std::make_shared<siconos::algebra::SiconosMatrix>(*(FOLDS.A()));
  if (FOLDS.b())
    _b = std::make_shared<siconos::algebra::SiconosVector>(*(FOLDS.b()));

  _hasConstantA = FOLDS.hasConstantA();
  _hasConstantB = FOLDS.hasConstantB();

  if (not _hasConstantA)
    _pluginA = std::make_shared<siconos::plugins::PluggedObject>(
        *(FOLDS.getPluginA()));
  if (not _hasConstantB)
    _pluginb = std::make_shared<siconos::plugins::PluggedObject>(
        *(FOLDS.getPluginB()));
}

void siconos::modeling::FirstOrderLinearDS::initRhs(double time) {
  DEBUG_PRINT("init Rhs in FirstOrderLinearDS");
  computeRhs(time);  // If necessary, this will also compute A and b.
  if (!_jacxRhs)     // if not allocated with a set or anything else
  {
    if (_A && !_M) {  // if M is not defined, then A = jacobianRhsx, no memory
                      // allocation for that one.
      std::cout << "NOT IMPLEMENTED\n";
      //////// _jacxRhs = _A; /// WARNING FP: constructor is missing
    } else if (_A && _M) {
      auto tmp = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
      tmp->setZero();
      _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(*tmp);
      // else no allocation, jacobian is equal to 0.
    }
    computeJacobianRhsx(time);
  }
}

void siconos::modeling::FirstOrderLinearDS::updatePlugins(double time) {
  if (_M) computeM(time);
  if (_A) computeA(time);
  if (_b) computeb(time);
}

void siconos::modeling::FirstOrderLinearDS::setComputeAFunction(
    const std::string &pluginPath, const std::string &functionName) {
  if (!_A) _A = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _pluginA->setComputeFunction(pluginPath, functionName);
  _hasConstantA = false;
}

void siconos::modeling::FirstOrderLinearDS::setComputeAFunction(
    LDSPtrFunction fct) {
  if (!_A) _A = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _pluginA->setComputeFunction((void *)fct);
  _hasConstantA = false;
}

void siconos::modeling::FirstOrderLinearDS::setComputebFunction(
    const std::string &pluginPath, const std::string &functionName) {
  if (!_b) _b = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _pluginb->setComputeFunction(pluginPath, functionName);
  _hasConstantB = false;
}

void siconos::modeling::FirstOrderLinearDS::setComputebFunction(
    LDSPtrFunction fct) {
  if (!_b) _b = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _pluginb->setComputeFunction((void *)fct);
  _hasConstantB = false;
}

void siconos::modeling::FirstOrderLinearDS::clearComputebFunction() {
  _pluginb = std::make_shared<siconos::plugins::PluggedObject>();
}

void siconos::modeling::FirstOrderLinearDS::computeA(double time) {
  if (_A && _pluginA->fPtr) {
    ((computeAfct)_pluginA->fPtr)(time, _n, _n, &(*_A)(0, 0), _z->size(),
                                  &(*_z)(0));
  }
}

void siconos::modeling::FirstOrderLinearDS::computeb(double time) {
  if (_b && _pluginb->fPtr)
    ((LDSPtrFunction)_pluginb->fPtr)(time, _n, &(*_b)(0), _z->size(),
                                     &(*_z)(0));
}
/*This function is called only by LsodarOSI and eventDriven*/
void siconos::modeling::FirstOrderLinearDS::computeRhs(double time) {
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianfx

  *_x[1] = *_r;

  if (_A) {
    computeA(time);
    siconos::algebra::prod(*_A, *_x[0], *_x[1], false);
  }

  // compute and add b if required
  if (_b) {
    computeb(time);
    *_x[1] += *_b;
  }

  if (_M) {
    computeM(time);
    // allocate invM at the first call of the present function
    if (!_invM) _invM = std::make_shared<siconos::algebra::SiconosMatrix>(*_M);

    // _invM->Solve(*_x[1]);
    algebra::solveInPlace(*_invM, *(_x[1]));
  }
}

void siconos::modeling::FirstOrderLinearDS::computeJacobianRhsx(double time) {
  if (_A) {
    computeA(time);
    if (_M) {
      computeM(time);
      _jacxRhs->copyBlock(0, 0, _A);
      if (!_invM)
        _invM = std::make_shared<siconos::algebra::SiconosMatrix>(*_M);
      else if (_pluginM->fPtr)  // if M is plugged, invM must be updated
        *_invM = *_M;
      // solve MjacobianRhsx = A
      algebra::solveInPlace(*_invM, *(_jacxRhs->block(0, 0)));
    }
  }
  // else 0
}

void siconos::modeling::FirstOrderLinearDS::display(bool brief) const {
  std::cout << "=== Linear system display, " << _number << std::endl;
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
  std::cout << "M :" << std::endl;
  if (_M) {
    _M->display();
  } else
    std::cout << "M is identity" << std::endl;
  std::cout << "A :" << std::endl;
  if (_A)
    _A->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "b :" << std::endl;
  if (_b)
    _b->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "r :" << std::endl;
  if (_r)
    _r->display();
  else
    std::cout << "-> nullptr" << std::endl;
  if (_hasConstantA) {
    std::cout << "A is a time invariant matrix" << std::endl;
  }
  if (_hasConstantB) {
    std::cout << "b is a time invariant vector" << std::endl;
  }
  if (_pluginA->fPtr) {
    std::cout << "Has a plugin for A" << std::endl;
  }
  if (_pluginb->fPtr) {
    std::cout << "Has a plugin for b" << std::endl;
  }
  std::cout << "=============================" << std::endl;
}

void siconos::modeling::FirstOrderLinearDS::_zeroPlugin() {
  _pluginA = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginb = std::make_shared<siconos::plugins::PluggedObject>();
}

void siconos::modeling::FirstOrderLinearDS::setA(
    const siconos::algebra::SiconosMatrix &newA) {
  if (_A)
    *_A = newA;
  else
    _A = std::make_shared<siconos::algebra::SiconosMatrix>(newA);
  _hasConstantA = true;
}

void siconos::modeling::FirstOrderLinearDS::setb(
    const siconos::algebra::SiconosVector &b) {
  if (_b)
    *_b = b;
  else
    _b = std::make_shared<siconos::algebra::SiconosVector>(b);
  _hasConstantB = true;
}
