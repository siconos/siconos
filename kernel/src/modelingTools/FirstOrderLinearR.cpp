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
#include "FirstOrderLinearR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <cassert>

#include "siconos_debug.h"

// Constructor with C and B plug-in names
siconos::modeling::FirstOrderLinearR::FirstOrderLinearR(const std::string &Cname,
                                                        const std::string &Bname)
    : FirstOrderR(RelationSubType::LinearR) {
  // Warning: we cannot allocate memory for C/D matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the
  // user-defined function.
  _pluginJachx->setComputeFunction(Cname);
  _pluginJacglambda->setComputeFunction(Bname);
}

// Constructor from a complete set of data (plugin)
siconos::modeling::FirstOrderLinearR::FirstOrderLinearR(const std::string &Cname,
                                                        const std::string &Dname,
                                                        const std::string &Fname,
                                                        const std::string &Ename,
                                                        const std::string &Bname)
    : FirstOrderR(RelationSubType::LinearR) {
  _pluginJachx->setComputeFunction(Cname);
  _pluginJachlambda->setComputeFunction(Dname);
  _pluginJacglambda->setComputeFunction(Bname);
  _pluginf->setComputeFunction(Fname);
  _plugine->setComputeFunction(Ename);
}

// Minimum data (C, B as pointers) constructor
siconos::modeling::FirstOrderLinearR::FirstOrderLinearR(
    std::shared_ptr<siconos::algebra::SimpleMatrix> C,
    std::shared_ptr<siconos::algebra::SimpleMatrix> B)
    : FirstOrderR(RelationSubType::LinearR) {
  _C = C;
  _B = B;
}

// // Constructor from a complete set of data
siconos::modeling::FirstOrderLinearR::FirstOrderLinearR(
    std::shared_ptr<siconos::algebra::SimpleMatrix> C,
    std::shared_ptr<siconos::algebra::SimpleMatrix> D,
    std::shared_ptr<siconos::algebra::SimpleMatrix> F,
    std::shared_ptr<siconos::algebra::SiconosVector> E,
    std::shared_ptr<siconos::algebra::SimpleMatrix> B)
    : FirstOrderR(RelationSubType::LinearR) {
  _C = C;
  _B = B;
  _D = D;
  _F = F;
  _e = E;
}

void siconos::modeling::FirstOrderLinearR::initialize(Interaction &inter) {
  FirstOrderR::initialize(inter);

  // get interesting size
  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();

  auto &DSlink = inter.linkToDSVariables();
  auto sizeZ = DSlink[FirstOrderR::z]->size();
  auto &relationMat = inter.relationMatrices();
  auto &relationVec = inter.relationVectors();

  if (!_C && _pluginJachx->fPtr)
    relationMat[FirstOrderR::mat_C] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeX);
  if (!_D && _pluginJachlambda->fPtr)
    relationMat[FirstOrderR::mat_D] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeY);
  if (!_B && _pluginJacglambda->fPtr)
    relationMat[FirstOrderR::mat_B] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeX, sizeY);
  if (!_F && _pluginf->fPtr)
    relationMat[FirstOrderR::mat_F] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeZ);
  if (!_e && _plugine->fPtr)
    relationVec[FirstOrderR::e] = std::make_shared<siconos::algebra::SiconosVector>(sizeY);

  checkSize(inter);
}

void siconos::modeling::FirstOrderLinearR::checkSize(Interaction &inter) {
  auto &DSlink = inter.linkToDSVariables();

  // get inter and ds sizes
  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();
  auto sizeZ = DSlink[FirstOrderR::z]->size();

  // Check if various operators sizes are consistent.
  // Reference: interaction.

  if (_C) {
    if (_C->size(0) == 0)
      _C->resize(sizeX, sizeY);
    else
      assert(
          (_C->size(0) == sizeY && _C->size(1) == sizeX) &&
          "siconos::modeling::FirstOrderLinearR::initialize , inconsistent size between C and "
          "Interaction.");
  }
  if (_B) {
    if (_B->size(0) == 0)
      _B->resize(sizeY, sizeX);
    else
      assert(
          (_B->size(1) == sizeY && _B->size(0) == sizeX) &&
          "siconos::modeling::FirstOrderLinearR::initialize , inconsistent size between B and "
          "interaction.");
  }

  if (_D) {
    if (_D->size(0) == 0)
      _D->resize(sizeY, sizeY);
    else
      assert((_D->size(0) == sizeY || _D->size(1) == sizeY) &&
             "siconos::modeling::FirstOrderLinearR::initialize , inconsistent size between C "
             "and D.");
  }

  if (_F) {
    if (_F->size(0) == 0)
      _F->resize(sizeY, sizeZ);
    else
      assert(((_F->size(0) == sizeY) && (_F->size(1) == sizeZ)) &&
             "siconos::modeling::FirstOrderLinearR::initialize , inconsistent size between C "
             "and F.");
  }

  if (_e) {
    if (_e->size() == 0)
      _e->resize(sizeY);
    else
      assert(_e->size() == sizeY &&
             "siconos::modeling::FirstOrderLinearR::initialize , inconsistent size between C "
             "and e.");
  }
}
void siconos::modeling::FirstOrderLinearR::computeC(double time,
                                                    siconos::algebra::BlockVector &z,
                                                    siconos::algebra::SimpleMatrix &C) {
  if (_pluginJachx->fPtr) {
    auto zp = z.toSiconosVector();
    ((FOMatPtr1)(_pluginJachx->fPtr))(time, C.size(0), C.size(1), &(C)(0, 0), zp->size(),
                                      &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderLinearR::computeD(double time,
                                                    siconos::algebra::BlockVector &z,
                                                    siconos::algebra::SimpleMatrix &D) {
  if (_pluginJachlambda->fPtr) {
    auto zp = z.toSiconosVector();
    ((FOMatPtr1)(_pluginJachlambda->fPtr))(time, D.size(0), D.size(1), &(D)(0, 0), zp->size(),
                                           &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderLinearR::computeF(double time,
                                                    siconos::algebra::BlockVector &z,
                                                    siconos::algebra::SimpleMatrix &F) {
  if (_pluginf->fPtr) {
    auto zp = z.toSiconosVector();
    ((FOMatPtr1)(_pluginf->fPtr))(time, F.size(0), F.size(1), &(F)(0, 0), zp->size(),
                                  &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderLinearR::computee(double time,
                                                    siconos::algebra::BlockVector &z,
                                                    siconos::algebra::SiconosVector &e) {
  if (_plugine->fPtr) {
    auto zp = z.toSiconosVector();
    ((FOVecPtr)_plugine->fPtr)(time, e.size(), &(e)(0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderLinearR::computeB(double time,
                                                    siconos::algebra::BlockVector &z,
                                                    siconos::algebra::SimpleMatrix &B) {
  if (_pluginJacglambda->fPtr) {
    auto zp = z.toSiconosVector();
    ((FOMatPtr1)_pluginJacglambda->fPtr)(time, B.size(0), B.size(1), &(B)(0, 0), zp->size(),
                                         &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderLinearR::computeh(
    double time, const siconos::algebra::BlockVector &x,
    const siconos::algebra::SiconosVector &lambda, siconos::algebra::BlockVector &z,
    siconos::algebra::SiconosVector &y) {
  if (_C) {
    computeC(time, z, *_C);
    siconos::algebra::prod(*_C, x, y, true);
  } else
    y.setZero();

  if (_D) {
    computeD(time, z, *_D);
    siconos::algebra::prod(*_D, lambda, y, false);
  }
  if (_F) {
    computeF(time, z, *_F);
    siconos::algebra::prod(*_F, z, y, false);
  }
  if (_e) {
    computee(time, z, *_e);
    y += *_e;
  }
}

void siconos::modeling::FirstOrderLinearR::computeOutput(double time, Interaction &inter,
                                                         unsigned int level) {
  DEBUG_BEGIN("siconos::modeling::FirstOrderLinearR::computeOutput \n");
  siconos::algebra::SiconosVector &y = *inter.y(level);
  siconos::algebra::SiconosVector &lambda = *inter.lambda(level);
  auto &DSlink = inter.linkToDSVariables();
  computeh(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], y);
  DEBUG_END("siconos::modeling::FirstOrderLinearR::computeOutput \n");
}

void siconos::modeling::FirstOrderLinearR::computeg(
    double time, const siconos::algebra::SiconosVector &lambda,
    siconos::algebra::BlockVector &z, siconos::algebra::BlockVector &r) {
  computeB(time, z, *_B);
  siconos::algebra::prod(*_B, lambda, r, false);
}

void siconos::modeling::FirstOrderLinearR::computeInput(double time, Interaction &inter,
                                                        unsigned int level) {
  auto &DSlink = inter.linkToDSVariables();
  siconos::algebra::BlockVector &z = *DSlink[FirstOrderR::z];
  computeg(time, *inter.lambda(level), z, *DSlink[FirstOrderR::r]);
}

void siconos::modeling::FirstOrderLinearR::display() const {
  std::cout << " ===== Linear relation display ===== "
            << "\n";
  std::cout << "| C "
            << "\n";
  if (_C)
    _C->display();
  else
    std::cout << "->nullptr\n";
  std::cout << "| D "
            << "\n";
  if (_D)
    _D->display();
  else
    std::cout << "->nullptr\n";
  std::cout << "| F "
            << "\n";
  if (_F)
    _F->display();
  else
    std::cout << "->nullptr\n";
  std::cout << "| e "
            << "\n";
  if (_e)
    _e->display();
  else
    std::cout << "->nullptr\n";
  std::cout << "| B "
            << "\n";
  if (_B)
    _B->display();
  else
    std::cout << "->nullptr\n";
  std::cout << " ================================================== \n";
}
