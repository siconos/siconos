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
#include "FirstOrderLinearTIR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// Minimum data (C, B as pointers) constructor
siconos::modeling::FirstOrderLinearTIR::FirstOrderLinearTIR(
    std::shared_ptr<siconos::algebra::SiconosMatrix> C,
    std::shared_ptr<siconos::algebra::SiconosMatrix> B)
    : FirstOrderR(RelationSubType::LinearTIR) {
  _C = C;
  _B = B;
}

// Constructor from a complete set of data
siconos::modeling::FirstOrderLinearTIR::FirstOrderLinearTIR(
    std::shared_ptr<siconos::algebra::SiconosMatrix> C,
    std::shared_ptr<siconos::algebra::SiconosMatrix> D,
    std::shared_ptr<siconos::algebra::SiconosMatrix> F,
    std::shared_ptr<siconos::algebra::SiconosVector> e,
    std::shared_ptr<siconos::algebra::SiconosMatrix> B)
    : FirstOrderR(RelationSubType::LinearTIR) {
  _C = C;
  _B = B;
  _D = D;
  _F = F;
  _e = e;
}

void siconos::modeling::FirstOrderLinearTIR::initialize(Interaction &inter) {
  DEBUG_PRINT("siconos::modeling::FirstOrderLinearTIR::initialize(Interaction & inter)\n");

  FirstOrderR::initialize(inter);  // ?

  if (!_C)
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderLinearTIR::initialize() C is null and is a required "
        "input.");
  if (!_B)
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderLinearTIR::initialize() B is null and is a required "
        "input.");

  checkSize(inter);
}

void siconos::modeling::FirstOrderLinearTIR::checkSize(Interaction &inter) {
  DEBUG_PRINT("siconos::modeling::FirstOrderLinearTIR::checkSize(Interaction & inter)\n");
  DEBUG_PRINTF("_C->size(0) = %i,\t inter.dimension() = %i\n ", _C->size(0),
               inter.dimension());
  DEBUG_PRINTF("_C->size(1) = %i,\t inter.getSizeOfDS() = %i\n ", _C->size(1),
               inter.getSizeOfDS());

  assert(
      (_C->size(0) == inter.dimension() && _C->size(1) == inter.getSizeOfDS()) &&
      "siconos::modeling::FirstOrderLinearTIR::initialize , inconsistent size between C and "
      "Interaction sizes.");

  assert(
      (_B->size(1) == inter.dimension() && _B->size(0) == inter.getSizeOfDS()) &&
      "siconos::modeling::FirstOrderLinearTIR::initialize , inconsistent size between B and "
      "interaction sizes.");

  // C and B are the minimum inputs. The others may remain null.

  if (_D)
    assert(
        (_D->size(0) == inter.dimension() && _D->size(1) == inter.dimension()) &&
        "siconos::modeling::FirstOrderLinearTIR::initialize , inconsistent size between D and "
        "interaction sizes");

  DEBUG_EXPR(if (_F) _F->display(); (inter.linkToDSVariables())[FirstOrderR::z]->display(););

  if (_F)
    assert(((_F->size(0) == inter.dimension()) &&
            (_F->size(1) == (inter.linkToDSVariables())[FirstOrderR::z]->size())) &&
           "siconos::modeling::FirstOrderLinearTIR::initialize , inconsistent size between F "
           "and z.");
  if (_e)
    assert(_e->size() == inter.dimension() &&
           "siconos::modeling::FirstOrderLinearTIR::initialize , inconsistent size between C "
           "and e.");
}

void siconos::modeling::FirstOrderLinearTIR::computeh(
    const siconos::algebra::BlockVector &x, const siconos::algebra::SiconosVector &lambda,
    siconos::algebra::BlockVector &z, siconos::algebra::SiconosVector &y) {
  // if (_C) C must be allocated. Checksize is there to ensure it.
  siconos::algebra::matrixBlockVector_prod(*_C, x, y, true);
  // else
  //   y.setZero();

  if (_D) siconos::algebra::prod(*_D, lambda, y, false);

  if (_e) y += *_e;

  if (_F) siconos::algebra::matrixBlockVector_prod(*_F, z, y, false);
}

void siconos::modeling::FirstOrderLinearTIR::computeOutput(double time, Interaction &inter,
                                                           unsigned int level) {
  // We get y and lambda of the interaction (pointers)
  siconos::algebra::SiconosVector &y = *inter.y(level);
  siconos::algebra::SiconosVector &lambda = *inter.lambda(level);
  auto &DSlink = inter.linkToDSVariables();
  computeh(*DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], y);
}

void siconos::modeling::FirstOrderLinearTIR::computeg(
    const siconos::algebra::SiconosVector &lambda, siconos::algebra::BlockVector &r) {
  siconos::algebra::matrixVector_prod_toBlock(*_B, lambda, r, false);
}

void siconos::modeling::FirstOrderLinearTIR::computeInput(double time, Interaction &inter,
                                                          unsigned int level) {
  DEBUG_BEGIN(
      "siconos::modeling::FirstOrderLinearTIR::computeInput(double time, Interaction& "
      "inter, unsigned int level)\n")
  auto &DSlink = inter.linkToDSVariables();
  DEBUG_EXPR(inter.lambda(level)->display(););
  DEBUG_EXPR(DSlink[FirstOrderR::r]->display(););
  computeg(*inter.lambda(level), *DSlink[FirstOrderR::r]);
  DEBUG_END(
      "siconos::modeling::FirstOrderLinearTIR::computeInput(double time, Interaction& "
      "inter, unsigned int level)\n")
}

void siconos::modeling::FirstOrderLinearTIR::display() const {
  std::cout << " ===== Linear Time Invariant relation display =====\n";
  std::cout << "| C\n";
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
  std::cout << "| F \n";
  if (_F)
    _F->display();
  else
    std::cout << "->nullptr\n";
  std::cout << "| e\n";
  if (_e)
    _e->display();
  else
    std::cout << "->nullptr\n";
  std::cout << "| B\n";
  if (_B)
    _B->display();
  else
    std::cout << "->nullptr\n";
  std::cout << " ==================================================\n";
}
