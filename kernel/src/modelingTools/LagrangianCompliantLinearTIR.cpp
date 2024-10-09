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
#include "LagrangianCompliantLinearTIR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"

// Minimum data (C as pointer) constructor
siconos::modeling::LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(
    std::shared_ptr<siconos::algebra::SiconosMatrix> C,
    std::shared_ptr<siconos::algebra::SiconosMatrix> D)
    : LagrangianR(RelationSubType::CompliantLinearTIR) {
  _jachq = C;
  _jachlambda = D;
}

// Constructor from a complete set of data
siconos::modeling::LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(
    std::shared_ptr<siconos::algebra::SiconosMatrix> C,
    std::shared_ptr<siconos::algebra::SiconosMatrix> D,
    std::shared_ptr<siconos::algebra::SiconosMatrix> F,
    std::shared_ptr<siconos::algebra::SiconosVector> e)
    : LagrangianR(RelationSubType::CompliantLinearTIR) {
  _jachq = C;
  _jachlambda = D;
  _F = F;
  _e = e;
}

// Minimum data (C, e as pointers) constructor
siconos::modeling::LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(
    std::shared_ptr<siconos::algebra::SiconosMatrix> C,
    std::shared_ptr<siconos::algebra::SiconosMatrix> D,
    std::shared_ptr<siconos::algebra::SiconosVector> e)
    : LagrangianR(RelationSubType::CompliantLinearTIR) {
  _jachq = C;
  _jachlambda = D;
  _e = e;
}

void siconos::modeling::LagrangianCompliantLinearTIR::initialize(Interaction& inter) {
  checkSize(inter);
}
void siconos::modeling::LagrangianCompliantLinearTIR::checkSize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto& DSlink = inter.linkToDSVariables();

  if (!(_jachq) || _jachq->size(1) != inter.getSizeOfDS() || _jachq->size(0) != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between H matrix and the interaction.");

  if ((_jachlambda) && (_jachlambda->size(0) != sizeY || _jachlambda->size(1) != sizeY))
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between D matrix and the interaction.");

  if ((_e) && _e->size() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between e vector and the dimension of the interaction.");

  auto sizeZ = DSlink[LagrangianR::z]->size();
  if ((_F) && (_F->size(0) != sizeZ || _F->size(1) != sizeZ))
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between F matrix and the interaction.");
}

void siconos::modeling::LagrangianCompliantLinearTIR::computeInput(double time,
                                                                   Interaction& inter,
                                                                   unsigned int level) {
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  auto& DSlink = inter.linkToDSVariables();
  // computation of p = Ht lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
}
void siconos::modeling::LagrangianCompliantLinearTIR::computeOutput(
    double time, Interaction& inter, unsigned int derivativeNumber) {
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  auto& lambda = *inter.lambda(derivativeNumber);
  auto& DSlink = inter.linkToDSVariables();

  siconos::algebra::matrixBlockVector_prod(*_jachq,
                                           *DSlink[LagrangianR::q0 + derivativeNumber], y);
  y += *_jachlambda * lambda;
  siconos::algebra::prod(*_jachlambda, lambda, y, false);

  if (derivativeNumber == 0) {
    if (_e) y += *_e;
    if (_F) siconos::algebra::matrixBlockVector_prod(*_F, *DSlink[LagrangianR::z], y, false);
  }
}

void siconos::modeling::LagrangianCompliantLinearTIR::display() const {
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== " << std::endl;
  std::cout << " C: " << std::endl;
  if (_jachq)
    _jachq->display();
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << " e: " << std::endl;
  if (_e)
    _e->display();
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << " D: " << std::endl;
  if (_jachlambda)
    _jachlambda->display();
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << " F: " << std::endl;
  if (_F)
    _F->display();
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << "===================================== " << std::endl;
}
