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
#include "LagrangianLinearTIDS.hpp"

#include "BlockMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "siconos_debug.h"

// --- Constructor from a initial conditions and matrix-operators
siconos::modeling::LagrangianLinearTIDS::LagrangianLinearTIDS(
    std::shared_ptr<siconos::algebra::SiconosVector> newQ0,
    std::shared_ptr<siconos::algebra::SiconosVector> newVelocity0,
    std::shared_ptr<siconos::algebra::SiconosMatrix> newMass,
    std::shared_ptr<siconos::algebra::SiconosMatrix> newK,
    std::shared_ptr<siconos::algebra::SiconosMatrix> newC)
    : LagrangianDS(newQ0, newVelocity0, newMass) {
  assert(
      (newK->size(0) == _ndof && newK->size(1) == _ndof) &&
      "LagrangianLinearTIDS - constructor from data, inconsistent size between K and _ndof");

  assert((newC->size(0) == _ndof && newC->size(1) == _ndof) &&
         "LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  _K = newK;
  _C = newC;
}

void siconos::modeling::LagrangianLinearTIDS::initRhs(double time) {

  LagrangianDS::initRhs(time);

  // jacobianRhsx
  if (_K) {
    //  bloc10 of jacobianX is solution of Mass*Bloc10 = K
    if (!_rhsMatrices[jacobianXBloc10_])
      _rhsMatrices[jacobianXBloc10_] =
          std::make_shared<siconos::algebra::SimpleMatrix>(-1 * *_K);
    algebra::solveInPlace(*_inverseMass, *_rhsMatrices[jacobianXBloc10_]);
  } else
    _rhsMatrices[jacobianXBloc10_] = _rhsMatrices[zeroMatrix_];

  if (_C) {
    //  bloc11 of jacobianX is solution of Mass*Bloc11 = C
    if (!_rhsMatrices[jacobianXBloc11_])
      _rhsMatrices[jacobianXBloc11_] =
          std::make_shared<siconos::algebra::SimpleMatrix>(-1 * *_C);
    algebra::solveInPlace(*_inverseMass, *_rhsMatrices[jacobianXBloc11_]);
  } else
    _rhsMatrices[jacobianXBloc11_] = _rhsMatrices[zeroMatrix_];

  if (_C || _K)
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[zeroMatrix_], _rhsMatrices[idMatrix_], _rhsMatrices[jacobianXBloc10_],
        _rhsMatrices[jacobianXBloc11_]);
}

void siconos::modeling::LagrangianLinearTIDS::setK(
    const siconos::algebra::SiconosMatrix& newValue) {
  if (newValue.size(0) != _ndof || newValue.size(1) != _ndof)
    THROW_EXCEPTION("LagrangianLinearTIDS - setK: inconsistent input matrix size ");

  if (!_K)
    _K = std::make_shared<siconos::algebra::SimpleMatrix>(newValue);
  else
    *_K = newValue;
}

void siconos::modeling::LagrangianLinearTIDS::setKPtr(
    std::shared_ptr<siconos::algebra::SiconosMatrix> newPtr) {
  if (newPtr->size(0) != _ndof || newPtr->size(1) != _ndof)
    THROW_EXCEPTION("LagrangianLinearTIDS - setKPtr: inconsistent input matrix size ");
  _K = newPtr;
}

void siconos::modeling::LagrangianLinearTIDS::setC(
    const siconos::algebra::SiconosMatrix& newValue) {
  if (newValue.size(0) != _ndof || newValue.size(1) != _ndof)
    THROW_EXCEPTION("LagrangianLinearTIDS - setC: inconsistent input matrix size ");

  if (!_C)
    _C = std::make_shared<siconos::algebra::SimpleMatrix>(newValue);
  else
    *_C = newValue;
}

void siconos::modeling::LagrangianLinearTIDS::setCPtr(
    std::shared_ptr<siconos::algebra::SiconosMatrix> newPtr) {
  if (newPtr->size(0) != _ndof || newPtr->size(1) != _ndof)
    THROW_EXCEPTION("LagrangianLinearTIDS - setCPtr: inconsistent input matrix size ");

  _C = newPtr;
}

void siconos::modeling::LagrangianLinearTIDS::display(bool brief) const {
  LagrangianDS::display();
  std::cout << "===== Lagrangian Linear Time Invariant System display ===== " << std::endl;
  std::cout << "- Mass Matrix M : " << std::endl;
  if (_mass)
    _mass->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- Stiffness Matrix K : " << std::endl;
  if (_K)
    _K->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- Viscosity Matrix C : " << std::endl;
  if (_C)
    _C->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "=========================================================== " << std::endl;
}

void siconos::modeling::LagrangianLinearTIDS::computeForces(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q2,
    std::shared_ptr<siconos::algebra::SiconosVector> v2) {
  DEBUG_PRINT(
      "siconos::modeling::LagrangianLinearTIDS::computeForces(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> q2, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v2) \n");
  if (!_forces) {
    _forces = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  } else
    _forces->zero();

  if (_fExt) {
    computeFExt(time);
    *_forces += *_fExt;
  }
  if (_K) *_forces -= siconos::algebra::prod(*_K, *q2);
  if (_C) *_forces -= siconos::algebra::prod(*_C, *v2);
}
