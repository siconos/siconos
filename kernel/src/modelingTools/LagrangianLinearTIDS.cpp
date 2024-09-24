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
#include "LagrangianLinearTIDS.hpp"

#include "BlockMatrix.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "siconos_debug.h"

siconos::modeling::LagrangianLinearTIDS::LagrangianLinearTIDS(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0,
    Eigen::Ref<siconos::algebra::SiconosMatrix> newmass)
    : LagrangianDS(q0, v0) {
  hasConstantMass_ = true;
  hasMass_ = true;
  computemass_ = nullptr;
  mass_view_ = std::make_shared<siconos::algebra::MapType>(newmass.data(), newmass.rows(),
                                                           newmass.cols());
};

void siconos::modeling::LagrangianLinearTIDS::initRhs(double time) {
  LagrangianDS::initRhs(time);

  // jacobianRhsx
  if (_K) {
    //  bloc10 of jacobianX is solution of Mass*Bloc10 = K
    if (!_rhsMatrices[jacobianXBloc10_])
      // We assume the stiffness matrix is constant
      _rhsMatrices[jacobianXBloc10_] =
          std::make_shared<siconos::algebra::SiconosMatrix>(-1 * *_K);
    algebra::solveInPlace(*_inverseMass, *_rhsMatrices[jacobianXBloc10_]);
  } else
    _rhsMatrices[jacobianXBloc10_] = _rhsMatrices[zeroMatrix_];

  if (_C) {
    //  bloc11 of jacobianX is solution of Mass*Bloc11 = C
    if (!_rhsMatrices[jacobianXBloc11_])
      // We assume the damping matrix is constant
      _rhsMatrices[jacobianXBloc11_] =
          std::make_shared<siconos::algebra::SiconosMatrix>(-1 * *_C);
    algebra::solveInPlace(*_inverseMass, *_rhsMatrices[jacobianXBloc11_]);
  } else
    _rhsMatrices[jacobianXBloc11_] = _rhsMatrices[zeroMatrix_];

  if (_C || _K)
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[zeroMatrix_], _rhsMatrices[idMatrix_], _rhsMatrices[jacobianXBloc10_],
        _rhsMatrices[jacobianXBloc11_]);
}

void siconos::modeling::LagrangianLinearTIDS::setK(const Matrix &newValue) {
  if (newValue.size(0) != ndof_ || newValue.size(1) != ndof_)
    THROW_EXCEPTION("LagrangianLinearTIDS - setK: inconsistent input matrix size ");

  if (!_K)
    _K = std::make_shared<Matrix>(newValue);
  else
    *_K = newValue;
}

void siconos::modeling::LagrangianLinearTIDS::setKPtr(std::shared_ptr<Matrix> newPtr) {
  if (newPtr->size(0) != ndof_ || newPtr->size(1) != ndof_)
    THROW_EXCEPTION("LagrangianLinearTIDS - setKPtr: inconsistent input matrix size ");
  _K = newPtr;
}

void siconos::modeling::LagrangianLinearTIDS::setC(const Matrix &newValue) {
  if (newValue.size(0) != ndof_ || newValue.size(1) != ndof_)
    THROW_EXCEPTION("LagrangianLinearTIDS - setC: inconsistent input matrix size ");

  if (!_C)
    _C = std::make_shared<Matrix>(newValue);
  else
    *_C = newValue;
}

void siconos::modeling::LagrangianLinearTIDS::setCPtr(std::shared_ptr<Matrix> newPtr) {
  if (newPtr->size(0) != ndof_ || newPtr->size(1) != ndof_)
    THROW_EXCEPTION("LagrangianLinearTIDS - setCPtr: inconsistent input matrix size ");

  _C = newPtr;
}

void siconos::modeling::LagrangianLinearTIDS::display(bool brief) const {
  LagrangianDS::display();
  std::cout << "===== Lagrangian Linear Time Invariant System display ===== " << std::endl;
  std::cout << "- Mass Matrix M : " << std::endl;
  if (mass_view_)
    std::cout << *mass_view_ << "\n";
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
    _forces = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  } else
    _forces->zero();

  if (fext_view_) {
    computeFext(time);
    *_forces += *fext_view_;
  }
  if (_K) *_forces -= siconos::algebra::prod(*_K, *q2);
  if (_C) *_forces -= siconos::algebra::prod(*_C, *v2);
}
