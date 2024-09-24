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
#include "LagrangianLinearDiagonalDS.hpp"

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "siconos_debug.h"

// --- Constructor for the complete system
siconos::modeling::LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity0,
    std::shared_ptr<siconos::algebra::SiconosVector> stiffness,
    std::shared_ptr<siconos::algebra::SiconosVector> damping,
    std::shared_ptr<siconos::algebra::SiconosVector> mass)
    : LagrangianDS(q0, velocity0) {
  _mass_diag = std::make_shared<siconos::algebra::MapVectorType>(mass->data(), dimension());
  _stiffness =
      std::make_shared<siconos::algebra::MapVectorType>(stiffness->data(), dimension());
  _damping = std::make_shared<siconos::algebra::MapVectorType>(damping->data(), dimension());
}

// --- Constructor for the complete system with identity mass matrix
siconos::modeling::LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity0,
    std::shared_ptr<siconos::algebra::SiconosVector> stiffness,
    std::shared_ptr<siconos::algebra::SiconosVector> damping)
    : LagrangianDS(q0, velocity0) {
  _stiffness =
      std::make_shared<siconos::algebra::MapVectorType>(stiffness->data(), dimension());
  _damping = std::make_shared<siconos::algebra::MapVectorType>(damping->data(), dimension());
}

// --- Constructor for the undamped system with identity mass matrix
siconos::modeling::LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity0,
    std::shared_ptr<siconos::algebra::SiconosVector> stiffness)
    : LagrangianDS(q0, velocity0) {
  _stiffness =
      std::make_shared<siconos::algebra::MapVectorType>(stiffness->data(), dimension());
}

void siconos::modeling::LagrangianLinearDiagonalDS::initRhs(double time) {
  THROW_EXCEPTION(
      "siconos::modeling::LagrangianLinearDiagonalDS::initRhs - not yet implemented for "
      "LagrangianLinearDiagonalDS.");
}

void siconos::modeling::LagrangianLinearDiagonalDS::computeForces(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q2,
    std::shared_ptr<siconos::algebra::SiconosVector> v2) {
  DEBUG_PRINT(
      "LagrangianLinearTIDS::computeForces(double time, "
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

  if (_stiffness)
    for (decltype(ndof_) i = 0; i < ndof_; ++i) (*_forces)(i) -= (*_stiffness)(i) * (*q2)(i);
  if (_damping)
    for (decltype(ndof_) i = 0; i < ndof_; ++i) (*_forces)(i) -= (*_damping)(i) * (*v2)(i);

  // if (_stiffness)
  // {
  //   DenseVect & stiff = *_stiffness->dense();
  //   noalias(*force->dense()) -= element_prod(stiff, *q2->dense());
  // }
  // if (_damping)
  // {
  //   DenseVect & damp = *_damping->dense();
  //   noalias(*force->dense()) -= element_prod(damp, *v2->dense());
  // }

  // Test this:
  // DenseVect& omega = *_stiffness.dense();
  // DenseVect& sigma = *_damping.dense();
  // DenseVect& res = *_forces.dense();
  // DenseVect& qvect = *q2.dense();
  // DenseVect& vvect = *v2.dense();
  // if (_stiffness)
  //   std::transform(omega.begin(), omega.end(),
  //                  qvect.begin(), res.begin(), std::multiplies<double>() );
  // if (_damping)
  //   std::transform(sigma.begin(), sigma.end(),
  //                  vvect.begin(), res.begin(), std::multiplies<double>() );
  // *forces *= -1.;
}

void siconos::modeling::LagrangianLinearDiagonalDS::display(bool brief) const {
  LagrangianDS::display();
  std::cout << "===== Lagrangian Linear Diagonal System display ===== \n ";
  std::cout << "- Mass Matrix M : \n";
  if (mass_view_)
    std::cout << *mass_view_ << "\n";
  else
    std::cout << "-> nullptr \n";
  std::cout << "- Stiffness Matrix K : \n";
  if (_stiffness)
    _stiffness->display();
  else
    std::cout << "-> nullptr\n";
  std::cout << "- Viscosity Matrix C : \n";
  if (_damping)
    _damping->display();
  else
    std::cout << "-> nullptr\n";
  std::cout << "=========================================================== \n";
}
