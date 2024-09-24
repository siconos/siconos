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
#include "SecondOrderDS.hpp"

#include "BoundaryCondition.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "siconos_debug.h"

void siconos::modeling::SecondOrderDS::setBoundaryConditions(
    std::shared_ptr<siconos::modeling::BoundaryCondition> newbd) {
  if (_boundaryConditions) {
    std::cout << "Warning : SecondOrderDS::setBoundaryConditions. old boundary conditions "
                 "were pre-existing. \n";
  }
  _boundaryConditions = newbd;
  _reactionToBoundaryConditions =
      std::make_shared<siconos::algebra::SiconosVector>(_boundaryConditions->size());
};

void siconos::modeling::SecondOrderDS::setConstantMass(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for mass
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  mass_internal_storage_ = nullptr;

  mass_view_ = std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasMass_ = true;
  hasConstantMass_ = true;
  computemass_ = nullptr;
}

void siconos::modeling::SecondOrderDS::setComputeMassFunction(MassFunction new_func) {
  // Ensure that memory is properly allocated for mass_
  if (!mass_internal_storage_) {
    mass_internal_storage_ = std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  mass_view_ = std::make_shared<siconos::algebra::MapType>(mass_internal_storage_->data(),
                                                           ndof_, ndof_);
  mass_view_->setZero();
  hasMass_ = true;
  hasConstantMass_ = false;
  computemass_ = new_func;
}

void siconos::modeling::SecondOrderDS::computeMass(
    Eigen::Ref<siconos::algebra::SiconosVector> position, double time) {
  if (computemass_) {
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    // std::span<double> sppos(position.data(), ndof_);
    // std::span<double> spmat(mass_internal_storage_->data(),
    //                         mass_internal_storage_->size());
    computemass_(position, time, *mass_view_);
  }
}
