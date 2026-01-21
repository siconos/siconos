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

#include "BoundaryCondition.hpp"

#include <iostream>

#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"

siconos::modeling::BoundaryCondition::BoundaryCondition(
    Indices newVelocityIndices,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& newVelocityValues)
    : velocityIndices_{std::move(newVelocityIndices)},
      prescribedVelocity_(newVelocityValues),
      prescribedVelocityOld_{newVelocityValues} {
  if (velocityIndices_.size() != prescribedVelocity_.size())
    THROW_EXCEPTION(
        "siconos::modeling::BoundaryCondition::BoundaryCondition  "
        "constructor. "
        "velocityIndices and prescribedVelocity must have the same size");
}

void siconos::modeling::BoundaryCondition::setComputePrescribedVelocityFunction(
    const siconos::modeling::func_prototypes::FunctionS_V& func) {
  assert(prescribedVelocity_.size() == velocityIndices_.size());
  computePrescribedVelocity_ = func;
}

void siconos::modeling::BoundaryCondition::computePrescribedVelocity(double time) {
  if (computePrescribedVelocity_) computePrescribedVelocity_(time, prescribedVelocity_);
}

void siconos::modeling::BoundaryCondition::display() const {
  std::cout << "=====  BoundaryCondition display ===== \n";
  tools::print("- Indices on which boundary conditions are applied:\n ", velocityIndices_);

  std::cout << "- velocities : \n";
  siconos::algebra::print(prescribedVelocity_);
  std::cout << "=========================================================== \n";
}

void siconos::modeling::BoundaryCondition::appendIndex(
    siconos::algebra::SiconosVector::Index ind, double imposedVelocity) {
  if (find(velocityIndices_.begin(), velocityIndices_.end(), ind) == velocityIndices_.end()) {
    velocityIndices_.push_back(ind);
    prescribedVelocity_.conservativeResize(size());
    prescribedVelocity_(size() - 1) = imposedVelocity;
    prescribedVelocityOld_.conservativeResize(size());
    prescribedVelocity_(size() - 1) = imposedVelocity;
  }
}
