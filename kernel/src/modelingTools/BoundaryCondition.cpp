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

siconos::modeling::BoundaryCondition::BoundaryCondition(const Indices& newVelocityIndices)
    : velocityIndices_{newVelocityIndices},
      prescribedVelocity_{
          std::make_shared<siconos::algebra::SiconosVector>(velocityIndices_.size())},
      prescribedVelocityOld_{
          std::make_shared<siconos::algebra::SiconosVector>(velocityIndices_.size())} {
  prescribedVelocity_->setZero();
  prescribedVelocityOld_->setZero();
}

siconos::modeling::BoundaryCondition::BoundaryCondition(Indices&& newVelocityIndices)
    : velocityIndices_{std::move(newVelocityIndices)},
      prescribedVelocity_{
          std::make_shared<siconos::algebra::SiconosVector>(velocityIndices_.size())},
      prescribedVelocityOld_{
          std::make_shared<siconos::algebra::SiconosVector>(velocityIndices_.size())} {
  prescribedVelocity_->setZero();
  prescribedVelocityOld_->setZero();
}

siconos::modeling::BoundaryCondition::BoundaryCondition(
    Indices&& newVelocityIndices,
    std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues)
    : velocityIndices_{std::move(newVelocityIndices)},
      prescribedVelocity_(newVelocityValues),
      prescribedVelocityOld_{
          std::make_shared<siconos::algebra::SiconosVector>(*newVelocityValues)} {
  if (newVelocityIndices.size() != (size_t)newVelocityValues->size())
    THROW_EXCEPTION(
        "siconos::modeling::BoundaryCondition::BoundaryCondition  "
        "constructor. "
        "velocityIndices and prescribedVelocity must have the same size");
}

siconos::modeling::BoundaryCondition::BoundaryCondition(
    const Indices& newVelocityIndices,
    std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues)
    : velocityIndices_{newVelocityIndices},
      prescribedVelocity_(newVelocityValues),
      prescribedVelocityOld_{
          std::make_shared<siconos::algebra::SiconosVector>(*newVelocityValues)} {
  if (newVelocityIndices.size() != (size_t)newVelocityValues->size())
    THROW_EXCEPTION(
        "siconos::modeling::BoundaryCondition::BoundaryCondition  "
        "constructor. "
        "velocityIndices and prescribedVelocity must have the same size");
}

void siconos::modeling::BoundaryCondition::setComputePrescribedVelocityFunction(
    const siconos::modeling::func_prototypes::FunctionS_V& func) {
  if (!prescribedVelocity_) {
    prescribedVelocity_ = std::make_shared<siconos::algebra::SiconosVector>(size());
  }
  computePrescribedVelocity_ = func;
}

void siconos::modeling::BoundaryCondition::computePrescribedVelocity(double time) {
  if (computePrescribedVelocity_) computePrescribedVelocity_(time, *prescribedVelocity_);
}

void siconos::modeling::BoundaryCondition::display() const {
  std::cout << "=====  BoundaryCondition display ===== " << std::endl;

  tools::print("- Indices on which boundary conditions are applied:\n ", velocityIndices_);

  std::cout << "- velocities : " << std::endl;
  if (prescribedVelocity_) siconos::algebra::print(*prescribedVelocity_);
  std::cout << "=========================================================== " << std::endl;
}

void siconos::modeling::BoundaryCondition::appendIndex(
    siconos::algebra::SiconosVector::Index ind) {
  if (find(velocityIndices_.begin(), velocityIndices_.end(), ind) == velocityIndices_.end()) {
    velocityIndices_.push_back(ind);
    prescribedVelocity_->resize(size());
    prescribedVelocityOld_->resize(size());
  }
}
