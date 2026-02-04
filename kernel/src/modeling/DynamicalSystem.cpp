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
#include "DynamicalSystem.hpp"

#include "SiconosException.hpp"
#include "SiconosMemory.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include <iostream>

#include "siconos_debug.h"

size_t siconos::modeling::DynamicalSystem::__count = 0;

siconos::modeling::DynamicalSystem::DynamicalSystem(const DynamicalSystem& ds)
    : x_size_(ds.x_size_), stepsInMemory_(ds.stepsInMemory_) {
  // The following data should always be initialize
  if (ds.hasX0()) {
    x0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ds.x0());
  }
  if (ds.r()) rVector_ = std::make_shared<siconos::algebra::SiconosVector>(*(ds.r()));

  state_x_.resize(2);
  if (ds.x()) state_x_[0] = std::make_shared<siconos::algebra::SiconosVector>(*(ds.x()));
  if (ds.rhs()) state_x_[1] = std::make_shared<siconos::algebra::SiconosVector>(*(ds.rhs()));
  jacobianRhsOver_x_ = ds.jacobianRhsOver_x();

  xMemory_ = ds.xMemory();
  stepsInMemory_ = ds.stepsInMemory();
}

void siconos::modeling::DynamicalSystem::resetToInitialState() {
  if (hasX0()) {
    *(state_x_[0]) = x0();  // Copy
  } else
    THROW_EXCEPTION(
        "siconos::modeling::DynamicalSystem::resetToInitialState() - initial state is not "
        "set");
}

void siconos::modeling::DynamicalSystem::update(double time) {
  computeRhs(time);
  computeJacobianRhsOver_x(time);
}

void siconos::modeling::DynamicalSystem::setX0(const siconos::algebra::SiconosVector& newValue,
                                               siconos::algebra::CopyTag tag) {
  if (newValue.size() != x_size_)
    throw std::invalid_argument("setX0(copy): input vector has wrong size");

  // Deep copy into Owned storage
  x0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(newValue);
}

void siconos::modeling::DynamicalSystem::setX0(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue, siconos::algebra::AliasTag tag) {
  if (newValue.size() != x_size_)
    throw std::invalid_argument("setX0(alias): input vector has wrong size");

  // Aliasing: shared_ptr to a Map (view over external memory)
  x0_storage_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void siconos::modeling::DynamicalSystem::initMemory(unsigned int steps) {
  DEBUG_BEGIN("void siconos::modeling::DynamicalSystem::initMemory(unsigned int steps)\n");
  if (steps > xMemory_.size()) {
    if (steps == 0)
      std::cout << "Warning : initMemory with size equal to zero" << std::endl;
    else {
      stepsInMemory_ = steps;
      xMemory_.setMemorySize(steps, x_size_);
    }
  }
  DEBUG_EXPR(siconos::algebra::print(xMemory_););

  DEBUG_END("void siconos::modeling::DynamicalSystem::initMemory(unsigned int steps)\n");
}
