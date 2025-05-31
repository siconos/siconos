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

#include "HarmonicBC.hpp"

#include <math.h>  // for cos

#include "BoundaryCondition.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"

siconos::modeling::HarmonicBC::HarmonicBC(Indices&& newVelocityIndices,
                                          Eigen::Ref<siconos::algebra::SiconosVector> newa,
                                          Eigen::Ref<siconos::algebra::SiconosVector> newb,
                                          Eigen::Ref<siconos::algebra::SiconosVector> omega,
                                          Eigen::Ref<siconos::algebra::SiconosVector> phi)
    : BoundaryCondition(std::move(newVelocityIndices)) {
  siconos::algebra::SiconosVector::Index bc_size = newVelocityIndices.size();
  assert(newa.size() == bc_size);
  assert(newb.size() == bc_size);
  assert(omega.size() == bc_size);
  assert(phi.size() == bc_size);

  a_view_ = std::make_unique<siconos::algebra::MapVectorType>(newa.data(), bc_size);
  b_view_ = std::make_unique<siconos::algebra::MapVectorType>(newb.data(), bc_size);
  omega_view_ = std::make_unique<siconos::algebra::MapVectorType>(omega.data(), bc_size);
  phi_view_ = std::make_unique<siconos::algebra::MapVectorType>(phi.data(), bc_size);
};

void siconos::modeling::HarmonicBC::computePrescribedVelocity(double time) {
  assert(prescribedVelocity_);  // Allocation must have been done during construction.

  if (!a_view_) {
    auto val = aCoeff_ + bCoeff_ * cos(omega_ * time + phi_);
    prescribedVelocity_->setConstant(val);
  }

  else {
    auto a = *a_view_;
    auto b = *b_view_;
    auto omega = *omega_view_;
    auto phi = *phi_view_;

    // Calculer le résultat de manière optimale
    *prescribedVelocity_ = a.array() + b.array() * (omega.array() * time + phi.array()).cos();
  }
}
