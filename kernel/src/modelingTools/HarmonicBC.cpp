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
#include "SiconosVector.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
// #include "siconos_debug.h"

siconos::modeling::HarmonicBC::HarmonicBC(
    Indices newVelocityIndices, const Eigen::Ref<const siconos::algebra::SiconosVector>& newa,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& newb,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& omega,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& phi)
    : BoundaryCondition(newVelocityIndices),
      a_vec_{newa},
      b_vec_{newb},
      omega_vec_{omega},
      phi_vec_{phi},
      has_vector_coeffs_{true} {
  siconos::algebra::SiconosVector::Index bc_size = newVelocityIndices.size();
  assert(newa.size() == bc_size);
  assert(newb.size() == bc_size);
  assert(omega.size() == bc_size);
  assert(phi.size() == bc_size);
};

void siconos::modeling::HarmonicBC::computePrescribedVelocity(double time) {
  if (!has_vector_coeffs_) {
    auto val = aCoeff_ + bCoeff_ * cos(omega_ * time + phi_);
    prescribedVelocity_.setConstant(val);
  }

  else {
    // .array() --> component wise operations
    prescribedVelocity_ =
        a_vec_.array() + b_vec_.array() * (omega_vec_.array() * time + phi_vec_.array()).cos();
  }
}
