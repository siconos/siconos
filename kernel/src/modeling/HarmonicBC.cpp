// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
      has_vector_coeffs_{true},
      a_vec_{newa},
      b_vec_{newb},
      omega_vec_{omega},
      phi_vec_{phi} {
  [[maybe_unused]] siconos::algebra::SiconosVector::Index bc_size =
      siconos::algebra::to_index(newVelocityIndices.size());
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
