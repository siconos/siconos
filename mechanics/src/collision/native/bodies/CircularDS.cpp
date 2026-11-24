// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "CircularDS.hpp"

#include "SiconosVector.hpp"
#include "StorageTools.hpp"

siconos::collision::native::bodies::CircularDS::CircularDS(
    double r, double m, const siconos::algebra::SiconosVector3& qinit,
    const siconos::algebra::SiconosVector3& vinit)
    : siconos::modeling::LagrangianDS{qinit, vinit, siconos::algebra::copy_t},
      radius_{r},
      massValue_{m} {
  ndof_ = 3;
  // Nothing is done regarding the mass matrix ...
}

double siconos::collision::native::bodies::CircularDS::getQ(siconos::algebra::Index pos) {
  assert(pos < ndof_);
  return (*state_q_[0])(pos);
};

double siconos::collision::native::bodies::CircularDS::getVelocity(
    siconos::algebra::Index pos) {
  assert(pos < ndof_);
  return (*state_q_[1])(pos);
};
