/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "Circle.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::collision::native::bodies::Circle::Circle(
    double r, double m, const siconos::algebra::SiconosVector3& qinit,
    const siconos::algebra::SiconosVector3& vinit)
    : CircularDS(r, m, qinit, vinit) {
  mass_storage_ = std::make_unique<siconos::algebra::SiconosMatrix>(3, 3);
  use_mass([&](auto& M) {
    M.setZero();
    M(0, 0) = massValue_;
    M(1, 1) = massValue_;
    M(2, 2) = massValue_ * radius_ * radius_;
  });
  hasMass_ = true;
  hasConstantMass_ = true;
  computemass_ = nullptr;
}
