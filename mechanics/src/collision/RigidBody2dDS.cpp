// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "RigidBody2dDS.hpp"

#include "SiconosContactor.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"

siconos::collision::RigidBody2dDS::RigidBody2dDS(
    const siconos::algebra::SiconosVector3& position,
    const siconos::algebra::SiconosVector3& velocity, double mass, double inertia)
    : LagrangianLinearTIDS(position, velocity, siconos::algebra::copy_t),
      scalarMass_(mass),
      contactors_(std::make_shared<siconos::collision::SiconosContactorSet>()) {
  mass_storage_ = std::make_unique<siconos::algebra::SiconosDenseMatrix>(ndof_, ndof_);
  use_mass([&](auto& M) {
    M.setZero();
    M(0, 0) = mass;
    M(1, 1) = mass;
    M(2, 2) = inertia;
  });
  hasMass_ = true;
  hasConstantMass_ = true;
  computemass_ = nullptr;
}
