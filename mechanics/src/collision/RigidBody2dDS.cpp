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
