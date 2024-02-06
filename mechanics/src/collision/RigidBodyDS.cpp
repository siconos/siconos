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

#include "RigidBodyDS.hpp"

#include "SiconosContactor.hpp"
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"

siconos::collision::RigidBodyDS::RigidBodyDS(
    std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity, double mass,
    std::shared_ptr<siconos::algebra::SiconosMatrix> inertia)
    : siconos::modeling::NewtonEulerDS(position, velocity, mass, inertia),
      _contactors(std::make_shared<siconos::collision::SiconosContactorSet>()) {}

void siconos::collision::RigidBodyDS::acceptSP(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
