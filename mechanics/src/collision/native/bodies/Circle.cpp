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
