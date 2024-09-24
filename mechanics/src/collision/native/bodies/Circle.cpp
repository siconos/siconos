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
    double r, double m, std::shared_ptr<siconos::algebra::SiconosVector> qinit,
    std::shared_ptr<siconos::algebra::SiconosVector> vinit)
    : CircularDS(r, m, qinit, vinit) {
  mass_internal_storage_ = std::make_unique<std::vector<double>>(ndof_ * ndof_);
  mass_view_ = std::make_shared<MapType>(mass_internal_storage_->data(), ndof_, ndof_);
  hasConstantMass_ = true;
  hasMass_ = true;
  computemass_ = nullptr;

  mass_view_->setZero();
  (*mass_view_)(0, 0) = (*mass_view_)(1, 1) = massValue;
  (*mass_view_)(2, 2) = massValue * radius * radius;
}
