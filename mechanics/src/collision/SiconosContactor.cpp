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

#include "SiconosContactor.hpp"

#include "SiconosVector.hpp"

siconos::collision::SiconosContactor::SiconosContactor(
    std::shared_ptr<SiconosShape> shape_in,
    std::shared_ptr<siconos::algebra::SiconosVector> offset_in, int collision_group_in)
    : shape(shape_in), offset(offset_in), collision_group(collision_group_in)
{
  // First strategy: fill offset with the identity if not provided
  // if(!offset)
  // {
  //   offset = std::make_shared<siconos::algebra::SiconosVector>(7);
  //   offset->setZero();
  //   (*offset)(3) = 1.0;
  // }

  // Second strategy: leave offset as a null pointer if identity is provided
  if (offset) {
    if ((*offset)(3) == 1.0 && (*offset)(0) == 0.0 && (*offset)(1) == 0.0 &&
        (*offset)(2) == 0.0 && (*offset)(4) == 0.0 && (*offset)(5) == 0.0 &&
        (*offset)(6) == 0.0) {
      offset.reset();
    }
  }
}

siconos::collision::SiconosContactor::SiconosContactor(
    std::shared_ptr<SiconosShape> shape_in,
    Eigen::Ref<siconos::algebra::SiconosVector> offset_in, int collision_group_in)
    : SiconosContactor(shape_in, std::make_shared<siconos::algebra::SiconosVector>(offset_in), collision_group_in) {}
