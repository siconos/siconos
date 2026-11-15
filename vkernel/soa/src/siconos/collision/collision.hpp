/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
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

#pragma once

#include "collision_head.hpp"

namespace siconos::collision {

/* 2D */
decltype(auto) distance_sq(match::vector auto& q1, match::vector auto& q2) {
  const auto dx = q2(0) - q1(0);
  const auto dy = q2(1) - q1(1);
  return dx * dx + dy * dy;
}

decltype(auto) distance(match::vector auto& q1, match::vector auto& q2) {
  return sqrt(distance_sq(q1, q2));
}

}  // namespace siconos::collision
