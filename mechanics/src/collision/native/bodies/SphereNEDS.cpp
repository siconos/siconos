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

#include "SphereNEDS.hpp"

#include <cmath>

#include "SiconosVector.hpp"

siconos::collision::native::bodies::SphereNEDS::SphereNEDS(
    double r, double m, Eigen::Ref<siconos::algebra::SiconosMatrix> inertia,
    Eigen::Ref<siconos::algebra::SiconosVector> qinit,
    Eigen::Ref<siconos::algebra::SiconosVector> vinit)
    : siconos::modeling::NewtonEulerDS{qinit, vinit, m, inertia}, radius{r} {
  // note : ndof_ = 3 in NewtonEuleurDS ? (=> ndof_ = 6 ?)

  assert(qinit.size() == qDim_);
  assert(vinit.size() == 6);  // == ndof_
}

double siconos::collision::native::bodies::SphereNEDS::getQ(siconos::algebra::Index pos) {
  assert(pos < 7);
  return ((*state_q_)(pos));
};

double siconos::collision::native::bodies::SphereNEDS::getTwist(siconos::algebra::Index pos) {
  assert(pos < 6);
  return ((*twist_)(pos));
};
