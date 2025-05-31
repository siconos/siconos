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

#include "CircularDS.hpp"

#include "SiconosVector.hpp"

siconos::collision::native::bodies::CircularDS::CircularDS(
    double r, double m, Eigen::Ref<siconos::algebra::SiconosVector> qinit,
    Eigen::Ref<siconos::algebra::SiconosVector> vinit)
    : siconos::modeling::LagrangianDS{qinit, vinit}, radius_{r}, massValue_{m} {
  ndof_ = 3;
  // Nothing is done regarding the mass matrix ...
}

double siconos::collision::native::bodies::CircularDS::getQ(unsigned int pos) {
  assert(pos < ndof_);
  return (*state_q_[0])(pos);
};

double siconos::collision::native::bodies::CircularDS::getVelocity(unsigned int pos) {
  assert(pos < ndof_);
  return (*state_q_[1])(pos);
};
