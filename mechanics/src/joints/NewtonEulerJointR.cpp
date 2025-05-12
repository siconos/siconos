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
/*! \file NewtonEulerR.cpp

*/

#include "NewtonEulerJointR.hpp"

#include "BlockVector.hpp"
#include "SiconosVector.hpp"

void siconos::joints::NewtonEulerJointR::setPoint(
    size_t index, const Eigen::Ref<siconos::algebra::SiconosVector3>& point) {
  // Should we copy point into points_[index] or shared memory?
  // - if shared memory -> see for example the handling of q0 in LagrangianDS
  // - if copy, see below. That's the way we choose for the moment.
  assert(index < points_.size());
  points_[index] = point;
}

void siconos::joints::NewtonEulerJointR::setAxis(
    size_t index, const Eigen::Ref<siconos::algebra::SiconosVector3>& axis) {
  assert(index < axes_.size());
  axes_[index] = axis;
  axes_[index].normalize();
}

siconos::algebra::SiconosVector3 siconos::joints::NewtonEulerJointR::projectVectorDoF(
    const siconos::algebra::SiconosVector& v, const siconos::algebra::BlockVector& q0,
    int axis, bool absoluteRef) {
  auto ax = normalDoF(q0, axis, absoluteRef);

  double L = v(0) * ax(0) + v(1) * ax(1) + v(2) * ax(2);
  return siconos::algebra::SiconosVector3{ax(0) * L, ax(1) * L, ax(2) * L};  // RVO
}
