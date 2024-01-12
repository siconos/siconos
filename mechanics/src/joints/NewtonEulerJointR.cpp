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
    unsigned int index, std::shared_ptr<siconos::algebra::SiconosVector> point) {
  _points[index] = point;
}

auto siconos::joints::NewtonEulerJointR::point(unsigned int index) { return _points[index]; }

void siconos::joints::NewtonEulerJointR::setAxis(
    unsigned int index, std::shared_ptr<siconos::algebra::SiconosVector> axis) {
  _axes[index] = axis;
}

auto siconos::joints::NewtonEulerJointR::axis(unsigned int index) { return _axes[index]; }

auto& siconos::joints::NewtonEulerJointR::axes() { return _axes; }

void siconos::joints::NewtonEulerJointR::projectVectorDoF(
    const siconos::algebra::SiconosVector& v, const siconos::algebra::BlockVector& q0,
    siconos::algebra::SiconosVector& ans, int axis, bool absoluteRef) {
  siconos::algebra::SiconosVector ax(3);
  normalDoF(ax, q0, axis, absoluteRef);

  double L = v(0) * ax(0) + v(1) * ax(1) + v(2) * ax(2);
  ans(0) = ax(0) * L;
  ans(1) = ax(1) * L;
  ans(2) = ax(2) * L;
}

std::shared_ptr<siconos::algebra::SiconosVector>
siconos::joints::NewtonEulerJointR::projectVectorDoF(const siconos::algebra::SiconosVector& v,
                                                     const siconos::algebra::BlockVector& q0,
                                                     int axis, bool absoluteRef) {
  auto ans(std::make_shared<siconos::algebra::SiconosVector>(3));
  projectVectorDoF(v, q0, *ans, axis, absoluteRef);
  return ans;
}

void siconos::joints::NewtonEulerJointR::normalDoF(siconos::algebra::SiconosVector& ans,
                                                   const siconos::algebra::BlockVector& q0,
                                                   int axis, bool absoluteRef) {
  _normalDoF(ans, q0, axis, absoluteRef);
}

std::shared_ptr<siconos::algebra::SiconosVector> siconos::joints::NewtonEulerJointR::normalDoF(
    const siconos::algebra::BlockVector& q0, int axis, bool absoluteRef) {

  auto ans = std::make_shared<siconos::algebra::SiconosVector>(3);
  _normalDoF(*ans, q0, axis, absoluteRef);
  return ans;
}
