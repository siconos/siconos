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

#include "BulletSiconosCommon.hpp"

#include <LinearMath/btVector3.h>

#include <boost/math/quaternion.hpp>

#include "SiconosVector.hpp"

void siconos::collision::bullet::copyQuatPos(const btVector3& from,
                                             boost::math::quaternion<double>& to) {
  to = boost::math::quaternion<double>{0, from.x(), from.y(), from.z()};
}

void siconos::collision::bullet::copyBtVector3(const btVector3& from,
                                               siconos::algebra::SiconosVector3& to) {
  to(0) = from.x();
  to(1) = from.y();
  to(2) = from.z();
}

void siconos::collision::bullet::copyQuatPos2d(const btVector3& from,
                                               boost::math::quaternion<double>& to) {
  to = boost::math::quaternion<double>{0, from.x(), from.y(), 0.0};
}

void siconos::collision::bullet::copyBtVector32d(const btVector3& from,
                                                 siconos::algebra::SiconosVector2& to) {
  to(0) = from.x();
  to(1) = from.y();
}
