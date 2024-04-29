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

#include "CircleCircleR.hpp"

#include <cmath>

#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"

double siconos::collision::native::bodies::CircleCircleR::distance(double x1, double y1,
                                                                   double r1, double x2,
                                                                   double y2, double r2)
{
  return (fabs(r1 - r2) - hypot(x1 - x2, y1 - y2));
}

void siconos::collision::native::bodies::CircleCircleR::computeh(
    const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z,
    siconos::algebra::SiconosVector& y)
{
  double q_0 = q(0);
  double q_1 = q(1);
  double q_3 = q(3);
  double q_4 = q(4);

  y(0) = distance(q_0, q_1, _r1, q_3, q_4, _r2);
}

void siconos::collision::native::bodies::CircleCircleR::computeJachq(
    const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z)
{
  double x1 = q(0);
  double y1 = q(1);
  double x2 = q(3);
  double y2 = q(4);

  double dx = x2 - x1;
  double dy = y2 - y1;

  double d = hypot(dx, dy);

  double dxsd = dx / d;
  double dysd = dy / d;

  _jachq->setValue(0, 0, dxsd);
  _jachq->setValue(1, 0, -dysd);
  _jachq->setValue(0, 1, dysd);
  _jachq->setValue(1, 1, dxsd);
  _jachq->setValue(0, 2, 0.);
  _jachq->setValue(1, 2, -_r1);
  _jachq->setValue(0, 3, -dxsd);
  _jachq->setValue(1, 3, dysd);
  _jachq->setValue(0, 4, -dysd);
  _jachq->setValue(1, 4, -dxsd);
  _jachq->setValue(0, 5, 0.);
  _jachq->setValue(1, 5, -_r2);
}
