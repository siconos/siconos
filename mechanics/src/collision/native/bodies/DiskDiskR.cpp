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

#include "DiskDiskR.hpp"

#include <cmath>  // for hypot

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::collision::native::bodies::DiskDiskR::DiskDiskR(double r, double rr)
    : CircularR(r, rr) {
  r1pr2 = r + rr;
}

double siconos::collision::native::bodies::DiskDiskR::distance(double x1, double y1, double r1,
                                                               double x2, double y2,
                                                               double r2) {
  return (hypot(x1 - x2, y1 - y2) - r1pr2);
}

void siconos::collision::native::bodies::DiskDiskR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  double q_0 = q(0);
  double q_1 = q(1);
  double q_3 = q(3);
  double q_4 = q(4);

  y(0) = distance(q_0, q_1, _r1, q_3, q_4, _r2);
}

void siconos::collision::native::bodies::DiskDiskR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& q) {
  assert(jacobianhOver_q_view_);

  double x1 = q(0);
  double y1 = q(1);
  double x2 = q(3);
  double y2 = q(4);

  double dx = x2 - x1;
  double dy = y2 - y1;

  double d = hypot(dx, dy);

  double dxsd = dx / d;
  double dysd = dy / d;

  /*
  [-dx  -dy       dx   dy     ]
  [---  ---   0   --   --   0 ]
  [ d    d        d    d      ]
  [                           ]
  [dy   -dx       -dy  dx     ]
  [--   ---  -r1  ---  --  -r2]
  [d     d         d   d      ]
  */

  jacobianhOver_q_view_->setValue(0, 0, -dxsd);
  jacobianhOver_q_view_->setValue(1, 0, dysd);
  jacobianhOver_q_view_->setValue(0, 1, -dysd);
  jacobianhOver_q_view_->setValue(1, 1, -dxsd);
  jacobianhOver_q_view_->setValue(0, 2, 0.);
  jacobianhOver_q_view_->setValue(1, 2, -_r1);
  jacobianhOver_q_view_->setValue(0, 3, dxsd);
  jacobianhOver_q_view_->setValue(1, 3, -dysd);
  jacobianhOver_q_view_->setValue(0, 4, dysd);
  jacobianhOver_q_view_->setValue(1, 4, dxsd);
  jacobianhOver_q_view_->setValue(0, 5, 0.);
  jacobianhOver_q_view_->setValue(1, 5, -_r2);
}
