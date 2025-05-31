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

#include "DiskMovingPlanR.hpp"

#include <cmath>

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::collision::native::bodies::DiskMovingPlanR::DiskMovingPlanR(double radius)
    : siconos::modeling::LagrangianRheonomousR{}, radius_{radius} {}

void siconos::collision::native::bodies::DiskMovingPlanR::init(double time) {
  if (time != _time) {
    _time = time;
    computeA(time);
    computeB(time);
    computeC(time);
    computeADot(time);
    computeBDot(time);
    computeCDot(time);

    sqrA2pB2_ = hypot(A_, B_);
    AAdot_ = A_ * Adot_;
    BBdot_ = B_ * Bdot_;
    cubsqrA2pB2_ = sqrA2pB2_ * sqrA2pB2_ * sqrA2pB2_;
  }
}

double siconos::collision::native::bodies::DiskMovingPlanR::distance(double x, double y,
                                                                     double rad) {
  return (fabs(A_ * x + B_ * y + C_) / sqrA2pB2_ - rad);
}

/* Called compute h, but only the gap function is needed! */
void siconos::collision::native::bodies::DiskMovingPlanR::computeh(
    const siconos::algebra::BlockVector &q, double time,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  init(time);

  double q_0 = q(0);
  double q_1 = q(1);

  y(0) = distance(q_0, q_1, radius_);
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector &q, double time) {
  init(time);

  double x = q(0);
  double y = q(1);

  double D1 = A_ * x + B_ * y + C_;
  double signD1 = copysign(1, D1);

  jacobianhOver_q_view_->setValue(0, 0, A_ * signD1 / sqrA2pB2_);
  jacobianhOver_q_view_->setValue(1, 0, -B_ * signD1 / sqrA2pB2_);
  jacobianhOver_q_view_->setValue(0, 1, B_ * signD1 / sqrA2pB2_);
  jacobianhOver_q_view_->setValue(1, 1, A_ * signD1 / sqrA2pB2_);
  jacobianhOver_q_view_->setValue(0, 2, 0);
  jacobianhOver_q_view_->setValue(1, 2, -radius_);
}

void siconos::collision::native::bodies::DiskMovingPlanR::computehdot(
    const siconos::algebra::BlockVector &q, double time) {
  init(time);

  double x = q(0);
  double y = q(1);

  double D1 = A_ * x + B_ * y + C_;
  double signD1 = copysign(1, D1);
  (*hdot_)(0) = (-AAdot_ - BBdot_) * fabs(D1) / cubsqrA2pB2_ +
                (Adot_ * x + Bdot_ * y + Cdot_) * signD1 / sqrA2pB2_;
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeAFunction(
    const siconos::modeling::func_prototypes::FunctionS_S &fct) {
  computeA_ = fct;
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeBFunction(
    const siconos::modeling::func_prototypes::FunctionS_S &fct) {
  computeB_ = fct;
}
void siconos::collision::native::bodies::DiskMovingPlanR::setComputeCFunction(
    const siconos::modeling::func_prototypes::FunctionS_S &fct) {
  computeC_ = fct;
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeAdotFunction(
    const siconos::modeling::func_prototypes::FunctionS_S &fct) {
  computeAdot_ = fct;
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeBdotFunction(
    const siconos::modeling::func_prototypes::FunctionS_S &fct) {
  computeBdot_ = fct;
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeCdotFunction(
    const siconos::modeling::func_prototypes::FunctionS_S &fct) {
  computeCdot_ = fct;
}

bool siconos::collision::native::bodies::DiskMovingPlanR::equal(
    const siconos::modeling::func_prototypes::FunctionS_S &pA,
    const siconos::modeling::func_prototypes::FunctionS_S &pB,
    const siconos::modeling::func_prototypes::FunctionS_S &pC, double pr) const {
  // Note FP How can we compare user-defined functions ???
  // And why ... ?
  return false;  // (computeC_ == pA && computeB_ == pB && computeC_ == pC && radius_ == pr);
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeA(double t) {
  if (computeA_)
    A_ = computeA_(t);
  else
    A_ = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeB(double t) {
  if (computeB_)
    B_ = computeB_(t);
  else
    B_ = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeC(double t) {
  if (computeC_)
    C_ = computeC_(t);
  else
    C_ = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeADot(double t) {
  if (computeAdot_)
    Adot_ = computeAdot_(t);
  else
    Adot_ = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeBDot(double t) {
  if (computeBdot_)
    Bdot_ = computeBdot_(t);
  else
    Bdot_ = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeCDot(double t) {
  if (computeCdot_)
    Cdot_ = computeCdot_(t);
  else
    Cdot_ = 0.;
}
