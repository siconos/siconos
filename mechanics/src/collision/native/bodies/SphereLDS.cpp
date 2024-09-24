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

#include "SphereLDS.hpp"

#include <cmath>    // fmod
#include <numbers>  // pi

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace {  // Anonymous

static void normalize(Eigen::Ref<siconos::algebra::MapVectorType> q, unsigned int i) {
  q(i) = fmod(q(i), 2.0 * std::numbers::pi);

  assert(fabs(q(i)) - std::numeric_limits<double>::epsilon() >= 0.);
  assert(fabs(q(i)) < 2.0 * std::numbers::pi);
}

}  // namespace

siconos::collision::native::bodies::SphereLDS::SphereLDS(
    double r, double m, std::shared_ptr<siconos::algebra::SiconosVector> qinit,
    std::shared_ptr<siconos::algebra::SiconosVector> vinit)
    : siconos::modeling::LagrangianDS{qinit, vinit}, radius{r}, massValue{m} {
  normalize(*q(), 3);
  normalize(*q(), 4);
  normalize(*q(), 5);
  ndof_ = 6;

  assert(qinit->size() == ndof_);
  assert(vinit->size() == ndof_);

  I = massValue * radius * radius * 2. / 5.;

  setComputeMassFunction([this](Eigen::Ref<siconos::algebra::MapVectorType> pos, double time,
                                Eigen::Ref<siconos::algebra::MapType> mass_result) {
    normalize(pos, 3);
    normalize(pos, 4);
    normalize(pos, 5);

    // // SS: Forcing modification of qold, is this necessary?
    // if (qMemory().nbVectorsInMemory() >= 1)
    // {
    //   SiconosVector& qold = qMemory().getsiconos::algebra::SiconosVector(0);
    //   normalize(qold, 3);
    //   normalize(qold, 4);
    //   normalize(qold, 5);
    // }

    double theta = pos(3);

    assert(fabs(theta) - std::numeric_limits<double>::epsilon() >= 0.);

    mass_result(4, 5) = mass_result(5, 4) = I * cos(theta);
  });

  mass_view_->setZero();
  (*mass_view_)(0, 0) = (*mass_view_)(1, 1) = (*mass_view_)(2, 2) = massValue;
  (*mass_view_)(3, 3) = (*mass_view_)(4, 4) = (*mass_view_)(5, 5) = I;
  computeMass(*q());

  _jacobianFGyrq = std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
  _jacobianFGyrqDot = std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);

  _fGyr = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  _fGyr->zero();

  computeFGyr();
}

void siconos::collision::native::bodies::SphereLDS::computeFGyr() {
  siconos::collision::native::bodies::SphereLDS::computeFGyr(q(), velocity());
}

void siconos::collision::native::bodies::SphereLDS::computeFGyr(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  assert(q->size() == 6);
  assert(v->size() == 6);

  //  normalize(q,3);
  // normalize(q,4);
  // normalize(q,5);

  double theta = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot = v->getValue(4);
  double psidot = v->getValue(5);

  double sintheta = sin(theta);

  (*_fGyr)(0) = (*_fGyr)(1) = (*_fGyr)(2) = 0;

  (*_fGyr)(3) = I * psidot * phidot * sintheta;
  (*_fGyr)(4) = -I * psidot * thetadot * sintheta;
  (*_fGyr)(5) = -I * phidot * thetadot * sintheta;
}

void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrq() {
  siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrq(q(), velocity());
}
void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrqDot() {
  siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrqDot(q(), velocity());
}

void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrq(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  double theta = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot = v->getValue(4);
  double psidot = v->getValue(5);

  double costheta = cos(theta);

  _jacobianFGyrq->zero();

  (*_jacobianFGyrq)(3, 3) = -I * psidot * phidot * costheta;
  (*_jacobianFGyrq)(4, 3) = I * psidot * thetadot * costheta;
  (*_jacobianFGyrq)(5, 3) = I * psidot * thetadot * costheta;
}
void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrqDot(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  double theta = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot = v->getValue(4);
  double psidot = v->getValue(5);

  double sintheta = sin(theta);

  _jacobianFGyrqDot->zero();

  (*_jacobianFGyrqDot)(3, 3) = 0;
  (*_jacobianFGyrqDot)(3, 4) = I * psidot * sintheta;
  (*_jacobianFGyrqDot)(3, 5) = I * phidot * sintheta;

  (*_jacobianFGyrqDot)(4, 3) = -I * psidot * sintheta;
  (*_jacobianFGyrqDot)(4, 4) = 0;
  (*_jacobianFGyrqDot)(4, 5) = -I * thetadot * sintheta;

  (*_jacobianFGyrqDot)(5, 3) = -I * phidot * sintheta;
  (*_jacobianFGyrqDot)(5, 4) = -I * thetadot * sintheta;
  (*_jacobianFGyrqDot)(5, 5) = 0;
}

double siconos::collision::native::bodies::SphereLDS::getQ(unsigned int pos) {
  assert(pos < ndof_);
  return (*_q[0])(pos);
};

double siconos::collision::native::bodies::SphereLDS::getVelocity(unsigned int pos) {
  assert(pos < ndof_);
  return (*_q[1])(pos);
};
