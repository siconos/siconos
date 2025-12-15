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
    double radius, double mass, const siconos::algebra::SiconosVector6& qinit,
    const siconos::algebra::SiconosVector6& vinit)
    : siconos::modeling::LagrangianDS{qinit, vinit, siconos::algebra::copy_t},
      radius_{radius},
      massValue_{mass} {
  normalize(*q(), 3);
  normalize(*q(), 4);
  normalize(*q(), 5);

  inertia_ = massValue_ * radius_ * radius_ * 2. / 5.;

  // setComputeMassFunction([this](const Eigen::Ref<const siconos::algebra::MapVectorType> pos,
  //                               Eigen::Ref<siconos::algebra::MapType> mass_result) {
  //   normalize(pos, 3);
  //   normalize(pos, 4);
  //   normalize(pos, 5);

  //   // // SS: Forcing modification of qold, is this necessary?
  //   // if (qMemory().nbVectorsInMemory() >= 1)
  //   // {
  //   //   SiconosVector& qold = qMemory().getsiconos::algebra::SiconosVector(0);
  //   //   normalize(qold, 3);
  //   //   normalize(qold, 4);
  //   //   normalize(qold, 5);
  //   // }

  //   double theta = pos(3);

  //   assert(fabs(theta) - std::numeric_limits<double>::epsilon() >= 0.);

  //   mass_result(4, 5) = mass_result(5, 4) = inertia_ * cos(theta);
  // });
  mass_storage_ = std::make_unique<siconos::algebra::SiconosDenseMatrix>(6, 6);
  use_mass([&](auto& M) {
    M.setZero();
    M(0, 0) = M(1, 1) = M(2, 2) = massValue_;
    M(3, 3) = M(4, 4) = M(5, 5) = inertia_;
  });
  hasMass_ = true;
  hasConstantMass_ = true;
  computemass_ = nullptr;

  computeMass(*q());

  // Set function to compute gyroscopic forces
  setComputeFgyrFunction(
      [this](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
             Eigen::Ref<siconos::algebra::MapVectorType> result) {
        assert(q.size() == 6);
        assert(velocity.size() == 6);

        double sintheta = sin(q(3));

        result(0) = result(1) = result(2) = 0.;

        result(3) = inertia_ * velocity(5) * velocity(4) * sintheta;
        result(4) = -inertia_ * velocity(5) * velocity(3) * sintheta;
        result(5) = -inertia_ * velocity(4) * velocity(3) * sintheta;
      });

  // and the jacobian of fgyr
  setComputeJacobianFgyrOver_qFunction(
      [this](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
             Eigen::Ref<siconos::algebra::MapVectorType> result) {
        double costheta = cos(q(3));

        result.setZero();

        result(3, 3) = -inertia_ * velocity(5) * velocity(4) * costheta;
        result(4, 3) = inertia_ * velocity(5) * velocity(3) * costheta;
        result(5, 3) = inertia_ * velocity(5) * velocity(3) * costheta;
      });

  setComputeJacobianFgyrOver_velocityFunction(
      [this](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
             Eigen::Ref<siconos::algebra::MapVectorType> result) {
        double sintheta = sin(q(3));
        result.setZero();

        result(3, 3) = 0;
        result(3, 4) = inertia_ * velocity(5) * sintheta;
        result(3, 5) = inertia_ * velocity(4) * sintheta;

        result(4, 3) = -inertia_ * velocity(5) * sintheta;
        result(4, 4) = 0;
        result(4, 5) = -inertia_ * velocity(3) * sintheta;

        result(5, 3) = -inertia_ * velocity(4) * sintheta;
        result(5, 4) = -inertia_ * velocity(3) * sintheta;
        result(5, 5) = 0;
      });

  computeFgyr(velocity_read(), q_read());
}
double siconos::collision::native::bodies::SphereLDS::getQ(unsigned int pos) {
  return (*state_q_[0])(pos);
};

double siconos::collision::native::bodies::SphereLDS::getVelocity(unsigned int pos) {
  return (*state_q_[1])(pos);
};

// FP : we need to override computeMass because it also update q which is supposed to be a
// constant input ...
// This is kind of ugly ...
void siconos::collision::native::bodies::SphereLDS::computeMass(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position) {
  state_q_[0]->tail(3) = position.tail(3);
  normalize(*state_q_[0], 3);
  normalize(*state_q_[0], 4);
  normalize(*state_q_[0], 5);

  // // SS: Forcing modification of qold, is this necessary?
  // if (qMemory().nbVectorsInMemory() >= 1)
  // {
  //   SiconosVector& qold = qMemory().getsiconos::algebra::SiconosVector(0);
  //   normalize(qold, 3);
  //   normalize(qold, 4);
  //   normalize(qold, 5);
  // }

  double theta = (*state_q_[0])(3);

  assert(fabs(theta) - std::numeric_limits<double>::epsilon() >= 0.);
  use_mass([&](auto& M) {
    M.setZero();
    M(4, 5) = M(5, 4) = inertia_ * cos(theta);
  });
}
