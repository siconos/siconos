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
/*! \file KneeJointR.cpp

*/

#include "KneeJointR.hpp"

#include <boost/math/quaternion.hpp>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "NewtonEulerDS.hpp"
#include "RotationQuaternion.hpp"  // rotquat, posquat ...
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #include <cfloat>
// #include <iostream>

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::joints::KneeJointR::KneeJointR(const Eigen::Ref<siconos::algebra::SiconosVector3>& P,
                                        bool absoluteRef,
                                        std::shared_ptr<siconos::modeling::NewtonEulerDS> d1,
                                        std::shared_ptr<siconos::modeling::NewtonEulerDS> d2) {
  points_.emplace_back(P);
  setAbsolute(absoluteRef);
  if (d1) {
    if (d2) {
      setBasePositions(d1->q_read(), d2->q_read());
      checkInitPos(d1->q_read(), d2->q_read());
    } else {
      setBasePositions(d1->q_read());
      checkInitPos(d1->q_read());
    }
  }

  setComputeH_NE_dotFunction([this](const siconos::algebra::BlockVector& q,
                                    const siconos::algebra::BlockVector& qdot,
                                    Eigen::Ref<siconos::algebra::MapType> jacob) {
    jacob.setZero();

    const auto& qdot1 = *qdot.vector(0);

    if (qdot.numberOfBlocks() == 2) {
      const auto& qdot2 = *qdot.vector(1);
      knee::computeH_dot_for2DS(qdot1.tail(4), G1P0_, qdot2.tail(4), G2P0_, *H_NE_dot_);
    } else
      knee::computeH_dot_for1DS(qdot1.tail(4), G1P0_, *H_NE_dot_);
  });
}

siconos::joints::KneeJointR::KneeJointR() {
  points_.resize(1);
  points_[0].setZero();
}

void siconos::joints::KneeJointR::checkInitPos(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2) {
  // Ensure that each component of the compute hfunction is smaller than epsilon
  siconos::algebra::SiconosVector3 tmp;
  if (q2) {
    knee::hfunction(q1, G1P0_, q2, G2P0_, tmp);
  } else
    knee::hfunction(q1, G1P0_, std::nullopt, G2P0_, tmp);
  assert((tmp.array().abs() < std::numeric_limits<double>::epsilon()).all());

  //       [this, q1, q2]() {
  //     Eigen::Vector3d vec;
  //     knee::hfunction(q1, G1P0_, q2, G2P0_, vec);
  //     return (vec.array().abs() < std::numeric_limits<double>::epsilon()).all();
  //   }

  //       if(q2)
  // {  assert(([this, q1, q2]() {
  //     Eigen::Vector3d vec;
  //     knee::hfunction(q1, G1P0_, q2, G2P0_, vec);
  //     return (vec.array().abs() < std::numeric_limits<double>::epsilon()).all();
  //   })());
  // }else{
  //   assert(([this, q1]() {
  //      Eigen::Vector3d vec;
  //     knee::hfunction(q1, G1P0_, std::nullopt, G2P0_, vec);
  //     return (vec.array().abs() < std::numeric_limits<double>::epsilon()).all();
  //   })());}
}

void siconos::joints::KneeJointR::setBasePositions(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2) {
  // Assumes that points_[0] is properly set
  boost::math::quaternion<double> rot1{siconos::geometry::rotquat(q1)}, quatBuff, quatP0_abs;

  /** Computation of G1P0_ and G2P0_ */

  /* Calculate G1P0 and P0_abs */
  if (absoluteRef_) {
    quatP0_abs = siconos::geometry::posquat(points_[0]);

    /* Move to q1 frame by unapplying q1 frame translation/rotation */
    quatBuff = (1.0 / rot1) * (quatP0_abs - siconos::geometry::posquat(q1)) * rot1;
    G1P0_ << quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4();
  } else {
    G1P0_ = points_[0];

    /* Move to abs frame by applying q1 frame rotation/translation */
    quatP0_abs = (rot1 * siconos::geometry::posquat(points_[0]) / rot1) +
                 siconos::geometry::posquat(q1);
  }

  /* Calculate G2P0, or set it to P0_abs (i.e. G2=absolute frame) */
  if (q2) {
    auto rot2{siconos::geometry::rotquat(*q2)};

    /* Move to q2 frame by unapplying q2 frame translation/rotation */
    quatBuff = (1.0 / rot2) * (quatP0_abs - siconos::geometry::posquat(*q2)) * rot2;
    G2P0_ << quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4();
  } else {
    /* q2 frame = absolute frame */
    G2P0_ << quatP0_abs.R_component_2(), quatP0_abs.R_component_3(),
        quatP0_abs.R_component_4();
  }
}

void siconos::joints::KneeJointR::computeH_NE_(double time,
                                               siconos::modeling::Interaction& inter,
                                               const siconos::algebra::BlockVector& q0) {
  H_NE_view_->setZero();
  auto q1 = q0.vector(0);
  // Only the quaternion part of q is required to compute H (last 4 components)
  if (q0.numberOfBlocks() > 1) {
    auto q2 = q0.vector(1);
    knee::computeH_for_2DS(q1->tail(4), G1P0_, q2->tail(4), G2P0_, *H_NE_view_);
  } else
    knee::computeH_for_1DS(q1->tail(4), G1P0_, *H_NE_view_);
}

void siconos::joints::KneeJointR::computeh(const siconos::algebra::BlockVector& q0,
                                           Eigen::Ref<siconos::algebra::SiconosVector> y) {
  if (q0.numberOfBlocks() == 2) {
    knee::hfunction(*q0.vector(0), G1P0_, *q0.vector(1), G2P0_, y);
  } else {
    knee::hfunction(*q0.vector(0), G1P0_, std::nullopt, G2P0_, y);
  }
}

void siconos::joints::KneeJointR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  if (q2) {
    knee::hfunction(q1, G1P0_, q2, G2P0_, y);
  } else {
    knee::hfunction(q1, G1P0_, std::nullopt, G2P0_, y);
  }
}

// Stand-alone functions
void siconos::joints::knee::hfunction(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const siconos::algebra::SiconosVector3& coords1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2opt,
    const siconos::algebra::SiconosVector3& coords2,
    Eigen::Ref<siconos::algebra::SiconosVector> result) {
  auto t1 = q1(4) * q1(4);
  auto t2 = q1(3) * q1(3);
  auto t3 = q1(6) * q1(6);
  auto t4 = q1(5) * q1(5);

  if (q2opt) {
    const auto& q2 = *q2opt;
    auto t17 = q2(4) * q2(4);
    auto t18 = q2(3) * q2(3);
    auto t19 = q2(6) * q2(6);
    auto t20 = q2(5) * q2(5);
    result(0) = q1(0) - q2(0) + (t1 + t2 - t3 - t4) * coords1(0) +
                (q1(4) * 2.0 * q1(5) - q1(3) * 2.0 * q1(6)) * coords1(1) +
                (q1(4) * 2.0 * q1(6) + q1(3) * 2.0 * q1(5)) * coords1(2) -
                (t17 + t18 - t19 - t20) * coords2(0) -
                (q2(4) * 2.0 * q2(5) - 2.0 * q2(3) * q2(6)) * coords2(1) -
                (q2(4) * 2.0 * q2(6) + 2.0 * q2(3) * q2(5)) * coords2(2);
    result(1) = q1(1) - q2(1) + (q1(4) * 2.0 * q1(5) + q1(3) * 2.0 * q1(6)) * coords1(0) +
                (t4 - t3 + t2 - t1) * coords1(1) +
                (q1(5) * 2.0 * q1(6) - q1(3) * 2.0 * q1(4)) * coords1(2) -
                (q2(4) * 2.0 * q2(5) + 2.0 * q2(3) * q2(6)) * coords2(0) -
                (t20 - t19 + t18 - t17) * coords2(1) -
                (q2(5) * 2.0 * q2(6) - 2.0 * q2(3) * q2(4)) * coords2(2);

    result(2) = q1(2) - q2(2) + (q1(4) * 2.0 * q1(6) - q1(3) * 2.0 * q1(5)) * coords1(0) +
                (q1(5) * 2.0 * q1(6) + q1(3) * 2.0 * q1(4)) * coords1(1) +
                (t3 - t4 - t1 + t2) * coords1(2) -
                (q2(4) * 2.0 * q2(6) - 2.0 * q2(3) * q2(5)) * coords2(0) -
                (q2(5) * 2.0 * q2(6) + 2.0 * q2(3) * q2(4)) * coords2(1) -
                (t19 - t20 - t17 + t18) * coords2(2);

  } else {
    result(0) = q1(0) + (t1 + t2 - t3 - t4) * coords1(0) +
                (q1(4) * 2.0 * q1(5) - q1(3) * 2.0 * q1(6)) * coords1(1) +
                (q1(4) * 2.0 * q1(6) + q1(3) * 2.0 * q1(5)) * coords1(2) - coords2(0);
    result(1) = q1(1) + (q1(4) * 2.0 * q1(5) + q1(3) * 2.0 * q1(6)) * coords1(0) +
                (t4 - t3 + t2 - t1) * coords1(1) +
                (q1(5) * 2.0 * q1(6) - q1(3) * 2.0 * q1(4)) * coords1(2) - coords2(1);
    result(2) = q1(2) + (q1(4) * 2.0 * q1(6) - q1(3) * 2.0 * q1(5)) * coords1(0) +
                (q1(5) * 2.0 * q1(6) + q1(3) * 2.0 * q1(4)) * coords1(1) +
                (t3 - t4 - t1 + t2) * coords1(2) - coords2(2);
  }
}

void siconos::joints::knee::computeH_for_2DS(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp1,
    const siconos::algebra::SiconosVector3& coords1,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp2,
    const siconos::algebra::SiconosVector3& coords2,
    Eigen::Ref<siconos::algebra::MapType> result) {
  assert(result.cols() == 14);
  assert(result.rows() == 3);
  assert(qp1.size() == 4);
  assert(qp2.size() == 4);

  result.setZero();
  result.setValue(0, 0, 1.0);
  result.setValue(
      0, 3, 2. * coords1(2) * qp1(2) - 2. * coords1(1) * qp1(3) + 2. * coords1(0) * qp1(0));
  result.setValue(
      0, 4, 2. * coords1(2) * qp1(3) + 2. * coords1(1) * qp1(2) + 2. * coords1(0) * qp1(1));
  result.setValue(
      0, 5, 2. * coords1(2) * qp1(0) + 2. * coords1(1) * qp1(1) - 2. * coords1(0) * qp1(2));
  result.setValue(
      0, 6, 2. * coords1(2) * qp1(1) - 2. * coords1(1) * qp1(0) - 2. * coords1(0) * qp1(3));
  result.setValue(0, 7, -1.0);
  result.setValue(
      0, 10, -2. * coords2(2) * qp2(2) + 2. * coords2(1) * qp2(3) - 2. * coords2(0) * qp2(0));
  result.setValue(
      0, 11, -2. * coords2(2) * qp2(3) - 2. * coords2(1) * qp2(2) - 2. * coords2(0) * qp2(1));
  result.setValue(
      0, 12, -2. * coords2(2) * qp2(0) - 2. * coords2(1) * qp2(1) + 2. * coords2(0) * qp2(2));
  result.setValue(
      0, 13, -2. * coords2(2) * qp2(1) + 2. * coords2(1) * qp2(0) + 2. * coords2(0) * qp2(3));
  result.setValue(1, 1, 1.0);
  result.setValue(
      1, 3, -2. * coords1(2) * qp1(1) + 2. * coords1(0) * qp1(3) + 2. * coords1(1) * qp1(0));
  result.setValue(
      1, 4, -2. * coords1(2) * qp1(0) + 2. * coords1(0) * qp1(2) - 2. * coords1(1) * qp1(1));
  result.setValue(
      1, 5, -2. * coords1(2) * qp1(3) + 2. * coords1(0) * qp1(1) + 2. * coords1(1) * qp1(2));
  result.setValue(
      1, 6, -2. * coords1(2) * qp1(2) + 2. * coords1(0) * qp1(0) - 2. * coords1(1) * qp1(3));
  result.setValue(1, 8, -1.0);
  result.setValue(
      1, 10, 2. * coords2(2) * qp2(1) - 2. * coords2(0) * qp2(3) - 2. * coords2(1) * qp2(0));
  result.setValue(
      1, 11, 2. * coords2(2) * qp2(0) - 2. * coords2(0) * qp2(2) + 2. * coords2(1) * qp2(1));
  result.setValue(
      1, 12, -2. * coords2(2) * qp2(3) - 2. * coords2(0) * qp2(1) - 2. * coords2(1) * qp2(2));
  result.setValue(
      1, 13, -2. * coords2(2) * qp2(2) - 2. * coords2(0) * qp2(0) + 2. * coords2(1) * qp2(3));
  result.setValue(2, 2, 1.0);
  result.setValue(
      2, 3, 2. * coords1(1) * qp1(1) - 2. * coords1(0) * qp1(2) + 2. * coords1(2) * qp1(0));
  result.setValue(
      2, 4, 2. * coords1(1) * qp1(0) + 2. * coords1(0) * qp1(3) - 2. * coords1(2) * qp1(1));
  result.setValue(
      2, 5, 2. * coords1(1) * qp1(3) - 2. * coords1(0) * qp1(0) - 2. * coords1(2) * qp1(2));
  result.setValue(
      2, 6, 2. * coords1(1) * qp1(2) + 2. * coords1(0) * qp1(1) + 2. * coords1(2) * qp1(3));
  result.setValue(2, 9, -1.0);
  result.setValue(
      2, 10, -2. * coords2(1) * qp2(1) + 2. * coords2(0) * qp2(2) - 2. * coords2(2) * qp2(0));
  result.setValue(
      2, 11, -2. * coords2(1) * qp2(0) - 2. * coords2(0) * qp2(3) + 2. * coords2(2) * qp2(1));
  result.setValue(
      2, 12, -2. * coords2(1) * qp2(3) + 2. * coords2(0) * qp2(0) + 2. * coords2(2) * qp2(2));
  result.setValue(
      2, 13, -2. * coords2(1) * qp2(2) - 2. * coords2(0) * qp2(1) - 2. * coords2(2) * qp2(3));
}

void siconos::joints::knee::computeH_for_1DS(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp1,
    const siconos::algebra::SiconosVector3& coords1,
    Eigen::Ref<siconos::algebra::MapType> result) {
  assert(result.cols() == 7);
  assert(result.rows() == 3);
  assert(qp1.size() == 4);
  result.setZero();

  result.setValue(0, 0, 1.);
  result.setValue(
      0, 3, 2.0 * coords1(2) * qp1(2) - 2.0 * coords1(1) * qp1(3) + 2.0 * coords1(0) * qp1(0));
  result.setValue(
      0, 4, 2.0 * coords1(2) * qp1(3) + 2.0 * coords1(1) * qp1(2) + 2.0 * coords1(0) * qp1(1));
  result.setValue(
      0, 5, 2.0 * coords1(2) * qp1(0) + 2.0 * coords1(1) * qp1(1) - 2.0 * coords1(0) * qp1(2));
  result.setValue(
      0, 6, 2.0 * coords1(2) * qp1(1) - 2.0 * coords1(1) * qp1(0) - 2.0 * coords1(0) * qp1(3));
  result.setValue(1, 1, 1.0);
  result.setValue(
      1, 3,
      -2.0 * coords1(2) * qp1(1) + 2.0 * coords1(0) * qp1(3) + 2.0 * coords1(1) * qp1(0));
  result.setValue(
      1, 4,
      -2.0 * coords1(2) * qp1(0) + 2.0 * coords1(0) * qp1(2) - 2.0 * coords1(1) * qp1(1));
  result.setValue(
      1, 5, 2.0 * coords1(2) * qp1(3) + 2.0 * coords1(0) * qp1(1) + 2.0 * coords1(1) * qp1(2));
  result.setValue(
      1, 6, 2.0 * coords1(2) * qp1(2) + 2.0 * coords1(0) * qp1(0) - 2.0 * coords1(1) * qp1(3));
  result.setValue(2, 0, 0.0);
  result.setValue(2, 1, 0.0);
  result.setValue(2, 2, 1.0);
  result.setValue(
      2, 3, 2.0 * coords1(1) * qp1(1) - 2.0 * coords1(0) * qp1(2) + 2.0 * coords1(2) * qp1(0));
  result.setValue(
      2, 4, 2.0 * coords1(1) * qp1(0) + 2.0 * coords1(0) * qp1(3) - 2.0 * coords1(2) * qp1(1));
  result.setValue(
      2, 5, 2.0 * coords1(1) * qp1(3) - 2.0 * coords1(0) * qp1(0) - 2.0 * coords1(2) * qp1(2));
  result.setValue(
      2, 6, 2.0 * coords1(1) * qp1(2) + 2.0 * coords1(0) * qp1(1) + 2.0 * coords1(2) * qp1(3));
}

void siconos::joints::knee::computeH_dot_for1DS(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qpdot1,
    const siconos::algebra::SiconosVector3& coords1,
    Eigen::Ref<siconos::algebra::MapType> result) {
  double t1 = coords1(1) * qpdot1(0);
  double t2 = coords1(0) * qpdot1(0);
  double t4 = -t1 + t2 + coords1(2) * qpdot1(2);
  double t7 = coords1(2) * qpdot1(0);
  double t8 = coords1(0) * qpdot1(1) + coords1(1) * qpdot1(2) + t7;
  double t11 = t7 - coords1(0) * qpdot1(2) + coords1(1) * qpdot1(1);
  double t13 = -t2 + coords1(2) * qpdot1(1) - t1;

  result.setZero();
  result.setValue(0, 3, 2.0 * t4);
  result.setValue(0, 4, 2.0 * t8);
  result.setValue(0, 5, 2.0 * t11);
  result.setValue(0, 6, 2.0 * t13);
  result.setValue(1, 3, -2.0 * t13);
  result.setValue(1, 4, -2.0 * t11);
  result.setValue(1, 5, 2.0 * t8);
  result.setValue(1, 6, 2.0 * t4);
  result.setValue(2, 3, 2.0 * t11);
  result.setValue(2, 4, -2.0 * t13);
  result.setValue(2, 5, -2.0 * t4);
  result.setValue(2, 6, 2.0 * t8);
}

void siconos::joints::knee::computeH_dot_for2DS(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qpdot1,
    const siconos::algebra::SiconosVector3& coords1,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qpdot2,
    const siconos::algebra::SiconosVector3& coords2,
    Eigen::Ref<siconos::algebra::MapType> result) {
  double t1 = coords1(1) * qpdot1(0);
  double t2 = coords1(0) * qpdot1(0);
  double t4 = -t1 + t2 + coords1(2) * qpdot1(2);
  double t7 = coords1(2) * qpdot1(0);
  double t8 = coords1(0) * qpdot1(1) + coords1(1) * qpdot1(2) + t7;
  double t11 = t7 - coords1(0) * qpdot1(2) + coords1(1) * qpdot1(1);
  double t13 = -t2 + coords1(2) * qpdot1(1) - t1;
  double t17 = -coords2(0) * qpdot2(0) + coords2(1) * qpdot2(3) - coords2(2) * qpdot2(2);
  double t21 = -coords2(0) * qpdot2(1) - coords2(1) * qpdot2(2) - coords2(2) * qpdot2(3);
  double t25 = -coords2(2) * qpdot2(0) + coords2(0) * qpdot2(2) - coords2(1) * qpdot2(1);
  double t29 = coords2(1) * qpdot2(0) - coords2(2) * qpdot2(1) + coords2(0) * qpdot2(3);
  result.setZero();
  result.setValue(0, 3, 2.0 * t4);
  result.setValue(0, 4, 2.0 * t8);
  result.setValue(0, 5, 2.0 * t11);
  result.setValue(0, 6, 2.0 * t13);
  result.setValue(0, 10, 2.0 * t17);
  result.setValue(0, 11, 2.0 * t21);
  result.setValue(0, 12, 2.0 * t25);
  result.setValue(0, 13, 2.0 * t29);
  result.setValue(1, 3, -2.0 * t13);
  result.setValue(1, 4, -2.0 * t11);
  result.setValue(1, 5, 2.0 * t8);
  result.setValue(1, 6, 2.0 * t4);
  result.setValue(1, 10, -2.0 * t29);
  result.setValue(1, 11, -2.0 * t25);
  result.setValue(1, 12, 2.0 * t21);
  result.setValue(1, 13, 2.0 * t17);
  result.setValue(2, 3, 2.0 * t11);
  result.setValue(2, 4, -2.0 * t13);
  result.setValue(2, 5, -2.0 * t4);
  result.setValue(2, 6, 2.0 * t8);
  result.setValue(2, 10, 2.0 * t25);
  result.setValue(2, 11, -2.0 * t29);
  result.setValue(2, 12, -2.0 * t17);
  result.setValue(2, 13, 2.0 * t21);
}
