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
/*! \file PivotJointR.cpp

*/

#include "PivotJointR.hpp"

#include <RotationQuaternion.hpp>  // for posquat, rotquat ...
#include <boost/math/quaternion.hpp>
#include <numbers>  // for pi

#include "BlockVector.hpp"
#include "NewtonEulerDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "op3x3.h"  // for orthoBaseFromVector

/*
 * This file contains some code generated using sympy.  The following
 * is the necessary predule:
 *
 * from sympy import Symbol
 * import numpy as np
 *
 * A1 = np.array([0, Symbol('axis1_(0)'), Symbol('axis1_(1)'), Symbol('axis1_(2)')])
 * A2 = np.array([0, Symbol('axis2_(0)'), Symbol('axis2_(1)'), Symbol('axis2_(2)')])
 * q1 = np.array([Symbol('q10'), Symbol('q11'), Symbol('q12'), Symbol('q13')])
 * q2 = np.array([Symbol('q20'), Symbol('q21'), Symbol('q22'), Symbol('q23')])
 * cq2q10 = np.array([Symbol('cq2q_(0)'),Symbol('cq2q_(1)'),
 *                    Symbol('cq2q_(2)'),Symbol('cq2q_(3)')])
 *
 * qinv = lambda q: np.array([q[0],-q[1],-q[2],-q[3]])
 * qmul = lambda a,b: np.array([
 *          a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
 *          a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
 *          a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
 *          a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]])
 */

namespace {  // anonymous namespace for local impl
double piwrap(double x) {
  return std::fmod(x + 3 * std::numbers::pi, 2 * std::numbers::pi) - std::numbers::pi;
}
}  // namespace

siconos::joints::PivotJointR::PivotJointR() : KneeJointR{} {
  axes_.resize(1);
  axes_[0].setZero();
}

siconos::joints::PivotJointR::PivotJointR(
    const Eigen::Ref<siconos::algebra::SiconosVector3>& P,
    const Eigen::Ref<siconos::algebra::SiconosVector3>& A, bool absoluteRef,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> d1,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> d2)
    : KneeJointR{P, absoluteRef, d1, d2} {
  axes_.emplace_back(A);
  axes_[0].normalize();
  if (d1) {
    if (d2)
      setBasePositions(d1->q_read(), d2->q_read());
    else {
      setBasePositions(d1->q_read());
    }
  }
}

void siconos::joints::PivotJointR::setBasePositions(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2) {
  KneeJointR::setBasePositions(q1, q2);

  // Assume that axes_[0] has been properly set (by constructor or setAxis)

  // The provided axis is the basis for the orthogonal plane
  // constraint relative to q1.
  // If provided in absolute coordinates, must be rotated to q1 frame.
  if (absoluteRef_) {
    boost::math::quaternion<double> rot1{siconos::geometry::rotquat(q1)}, quatBuff;

    // Move to q1 frame by unapplying q1 frame rotation
    quatBuff = (1.0 / rot1) * siconos::geometry::posquat(axes_[0]) * rot1;
    axes_[0] << quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4();
  }
  boost::math::quaternion<double> cq2q10;
  if (q2)
    // Initial orientation offset between q2 and q1.
    cq2q10 = 1.0 / siconos::geometry::rotquat(*q2) * siconos::geometry::rotquat(q1);
  else
    cq2q10 = 1. / boost::math::quaternion<double>(1, 0, 0, 0) * siconos::geometry::rotquat(q1);

  cq2q_ << cq2q10.R_component_1(), cq2q10.R_component_2(), cq2q10.R_component_3(),
      cq2q10.R_component_4();

  // Initialize the two orthogonal vectors that define the constraint plane.
  buildOrthonormalBase();

  // Get initial offsets relative to plane constraints.
  siconos::algebra::SiconosVector rot221{4};
  if (q2)
    pivot::rot2to1(q1.tail(4), q2->tail(4), cq2q_, rot221);
  else
    pivot::rot2to1(q1.tail(4), Eigen::Vector4d(1, 0, 0, 0), cq2q_, rot221);

  _initial_AscalA1 = axis1_.dot(rot221.head(3));
  _initial_AscalA2 = axis2_.dot(rot221.head(3));

  // In case of joint constraints, it's okay to use dot product=0, but
  // in the case of the free axis we must measure the actual angle
  // using atan2 so that stops can be placed correctly.
  double Adot2to1 = axes_[0].dot(rot221.head(3));
  _initial_AscalA = 2 * atan2(rot221(3), Adot2to1);

  _twistCount = 0;
  _previousAngle = 0;
}

void siconos::joints::PivotJointR::buildOrthonormalBase() {
  // Update axes_[0], A1 and A2 and check

  // Compute orthonormal basis (axes_[0],A1,A2) from axes_[0]
  auto base_ok = siconos::geometry::orthoBaseFromVector(axes_[0], axis1_, axis2_);
  if (!base_ok)
    THROW_EXCEPTION(
        "siconos::joints::PivotJointR::buildOrthonormalBase. Can not compute orthonormal "
        "vectors.")
}

void siconos::joints::PivotJointR::computeH_NE_(double time,
                                                siconos::modeling::Interaction& inter,
                                                const siconos::algebra::BlockVector& q0) {
  H_NE_view_->setZero();
  auto q1 = q0.vector(0);
  // Only the quaternion part of q is required to compute H (last 4 components)
  if (q0.numberOfBlocks() > 1) {
    auto q2 = q0.vector(1);
    pivot::computeH_for_2DS(q1->tail(4), G1P0_, q2->tail(4), G2P0_, axis1_, axis2_, cq2q_,
                            *H_NE_view_);
  } else
    pivot::computeH_for_1DS(q1->tail(4), G1P0_, axis1_, axis2_, cq2q_, *H_NE_view_);
}

void siconos::joints::PivotJointR::computeh(const siconos::algebra::BlockVector& q0,
                                            Eigen::Ref<siconos::algebra::SiconosVector> y) {
  KneeJointR::computeh(q0, y);

  auto q1 = q0.vector(0);
  siconos::algebra::SiconosVector3 rot221;
  if (q0.numberOfBlocks() == 2) {
    auto q2 = q0.vector(1);
    pivot::rot2to1(q1->tail(4), q2->tail(4), cq2q_, rot221);
  } else
    pivot::rot2to1(q1->tail(4), Eigen::Vector4d(1, 0, 0, 0), cq2q_, rot221);

  y(3) = axis1_.dot(rot221) - _initial_AscalA1;
  y(4) = axis2_.dot(rot221) - _initial_AscalA2;
}

/** Compute the vector of linear and angular positions of the free axes */
void siconos::joints::PivotJointR::computehDoF(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y, unsigned int axis) {
  // Normally we fill y starting at axis up to the number of columns,
  // but in this case there is only one, so just don't do anything if
  // it doesn't match.
  if (axis != 0) return;

  siconos::algebra::SiconosVector rot221{4};
  if (q2) {
    pivot::rot2to1(q1.tail(4), q2->tail(4), cq2q_, rot221);
  } else
    pivot::rot2to1(q1.tail(4), Eigen::Vector4d(1, 0, 0, 0), cq2q_, rot221);

  // In case of joint constraints, it's okay to use dot product=0, but
  // in the case of the free axis we must measure the actual angle
  // using atan2 so that stops can be placed correctly.
  double Adot2to1 = axes_[0].dot(rot221.head(3));
  double wrappedAngle = piwrap(2 * atan2(rot221(3), Adot2to1) - _initial_AscalA);

  // Count the number of twists around the angle, and report the
  // unwrapped angle.  Needed to implement joint stops near pi.
  if (wrappedAngle < -std::numbers::pi * 3 / 4 && _previousAngle > std::numbers::pi * 3 / 4)
    _twistCount++;
  else if (wrappedAngle > std::numbers::pi * 3 / 4 &&
           _previousAngle < -std::numbers::pi * 3 / 4)
    _twistCount--;
  _previousAngle = wrappedAngle;
  double unwrappedAngle = wrappedAngle + 2 * std::numbers::pi * _twistCount;

  y(0) = unwrappedAngle;
}

/** Compute the jacobian of linear and angular DoF with respect to some q */
void siconos::joints::PivotJointR::computeJachqDoF(
    siconos::modeling::Interaction& inter, const siconos::algebra::BlockVector& q0,
    Eigen::Ref<siconos::algebra::SiconosMatrix> jachq, unsigned int axis) {
  // Normally we fill jachq starting at axis up to the number of rows,
  // but in this case there is only one, so just don't do anything if
  // it doesn't match.
  if (axis != 0) return;

  auto q1 = (q0.getAllVect())[0];
  double q10 = (*q1)(3);
  double q11 = (*q1)(4);
  double q12 = (*q1)(5);
  double q13 = (*q1)(6);

  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;

  if (q0.numberOfBlocks() > 1) {
    auto q2 = (q0.getAllVect())[1];
    q20 = (*q2)(3);
    q21 = (*q2)(4);
    q22 = (*q2)(5);
    q23 = (*q2)(6);
  }

  jachq.setZero();

  /*
   * sympy expression:
   *
   * rot2to1 = qmul(qinv(qmul(q2,cq2q10)),q1)
   * Adot2to1 = np.dot(A, rot2to1)
   * angle = atan2(rot2to1[0], Adot2to1)
   *
   * [angle.diff(x) for x in q1]
   * r, e = cse(exprs=[angle.diff(x) for x in q1]
   *                 +[angle.diff(x) for x in q2],
   *            order='canonical')
   * for v,x in r: print('double {} = {};'.format(v,ccode(x)))
   * for i in range(4): print('jachq.setValue(0, {}, {});'.format(i+3,e[i]))
   */

  double x0 = cq2q_(2) * q23;
  double x1 = cq2q_(0) * q21;
  double x2 = cq2q_(1) * q20;
  double x3 = cq2q_(3) * q22;
  double x4 = x0 - x1 - x2 - x3;
  double x5 = cq2q_(3) * q21;
  double x6 = cq2q_(0) * q22;
  double x7 = cq2q_(1) * q23;
  double x8 = cq2q_(2) * q20;
  double x9 = x5 - x6 - x7 - x8;
  double x10 = cq2q_(1) * q22;
  double x11 = cq2q_(0) * q23;
  double x12 = cq2q_(2) * q21;
  double x13 = cq2q_(3) * q20;
  double x14 = x10 - x11 - x12 - x13;
  double x15 = q11 * x4;
  double x16 = q12 * x9;
  double x17 = q13 * x14;
  double x18 = cq2q_(0) * q20 - cq2q_(1) * q21 - cq2q_(2) * q22 - cq2q_(3) * q23;
  double x19 = q10 * x18;
  double x20 = axes_[0](0) * (q10 * x4 + q11 * x18 - q12 * x14 + q13 * x9) +
               axes_[0](1) * (q10 * x9 + q11 * x14 + q12 * x18 - q13 * x4) +
               axes_[0](2) * (q10 * x14 - q11 * x9 + q12 * x4 + q13 * x18);
  double x21 = 1.0 / (pow(x20, 2) + pow(-x15 - x16 - x17 + x19, 2));
  double x22 = 2 * x21 * (x15 + x16 + x17 - x19);
  double x23 = 2 * x20 * x21;
  double x24 = -x5 + x6 + x7 + x8;
  double x25 = -x0 + x1 + x2 + x3;
  double x26 = -x10 + x11 + x12 + x13;
  double x27 = cq2q_(0) * q11;
  double x28 = cq2q_(1) * q10;
  double x29 = -x28;
  double x30 = x27 + x29;
  double x31 = cq2q_(3) * q12;
  double x32 = cq2q_(2) * q13;
  double x33 = -x32;
  double x34 = cq2q_(0) * q12;
  double x35 = cq2q_(1) * q13;
  double x36 = cq2q_(2) * q10;
  double x37 = -x36;
  double x38 = cq2q_(3) * q11;
  double x39 = -x38;
  double x40 = x37 + x39;
  double x41 = cq2q_(0) * q13;
  double x42 = cq2q_(1) * q12;
  double x43 = -x42;
  double x44 = x41 + x43;
  double x45 = cq2q_(2) * q11;
  double x46 = cq2q_(3) * q10;
  double x47 = -x46;
  double x48 = cq2q_(0) * q10;
  double x49 = cq2q_(1) * q11;
  double x50 = cq2q_(2) * q12;
  double x51 = cq2q_(3) * q13;
  double x52 = x50 + x51;
  double x53 = -x48;
  double x54 = -x45;
  double x55 = -x35;
  double x56 = -x31;
  double x57 = x49 + x53;
  double x58 = x33 + x56;
  double x59 = x47 + x54;
  double x60 = x34 + x55;

  jachq.setValue(0, 3,
                 x18 * x23 + x22 * (axes_[0](0) * x4 + axes_[0](1) * x9 + axes_[0](2) * x14));
  jachq.setValue(
      0, 4, x22 * (axes_[0](0) * x18 + axes_[0](1) * x14 + axes_[0](2) * x24) + x23 * x25);
  jachq.setValue(0, 5,
                 x22 * (axes_[0](0) * x26 + axes_[0](1) * x18 + axes_[0](2) * x4) + x23 * x24);
  jachq.setValue(0, 6,
                 x22 * (axes_[0](0) * x9 + axes_[0](1) * x25 + axes_[0](2) * x18) + x23 * x26);

  if (q0.numberOfBlocks() < 2) return;

  /*
   * sympy expression:
   *
   * for i in range(4): print('jachq.setValue(0, {}, {});'.format(i+10,e[i+4]))
   */

  jachq.setValue(0, 10,
                 x22 * (axes_[0](0) * (x30 + x31 + x33) + axes_[0](1) * (x34 + x35 + x40) +
                        axes_[0](2) * (x44 + x45 + x47)) +
                     x23 * (x48 + x49 + x52));
  jachq.setValue(0, 11,
                 x22 * (axes_[0](0) * (-x49 + x52 + x53) + axes_[0](1) * (x44 + x46 + x54) +
                        axes_[0](2) * (-x34 + x40 + x55)) +
                     x23 * (x30 + x32 + x56));
  jachq.setValue(0, 12,
                 x22 * (axes_[0](0) * (-x41 + x43 + x59) + axes_[0](1) * (-x50 + x51 + x57) +
                        axes_[0](2) * (x27 + x28 + x58)) +
                     x23 * (x37 + x38 + x60));
  jachq.setValue(0, 13,
                 x22 * (axes_[0](0) * (x36 + x39 + x60) + axes_[0](1) * (-x27 + x29 + x58) +
                        axes_[0](2) * (x50 - x51 + x57)) +
                     x23 * (x41 + x42 + x59));
}

siconos::algebra::SiconosVector3 siconos::joints::PivotJointR::normalDoF(
    const siconos::algebra::SiconosVector& q0,
    const std::optional<Eigen::Ref<siconos::algebra::SiconosVector>>& q1, int axis,
    bool absoluteRef) {
  assert(axis == 0);
  if (axis != 0) return siconos::algebra::SiconosVector3{};

  // We assume that A is normalized.
  auto result = axes_[0];

  if (absoluteRef) siconos::geometry::rewriteVectorFromBodyToAbsoluteFrame(q0, result);
  return result;  // RVO
}

// Free functions

void siconos::joints::pivot::computeH_for_2DS(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp1,
    const siconos::algebra::SiconosVector3& coords1,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp2,
    const siconos::algebra::SiconosVector3& coords2,
    const siconos::algebra::SiconosVector3& A1, const siconos::algebra::SiconosVector3& A2,
    const siconos::algebra::SiconosVector& cq2q,
    Eigen::Ref<siconos::algebra::MapType> result) {
  knee::computeH_for_2DS(qp1, coords1, qp2, coords2, result.topRows(3));

  // sympy expression: [AscalA1.diff(x) for x in q1]
  result.setValue(
      3, 3,
      A1(0) * (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)) +
          A1(1) *
              (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)) +
          A1(2) *
              (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)));
  result.setValue(
      3, 4,
      A1(0) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)) +
          A1(1) *
              (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)) +
          A1(2) * (cq2q(0) * qp2(2) + cq2q(1) * qp2(3) + cq2q(2) * qp2(0) - cq2q(3) * qp2(1)));
  result.setValue(
      3, 5,
      A1(0) * (cq2q(0) * qp2(3) - cq2q(1) * qp2(2) + cq2q(2) * qp2(1) + cq2q(3) * qp2(0)) +
          A1(1) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)) +
          A1(2) *
              (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)));
  result.setValue(
      3, 6,
      A1(0) * (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)) +
          A1(1) * (cq2q(0) * qp2(1) + cq2q(1) * qp2(0) - cq2q(2) * qp2(3) + cq2q(3) * qp2(2)) +
          A1(2) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)));

  result.setValue(3, 7, 0);
  result.setValue(3, 8, 0);
  result.setValue(3, 9, 0);

  // sympy expression: [AscalA1.diff(x) for x in q2]
  result.setValue(
      3, 10,
      A1(0) * (cq2q(0) * qp1(1) - cq2q(1) * qp1(0) - cq2q(2) * qp1(3) + cq2q(3) * qp1(2)) +
          A1(1) * (cq2q(0) * qp1(2) + cq2q(1) * qp1(3) - cq2q(2) * qp1(0) - cq2q(3) * qp1(1)) +
          A1(2) * (cq2q(0) * qp1(3) - cq2q(1) * qp1(2) + cq2q(2) * qp1(1) - cq2q(3) * qp1(0)));
  result.setValue(
      3, 11,
      A1(0) * (-cq2q(0) * qp1(0) - cq2q(1) * qp1(1) + cq2q(2) * qp1(2) + cq2q(3) * qp1(3)) +
          A1(1) * (cq2q(0) * qp1(3) - cq2q(1) * qp1(2) - cq2q(2) * qp1(1) + cq2q(3) * qp1(0)) +
          A1(2) *
              (-cq2q(0) * qp1(2) - cq2q(1) * qp1(3) - cq2q(2) * qp1(0) - cq2q(3) * qp1(1)));
  result.setValue(
      3, 12,
      A1(0) * (-cq2q(0) * qp1(3) - cq2q(1) * qp1(2) - cq2q(2) * qp1(1) - cq2q(3) * qp1(0)) +
          A1(1) *
              (-cq2q(0) * qp1(0) + cq2q(1) * qp1(1) - cq2q(2) * qp1(2) + cq2q(3) * qp1(3)) +
          A1(2) * (cq2q(0) * qp1(1) + cq2q(1) * qp1(0) - cq2q(2) * qp1(3) - cq2q(3) * qp1(2)));
  result.setValue(
      3, 13,
      A1(0) * (cq2q(0) * qp1(2) - cq2q(1) * qp1(3) + cq2q(2) * qp1(0) - cq2q(3) * qp1(1)) +
          A1(1) *
              (-cq2q(0) * qp1(1) - cq2q(1) * qp1(0) - cq2q(2) * qp1(3) - cq2q(3) * qp1(2)) +
          A1(2) *
              (-cq2q(0) * qp1(0) + cq2q(1) * qp1(1) + cq2q(2) * qp1(2) - cq2q(3) * qp1(3)));

  // sympy expression: [AscalA2.diff(x) for x in q1]
  result.setValue(
      4, 3,
      A2(0) * (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)) +
          A2(1) *
              (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)) +
          A2(2) *
              (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)));
  result.setValue(
      4, 4,
      A2(0) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)) +
          A2(1) *
              (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)) +
          A2(2) * (cq2q(0) * qp2(2) + cq2q(1) * qp2(3) + cq2q(2) * qp2(0) - cq2q(3) * qp2(1)));
  result.setValue(
      4, 5,
      A2(0) * (cq2q(0) * qp2(3) - cq2q(1) * qp2(2) + cq2q(2) * qp2(1) + cq2q(3) * qp2(0)) +
          A2(1) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)) +
          A2(2) *
              (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)));
  result.setValue(
      4, 6,
      A2(0) * (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)) +
          A2(1) * (cq2q(0) * qp2(1) + cq2q(1) * qp2(0) - cq2q(2) * qp2(3) + cq2q(3) * qp2(2)) +
          A2(2) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)));

  // sympy expression: [AscalA2.diff(x) for x in q1]
  result.setValue(
      4, 10,
      A2(0) * (cq2q(0) * qp1(1) - cq2q(1) * qp1(0) - cq2q(2) * qp1(3) + cq2q(3) * qp1(2)) +
          A2(1) * (cq2q(0) * qp1(2) + cq2q(1) * qp1(3) - cq2q(2) * qp1(0) - cq2q(3) * qp1(1)) +
          A2(2) * (cq2q(0) * qp1(3) - cq2q(1) * qp1(2) + cq2q(2) * qp1(1) - cq2q(3) * qp1(0)));
  result.setValue(
      4, 11,
      A2(0) * (-cq2q(0) * qp1(0) - cq2q(1) * qp1(1) + cq2q(2) * qp1(2) + cq2q(3) * qp1(3)) +
          A2(1) * (cq2q(0) * qp1(3) - cq2q(1) * qp1(2) - cq2q(2) * qp1(1) + cq2q(3) * qp1(0)) +
          A2(2) *
              (-cq2q(0) * qp1(2) - cq2q(1) * qp1(3) - cq2q(2) * qp1(0) - cq2q(3) * qp1(1)));
  result.setValue(
      4, 12,
      A2(0) * (-cq2q(0) * qp1(3) - cq2q(1) * qp1(2) - cq2q(2) * qp1(1) - cq2q(3) * qp1(0)) +
          A2(1) *
              (-cq2q(0) * qp1(0) + cq2q(1) * qp1(1) - cq2q(2) * qp1(2) + cq2q(3) * qp1(3)) +
          A2(2) * (cq2q(0) * qp1(1) + cq2q(1) * qp1(0) - cq2q(2) * qp1(3) - cq2q(3) * qp1(2)));
  result.setValue(
      4, 13,
      A2(0) * (cq2q(0) * qp1(2) - cq2q(1) * qp1(3) + cq2q(2) * qp1(0) - cq2q(3) * qp1(1)) +
          A2(1) *
              (-cq2q(0) * qp1(1) - cq2q(1) * qp1(0) - cq2q(2) * qp1(3) - cq2q(3) * qp1(2)) +
          A2(2) *
              (-cq2q(0) * qp1(0) + cq2q(1) * qp1(1) + cq2q(2) * qp1(2) - cq2q(3) * qp1(3)));

  /*proj_with_q
  for (siconos::algebra::Index ii=0; ii <result.rows(); ii++)
    for (siconos::algebra::Index jj=0; jj <result.cols(); jj++)
  H_NE_view_Proj->setValue(ii,jj,result(ii, jj));

  H_NE_view_Proj->setValue(5,0,0);
  H_NE_view_Proj->setValue(5,1,0);
  H_NE_view_Proj->setValue(5,2,0);
  H_NE_view_Proj->setValue(5,3,2.0*q10);
  H_NE_view_Proj->setValue(5,4,2.0*q11);
  H_NE_view_Proj->setValue(5,5,2.0*q12);
  H_NE_view_Proj->setValue(5,6,2.0*q13);
  H_NE_view_Proj->setValue(6,0+7,0);
  H_NE_view_Proj->setValue(6,1+7,0);
  H_NE_view_Proj->setValue(6,2+7,0);
  H_NE_view_Proj->setValue(6,3+7,2.0*q20);
  H_NE_view_Proj->setValue(6,4+7,2.0*q21);
  H_NE_view_Proj->setValue(6,5+7,2.0*q22);
  H_NE_view_Proj->setValue(6,6+7,2.0*q23);
  */

  // siconos::algebra::print(result);
}

void siconos::joints::pivot::computeH_for_1DS(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp1,
    const siconos::algebra::SiconosVector3& coords1,
    const siconos::algebra::SiconosVector3& A1, const siconos::algebra::SiconosVector3& A2,
    const siconos::algebra::SiconosVector& cq2q,
    Eigen::Ref<siconos::algebra::MapType> result) {
  knee::computeH_for_1DS(qp1, coords1, result.topRows(3));
  // sympy expression: [AscalA1.diff(x) for x in q1]
  result.setValue(3, 3, A1(0) * (-cq2q(1)) + A1(1) * (-cq2q(2)) + A1(2) * (-cq2q(3)));
  result.setValue(3, 4, A1(0) * (cq2q(0)) + A1(1) * (-cq2q(3)) + A1(2) * (cq2q(2)));
  result.setValue(3, 5, A1(0) * (cq2q(3)) + A1(1) * (cq2q(0)) + A1(2) * (-cq2q(1)));
  result.setValue(3, 6, A1(0) * (-cq2q(2)) + A1(1) * (cq2q(1)) + A1(2) * (cq2q(0)));
  // sympy expression: [AscalA2.diff(x) for x in q1]
  result.setValue(4, 3, A2(0) * (-cq2q(1)) + A2(1) * (-cq2q(2)) + A2(2) * (-cq2q(3)));
  result.setValue(4, 4, A2(0) * (cq2q(0)) + A2(1) * (-cq2q(3)) + A2(2) * (cq2q(2)));
  result.setValue(4, 5, A2(0) * (cq2q(3)) + A2(1) * (cq2q(0)) + A2(2) * (-cq2q(1)));
  result.setValue(4, 6, A2(0) * (-cq2q(2)) + A2(1) * (cq2q(1)) + A2(2) * (cq2q(0)));

  /*proj_with_q
      for (siconos::algebra::Index ii=0; ii <result.rows(); ii++)
        for (siconos::algebra::Index jj=0; jj <result.cols(); jj++)
    H_NE_view_Proj->setValue(ii,jj,result(ii, jj));

      H_NE_view_Proj->setValue(5,0,0);
      H_NE_view_Proj->setValue(5,1,0);
      H_NE_view_Proj->setValue(5,2,0);
      H_NE_view_Proj->setValue(5,3,2.0*q10);
      H_NE_view_Proj->setValue(5,4,2.0*q11);
      H_NE_view_Proj->setValue(5,5,2.0*q12);
      H_NE_view_Proj->setValue(5,6,2.0*q13);
  */
}

void siconos::joints::pivot::rot2to1(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp1,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& qp2,
    const siconos::algebra::SiconosVector& cq2q,
    Eigen::Ref<siconos::algebra::SiconosVector> result) {
  /*
   * The current rotation vector taking into account initial rotation
   * difference.
   *
   * sympy expression:
   * rot2to1 = qmul(qinv(qmul(q2,cq2q10)),q1)
   */

  result(0) =
      (qp1(0) * (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)) +
       qp1(1) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)) -
       qp1(2) * (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)) +
       qp1(3) * (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)));
  result(1) =
      (qp1(0) * (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)) +
       qp1(1) * (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)) +
       qp1(2) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)) -
       qp1(3) * (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)));
  result(2) =
      (qp1(0) * (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)) -
       qp1(1) * (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)) +
       qp1(2) * (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)) +
       qp1(3) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)));
  if (result.size() == 4)
    result(3) =
        (qp1(0) * (cq2q(0) * qp2(0) - cq2q(1) * qp2(1) - cq2q(2) * qp2(2) - cq2q(3) * qp2(3)) -
         qp1(1) *
             (-cq2q(0) * qp2(1) - cq2q(1) * qp2(0) + cq2q(2) * qp2(3) - cq2q(3) * qp2(2)) -
         qp1(2) *
             (-cq2q(0) * qp2(2) - cq2q(1) * qp2(3) - cq2q(2) * qp2(0) + cq2q(3) * qp2(1)) -
         qp1(3) *
             (-cq2q(0) * qp2(3) + cq2q(1) * qp2(2) - cq2q(2) * qp2(1) - cq2q(3) * qp2(0)));
}
