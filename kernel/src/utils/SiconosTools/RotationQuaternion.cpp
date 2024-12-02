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
// #define DEBUG_NOCOLOR
// #define DEBUG_BEGIN_END_ONLY
//  #define DEBUG_STDOUT
//  #define DEBUG_MESSAGES
#include "RotationQuaternion.hpp"

#include <boost/math/quaternion.hpp>

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "siconos_debug.h"

void siconos::geometry::computeRotationMatrix(
    double q0, double q1, double q2, double q3,
    std::shared_ptr<siconos::algebra::SimpleMatrix> rotationMatrix) {
  /* Brute force version by multiplication of quaternion
   */
  // boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  // boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  // boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  // boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  // boost::math::quaternion<double>    quatz(0, 0, 0, 1);
  // boost::math::quaternion<double>    quatBuff;
  // quatBuff = quatQ * quatx * quatcQ;
  // rotationMatrix->setValue(0, 0, quatBuff.R_component_2());
  // rotationMatrix->setValue(1, 0, quatBuff.R_component_3());
  // rotationMatrix->setValue(2, 0, quatBuff.R_component_4());
  // quatBuff = quatQ * quaty * quatcQ;
  // rotationMatrix->setValue(0, 1, quatBuff.R_component_2());
  // rotationMatrix->setValue(1, 1, quatBuff.R_component_3());
  // rotationMatrix->setValue(2, 1, quatBuff.R_component_4());
  // quatBuff = quatQ * quatz * quatcQ;
  // rotationMatrix->setValue(0, 2, quatBuff.R_component_2());
  // rotationMatrix->setValue(1, 2, quatBuff.R_component_3());
  // rotationMatrix->setValue(2, 2, quatBuff.R_component_4());

  /* direct computation https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation */
  rotationMatrix->setValue(0, 0, q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3);
  rotationMatrix->setValue(0, 1, 2.0 * (q1 * q2 - q0 * q3));
  rotationMatrix->setValue(0, 2, 2.0 * (q1 * q3 + q0 * q2));

  rotationMatrix->setValue(1, 0, 2.0 * (q1 * q2 + q0 * q3));
  rotationMatrix->setValue(1, 1, q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3);
  rotationMatrix->setValue(1, 2, 2.0 * (q2 * q3 - q0 * q1));

  rotationMatrix->setValue(2, 0, 2.0 * (q1 * q3 - q0 * q2));
  rotationMatrix->setValue(2, 1, 2.0 * (q2 * q3 + q0 * q1));
  rotationMatrix->setValue(2, 2, q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3);
}

void siconos::geometry::quaternionRotate(double q0, double q1, double q2, double q3,
                                         siconos::algebra::SiconosVector& v) {
  DEBUG_BEGIN(
      "::quaternionRotate(double q0, double q1, double q2, double q3, "
      "siconos::algebra::SiconosVector& v )\n");
  DEBUG_EXPR(v.display(););
  DEBUG_PRINTF("( q0 = %16.12e,  q1 = %16.12e,  q2= %16.12e,  q3= %16.12e )\n", q0, q1, q2,
               q3);
  assert(v.size() == 3);

  // First way. Using the rotation matrix
  // std::shared_ptr<siconos::algebra::SimpleMatrix> rotationMatrix(new
  // siconos::algebra::SimpleMatrix(3,3)); siconos::algebra::SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,q1,q2,q3, rotationMatrix);
  // prod(*rotationMatrix, v, tmp);
  // v = tmp;
  // return;

  // Second way. Using the transpose of the rotation matrix
  // std::shared_ptr<siconos::algebra::SimpleMatrix> rotationMatrix(new
  // siconos::algebra::SimpleMatrix(3,3)); siconos::algebra::SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,-q1,-q2,-q3, rotationMatrix);
  // prod(v, *rotationMatrix, tmp);
  // v = tmp;

  // Third way. cross product and axis angle
  // see http://www.geometrictools.com/Documentation/RotationIssues.pdf
  // std::shared_ptr<siconos::algebra::SiconosVector> axis(new
  // siconos::algebra::SiconosVector(3)); double angle = ::axisAngleFromQuaternion(q0,q1,q2,q3,
  // axis); siconos::algebra::SiconosVector t(3), tmp(3); cross_product(*axis,v,t);
  // cross_product(*axis,t,tmp);
  // v += sin(angle)*t + (1.0-cos(angle))*tmp;

  // Direct computation with cross product
  // Works only with unit quaternion
  siconos::algebra::SiconosVector t(3), tmp(3);
  siconos::algebra::SiconosVector qvect(3);
  qvect(0) = q1;
  qvect(1) = q2;
  qvect(2) = q3;
  cross_product(qvect, v, t);
  t *= 2.0;
  cross_product(qvect, t, tmp);
  v += tmp;
  v += q0 * t;
  DEBUG_EXPR(v.display(););
  DEBUG_END(
      "::quaternionRotate(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v )\n");
}

void siconos::geometry::quaternionRotate(double q0, double q1, double q2, double q3,
                                         std::shared_ptr<siconos::algebra::SiconosVector> v) {
  siconos::geometry::quaternionRotate(q0, q1, q2, q3, *v);
}

void siconos::geometry::quaternionRotate(double q0, double q1, double q2, double q3,
                                         std::shared_ptr<siconos::algebra::SimpleMatrix> m) {
  DEBUG_BEGIN(
      "::quaternionRotate(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m )\n");
  DEBUG_EXPR(m->display(););
  DEBUG_PRINTF("( q0 = %16.12e,  q1 = %16.12e,  q2= %16.12e,  q3= %16.12e )\n", q0, q1, q2,
               q3);

  // Direct computation with cross product for each column
  assert(m->size(0) == 3 &&
         "::quaternionRotate(double q0, double q1, double q2, double q3, "
         "std::shared_ptr<siconos::algebra::SimpleMatrix> m ) m must have 3 rows");
  siconos::algebra::SiconosVector v(3);
  siconos::algebra::SiconosVector t(3), tmp(3);
  siconos::algebra::SiconosVector qvect(3);
  qvect(0) = q1;
  qvect(1) = q2;
  qvect(2) = q3;
  for (unsigned int j = 0; j < m->size(1); j++) {
    v(0) = m->getValue(0, j);
    v(1) = m->getValue(1, j);
    v(2) = m->getValue(2, j);
    cross_product(qvect, v, t);
    t *= 2.0;
    cross_product(qvect, t, tmp);
    v += tmp;
    v += q0 * t;
    m->setValue(0, j, v(0));
    m->setValue(1, j, v(1));
    m->setValue(2, j, v(2));
  }
  DEBUG_EXPR(m->display(););
  DEBUG_END(
      "::quaternionRotate(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m )\n");
}

void siconos::geometry::quaternionRotate(std::shared_ptr<siconos::algebra::SiconosVector> q,
                                         std::shared_ptr<siconos::algebra::SiconosVector> v) {
  DEBUG_BEGIN(
      "::quaternionRotate(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v )\n");
  siconos::geometry::quaternionRotate(q->getValue(3), q->getValue(4), q->getValue(5),
                                      q->getValue(6), v);
  DEBUG_END(
      "::quaternionRotate(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v )\n");
}

void siconos::geometry::quaternionRotate(std::shared_ptr<siconos::algebra::SiconosVector> q,
                                         std::shared_ptr<siconos::algebra::SimpleMatrix> m) {
  DEBUG_BEGIN(
      "::quaternionRotate(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m )\n");
  siconos::geometry::quaternionRotate(q->getValue(3), q->getValue(4), q->getValue(5),
                                      q->getValue(6), m);
  DEBUG_END(
      "::quaternionRotate(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m)\n");
}

void siconos::geometry::changeFrameAbsToBody(const siconos::algebra::SiconosVector& q,
                                             siconos::algebra::SiconosVector& v) {
  DEBUG_BEGIN(
      "::changeFrameAbsToBody(const siconos::algebra::SiconosVector& q, "
      "siconos::algebra::SiconosVector& v )\n");
  siconos::geometry::quaternionRotate(q.getValue(3), -q.getValue(4), -q.getValue(5),
                                      -q.getValue(6), v);
  DEBUG_END(
      "::changeFrameAbsToBody(const siconos::algebra::SiconosVector& q, "
      "siconos::algebra::SiconosVector& v )\n");
}
void siconos::geometry::changeFrameAbsToBody(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  DEBUG_BEGIN(
      "::changeFrameAbsToBody(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v )\n");
  siconos::geometry::quaternionRotate(q->getValue(3), -q->getValue(4), -q->getValue(5),
                                      -q->getValue(6), v);
  DEBUG_END(
      "::changeFrameAbsToBody(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v )\n");
}
void siconos::geometry::changeFrameAbsToBody(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SimpleMatrix> m) {
  DEBUG_BEGIN(
      "::changeFrameAbsToBody(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m )\n");
  siconos::geometry::quaternionRotate(q->getValue(3), -q->getValue(4), -q->getValue(5),
                                      -q->getValue(6), m);
  DEBUG_END(
      "::changeFrameAbsToBody(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m )\n");
}

void siconos::geometry::changeFrameBodyToAbs(const siconos::algebra::SiconosVector& q,
                                             siconos::algebra::SiconosVector& v) {
  DEBUG_BEGIN(
      "::changeFrameBodyToAbs(const siconos::algebra::SiconosVector& q, "
      "siconos::algebra::SiconosVector& v )\n");
  siconos::geometry::quaternionRotate(q.getValue(3), q.getValue(4), q.getValue(5),
                                      q.getValue(6), v);
  DEBUG_END(
      "::changeFrameBodyToAbs(const siconos::algebra::SiconosVector& q, "
      "siconos::algebra::SiconosVector& v )\n");
}
void siconos::geometry::changeFrameBodyToAbs(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  DEBUG_BEGIN(
      "::changeFrameBodyToAbs(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v )\n");
  siconos::geometry::quaternionRotate(q->getValue(3), q->getValue(4), q->getValue(5),
                                      q->getValue(6), *v);
  DEBUG_END(
      "::changeFrameBodyToAbs(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v )\n");
}
void siconos::geometry::changeFrameBodyToAbs(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SimpleMatrix> m) {
  DEBUG_BEGIN(
      "::changeFrameBodyToAbs(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m )\n");
  siconos::geometry::quaternionRotate(q->getValue(3), q->getValue(4), q->getValue(5),
                                      q->getValue(6), m);
  DEBUG_END(
      "::changeFrameBodyToAbs(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SimpleMatrix> m )\n");
}

void siconos::geometry::computeRotationMatrix(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SimpleMatrix> rotationMatrix) {
  siconos::geometry::computeRotationMatrix(q->getValue(3), q->getValue(4), q->getValue(5),
                                           q->getValue(6), rotationMatrix);
}
void siconos::geometry::computeRotationMatrixTransposed(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SimpleMatrix> rotationMatrix) {
  siconos::geometry::computeRotationMatrix(q->getValue(3), -q->getValue(4), -q->getValue(5),
                                           -q->getValue(6), rotationMatrix);
}

double siconos::geometry::axisAngleFromQuaternion(
    double q0, double q1, double q2, double q3,
    std::shared_ptr<siconos::algebra::SiconosVector> axis) {
  DEBUG_BEGIN(
      "axisAngleFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> axis )\n");
  double angle = acos(q0) * 2.0;
  // double f = sin( angle *0.5);
  double f = sqrt(1 - q0 * q0);  // cheaper than sin ?
  if (f != 0.0) {
    axis->setValue(0, q1 / f);
    axis->setValue(1, q2 / f);
    axis->setValue(2, q3 / f);
  } else {
    axis->zero();
  }
  DEBUG_PRINTF("angle= %12.8e\n", angle);
  DEBUG_EXPR(axis->display(););
  DEBUG_END(
      "axisAngleFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> axis )\n");
  return angle;
}

double siconos::geometry::axisAngleFromConfiguration(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> axis) {
  double angle = siconos::geometry::axisAngleFromQuaternion(
      q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6), axis);
  return angle;
}

void siconos::geometry::rotationVectorFromQuaternion(
    double q0, double q1, double q2, double q3,
    std::shared_ptr<siconos::algebra::SiconosVector> rotationVector) {
  DEBUG_BEGIN(
      "rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> rotationVector )\n");

  rotationVector->setValue(0, q1);
  rotationVector->setValue(1, q2);
  rotationVector->setValue(2, q3);

  double norm_v = sqrt(q1 * q1 + q2 * q2 + q3 * q3);
  assert(norm_v <= M_PI); /* it should be called for a unit quaternion */
  if (norm_v < 1e-12) {
    rotationVector->setValue(0, 0.0);
    rotationVector->setValue(1, 0.0);
    rotationVector->setValue(2, 0.0);
  } else {
    *rotationVector *= 2.0 * asin(norm_v) / norm_v;
  }
  DEBUG_EXPR(rotationVector->display(););
  DEBUG_END(
      "rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> rotationVector )\n");
}

void siconos::geometry::rotationVectorFromConfiguration(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> rotationVector) {
  siconos::geometry::rotationVectorFromQuaternion(
      q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6), rotationVector);
}

void siconos::geometry::quaternionFromAxisAngle(
    std::shared_ptr<siconos::algebra::SiconosVector> axis, double angle,
    std::shared_ptr<siconos::algebra::SiconosVector> q) {
  q->setValue(3, cos(angle / 2.0));
  q->setValue(4, axis->getValue(0) * sin(angle * 0.5));
  q->setValue(5, axis->getValue(1) * sin(angle * 0.5));
  q->setValue(6, axis->getValue(2) * sin(angle * 0.5));
}

double siconos::geometry::sin_x(double x) {
  if (std::abs(x) <= 1e-3) {
    return 1.0 + x * x / 3.0 + pow(x, 4) * 2.0 / 15.0 + pow(x, 6) * 17.0 / 315.0 +
           pow(x, 8) * 62.0 / 2835.0;
  } else {
    return sin(x) / x;
  }
}

void siconos::geometry::quaternionFromRotationVector(
    std::shared_ptr<siconos::algebra::SiconosVector> rotationVector,
    std::shared_ptr<siconos::algebra::SiconosVector> q) {
  double angle = sqrt(rotationVector->getValue(0) * rotationVector->getValue(0) +
                      rotationVector->getValue(1) * rotationVector->getValue(1) +
                      rotationVector->getValue(2) * rotationVector->getValue(2));

  double f = 0.5 * sin_x(angle * 0.5);

  q->setValue(3, cos(angle / 2.0));
  q->setValue(4, rotationVector->getValue(0) * f);
  q->setValue(5, rotationVector->getValue(1) * f);
  q->setValue(6, rotationVector->getValue(2) * f);
}

double siconos::geometry::quaternionNorm(const siconos::algebra::SiconosVector& q) {
  double normq = sqrt(q.getValue(3) * q.getValue(3) + q.getValue(4) * q.getValue(4) +
                      q.getValue(5) * q.getValue(5) + q.getValue(6) * q.getValue(6));
  return normq;
}

void siconos::geometry::normalizeq(std::shared_ptr<siconos::algebra::SiconosVector> q) {
  double normq = sqrt(q->getValue(3) * q->getValue(3) + q->getValue(4) * q->getValue(4) +
                      q->getValue(5) * q->getValue(5) + q->getValue(6) * q->getValue(6));
  assert(normq > 0);
  normq = 1.0 / normq;
  q->setValue(3, q->getValue(3) * normq);
  q->setValue(4, q->getValue(4) * normq);
  q->setValue(5, q->getValue(5) * normq);
  q->setValue(6, q->getValue(6) * normq);
}

void siconos::geometry::normalizeq(siconos::algebra::SiconosVector& q) {
  double normq = sqrt(q.getValue(3) * q.getValue(3) + q.getValue(4) * q.getValue(4) +
                      q.getValue(5) * q.getValue(5) + q.getValue(6) * q.getValue(6));
  assert(normq > 0);
  normq = 1.0 / normq;
  q.setValue(3, q.getValue(3) * normq);
  q.setValue(4, q.getValue(4) * normq);
  q.setValue(5, q.getValue(5) * normq);
  q.setValue(6, q.getValue(6) * normq);
}

void siconos::geometry::quaternionFromTwistVector(siconos::algebra::SiconosVector& twist,
                                                  siconos::algebra::SiconosVector& q) {
  assert(twist.size() == 6);
  assert(q.size() == 7);
  double angle =
      sqrt(twist.getValue(3) * twist.getValue(3) + twist.getValue(4) * twist.getValue(4) +
           twist.getValue(5) * twist.getValue(5));

  double f = 0.5 * sin_x(angle * 0.5);

  q.setValue(3, cos(angle / 2.0));
  q.setValue(4, twist.getValue(3) * f);
  q.setValue(5, twist.getValue(4) * f);
  q.setValue(6, twist.getValue(5) * f);
}
void siconos::geometry::compositionLawLieGroup(const siconos::algebra::SiconosVector& a,
                                               siconos::algebra::SiconosVector& b,
                                               siconos::algebra::SiconosVector& ab) {
  assert(a.size() == 7);
  assert(b.size() == 7);
  assert(ab.size() == 7);

  // For the translational component, the composition law is the addition
  ab.setValue(0, a.getValue(0) + b.getValue(0));
  ab.setValue(1, a.getValue(1) + b.getValue(1));
  ab.setValue(2, a.getValue(2) + b.getValue(2));

  // For the quaternion that encodes rotation, the composition law is the quaternion product.
  boost::math::quaternion<double> quat_a(a.getValue(3), a.getValue(4), a.getValue(5),
                                         a.getValue(6));
  boost::math::quaternion<double> quat_b(b.getValue(3), b.getValue(4), b.getValue(5),
                                         b.getValue(6));
  boost::math::quaternion<double> quat_ab = quat_a * quat_b;
  ab.setValue(3, quat_ab.R_component_1());
  ab.setValue(4, quat_ab.R_component_2());
  ab.setValue(5, quat_ab.R_component_3());
  ab.setValue(6, quat_ab.R_component_4());
}

void siconos::geometry::compositionLawLieGroup(const siconos::algebra::SiconosVector& a,
                                               siconos::algebra::SiconosVector& b) {
  assert(a.size() == 7);
  assert(b.size() == 7);

  // For the translational component, the composition law is the addition
  b.setValue(0, a.getValue(0) + b.getValue(0));
  b.setValue(1, a.getValue(1) + b.getValue(1));
  b.setValue(2, a.getValue(2) + b.getValue(2));

  // For the quaternion that encodes rotation, the composition law is the quaternion product.
  boost::math::quaternion<double> quat_a(a.getValue(3), a.getValue(4), a.getValue(5),
                                         a.getValue(6));
  boost::math::quaternion<double> quat_b(b.getValue(3), b.getValue(4), b.getValue(5),
                                         b.getValue(6));
  boost::math::quaternion<double> quat_ab = quat_a * quat_b;
  b.setValue(3, quat_ab.R_component_1());
  b.setValue(4, quat_ab.R_component_2());
  b.setValue(5, quat_ab.R_component_3());
  b.setValue(6, quat_ab.R_component_4());
  // normalizeq(b);
}

void siconos::geometry::copyQuatRot(const siconos::algebra::SiconosVector& from,
                                    boost::math::quaternion<double>& to) {
  to = boost::math::quaternion<double>{from(3), from(4), from(5), from(6)};
}

void siconos::geometry::copyQuatPos(const boost::math::quaternion<double>& from,
                                    siconos::algebra::SiconosVector& to) {
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
  to(2) = from.R_component_4();
}

void siconos::geometry::copyQuatPos(const siconos::algebra::SiconosVector& from,
                                    boost::math::quaternion<double>& to) {
  to = boost::math::quaternion<double>{0, from(0), from(1), from(2)};
}

void siconos::geometry::copyQuatRot2d(const siconos::algebra::SiconosVector& from,
                                      boost::math::quaternion<double>& to) {
  double half_angle = from(2) / 2.0;

  to = boost::math::quaternion<double>{cos(half_angle), 0.0, 0.0, sin(half_angle)};
}

void siconos::geometry::copyQuatPos2d(const boost::math::quaternion<double>& from,
                                      siconos::algebra::SiconosVector& to) {
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
}

void siconos::geometry::copyQuatPos2d(const siconos::algebra::SiconosVector& from,
                                      boost::math::quaternion<double>& to) {
  to = boost::math::quaternion<double>{0, from(0), from(1), 0.0};
}

boost::math::quaternion<double> siconos::geometry::rotquat(
    const std::shared_ptr<siconos::algebra::SiconosVector>& v) {
  if (v)
    return boost::math::quaternion<double>(v->getValue(3), v->getValue(4), v->getValue(5),
                                           v->getValue(6));
  else
    return boost::math::quaternion<double>(1, 0, 0, 0);
}

boost::math::quaternion<double> siconos::geometry::posquat(
    const std::shared_ptr<siconos::algebra::SiconosVector>& v) {
  return boost::math::quaternion<double>{0, v->getValue(0), v->getValue(1), v->getValue(2)};
}
