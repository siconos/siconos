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

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "siconos_debug.h"

void siconos::geometry::computeRotationMatrix(
    double q0, double q1, double q2, double q3,
    siconos::algebra::SiconosMatrix33 &rotationMatrix) {
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
  rotationMatrix.setValue(0, 0, q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3);
  rotationMatrix.setValue(0, 1, 2.0 * (q1 * q2 - q0 * q3));
  rotationMatrix.setValue(0, 2, 2.0 * (q1 * q3 + q0 * q2));
  rotationMatrix.setValue(1, 0, 2.0 * (q1 * q2 + q0 * q3));
  rotationMatrix.setValue(1, 1, q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3);
  rotationMatrix.setValue(1, 2, 2.0 * (q2 * q3 - q0 * q1));
  rotationMatrix.setValue(2, 0, 2.0 * (q1 * q3 - q0 * q2));
  rotationMatrix.setValue(2, 1, 2.0 * (q2 * q3 + q0 * q1));
  rotationMatrix.setValue(2, 2, q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3);
}

void siconos::geometry::quaternionRotateVector(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  DEBUG_BEGIN("::quaternionRotateVector(q,v)\n");
  DEBUG_EXPR(siconos::algebra::print(v););
  DEBUG_EXPR(std::cout << std::scientific << std::setprecision(12) << std::setw(16)
                       << "q[3:6] " << q.tail(4) << "\n";);

  // First way. Using the rotation matrix
  // std::shared_ptr<siconos::algebra::SiconosMatrix> rotationMatrix(new
  // siconos::algebra::SiconosMatrix(3,3)); siconos::algebra::SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,q1,q2,q3, rotationMatrix);
  // tmp.noalias() = *rotationMatrix * v;
  // v = tmp;
  // return;

  // Second way. Using the transpose of the rotation matrix
  // std::shared_ptr<siconos::algebra::SiconosMatrix> rotationMatrix(new
  // siconos::algebra::SiconosMatrix(3,3)); siconos::algebra::SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,-q1,-q2,-q3, rotationMatrix);
  // tmp.noalias() = *rotationMatrix * v;
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
  siconos::algebra::ConstMapVector3Type qvect(q.data() + 4);  // view onto  q4, 5, 6
  auto q0 = q(3);
  siconos::algebra::SiconosVector3 t = 2 * qvect.cross(v);
  v += qvect.cross(t);
  v += q0 * t;
  DEBUG_EXPR(std::cout << v << "\n";);
  DEBUG_END("::quaternionRotateVector(q,v)\n");
}

void siconos::geometry::quaternionRotateMatrix(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  DEBUG_BEGIN("::quaternionRotateMatrix(q,m)\n");
  DEBUG_EXPR(std::cout << m << "\n";);
  DEBUG_EXPR(std::cout << std::scientific << std::setprecision(12) << std::setw(16)
                       << "q[3:6] " << q.tail(4) << "\n";);

  // Direct computation with cross product for each column
  siconos::algebra::ConstMapVector3Type qvect(q.data() + 4);  // view onto  q4, 5, 6
  auto q0 = q(3);
  for (unsigned int j = 0; j < m.cols(); j++) {
    Eigen::Map<Eigen::Vector3d> mcol(m.col(j).data());
    auto t = 2 * qvect.cross(mcol);
    mcol += qvect.cross(t) + q0 * t;
  }
  DEBUG_EXPR(std::cout << m << "\n";);
  DEBUG_END("::quaternionRotateMatrix(q,m)\n");
}

void siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  DEBUG_BEGIN("::rewriteVectorFromAbsoluteToBodyFrame(q,v)\n");
  siconos::algebra::SiconosVector7 qbis;
  qbis << q(0), q(1), q(2), q(3), -q(4), -q(5), -q(6);
  siconos::geometry::quaternionRotateVector(qbis, v);
  DEBUG_END("::rewriteVectorFromAbsoluteToBodyFrame(q,v)\n");
}

void siconos::geometry::rewriteMatrixFromAbsoluteToBodyFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  DEBUG_BEGIN("::rewriteMatrixFromAbsoluteToBodyFrame(q,m)\n");
  siconos::algebra::SiconosVector7 qbis;
  qbis << q(0), q(1), q(2), q(3), -q(4), -q(5), -q(6);
  siconos::geometry::quaternionRotateMatrix(qbis, m);
  DEBUG_END("::rewriteMatrixFromAbsoluteToBodyFrame(q,m)\n");
}

void siconos::geometry::rewriteVectorFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  DEBUG_BEGIN("::rewriteVectorFromBodyToAbsoluteFrame(q,v)\n");
  siconos::geometry::quaternionRotateVector(q, v);
  DEBUG_END("::rewriteVectorFromBodyToAbsoluteFrame(q,v)\n");
}

void siconos::geometry::rewriteMatrixFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  DEBUG_BEGIN("::rewriteMatrixFromBodyToAbsoluteFrame(q,m )\n");
  siconos::geometry::quaternionRotateMatrix(q, m);
  DEBUG_END("::rewriteMatrixFromBodyToAbsoluteFrame(q,m )\n");
}

void siconos::geometry::computeRotationMatrix(
    const siconos::algebra::SiconosVector7 &q,
    siconos::algebra::SiconosMatrix33 &rotationMatrix) {
  siconos::geometry::computeRotationMatrix(q(3), q(4), q(5), q(6), rotationMatrix);
}
void siconos::geometry::computeRotationMatrixTransposed(
    const siconos::algebra::SiconosVector7 &q,
    siconos::algebra::SiconosMatrix33 &rotationMatrix) {
  siconos::geometry::computeRotationMatrix(q(3), -q(4), -q(5), -q(6), rotationMatrix);
}

double siconos::geometry::axisAngleFromQuaternion(
    double q0, double q1, double q2, double q3,
    Eigen::Ref<siconos::algebra::SiconosVector3> &axis) {
  DEBUG_BEGIN(
      "axisAngleFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> axis )\n");
  double angle = acos(q0) * 2.0;
  // double f = sin( angle *0.5);
  double f = sqrt(1 - q0 * q0);  // cheaper than sin ?
  if (f != 0.0) {
    axis(0) = q1 / f;
    axis(1) = q2 / f;
    axis(2) = q3 / f;
  } else {
    axis.setZero();
  }
  DEBUG_PRINTF("angle= %12.8e\n", angle);
  DEBUG_EXPR(siconos::algebra::print(*axis);)
  DEBUG_END(
      "axisAngleFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> axis )\n");
  return angle;
}

double siconos::geometry::axisAngleFromConfiguration(
    const Eigen::Ref<siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> axis) {
  double angle = siconos::geometry::axisAngleFromQuaternion(q(3), q(4), q(5), q(6), axis);
  return angle;
}

void siconos::geometry::rotationVectorFromQuaternion(
    double q0, double q1, double q2, double q3,
    siconos::algebra::SiconosVector &rotationVector) {
  DEBUG_BEGIN(
      "rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> rotationVector )\n");

  rotationVector(0) = q1;
  rotationVector(1) = q2;
  rotationVector(2) = q3;

  double norm_v = sqrt(q1 * q1 + q2 * q2 + q3 * q3);
  assert(norm_v <= M_PI); /* it should be called for a unit quaternion */
  if (norm_v < 1e-12) {
    rotationVector(0) = 0.0;
    rotationVector(1) = 0.0;
    rotationVector(2) = 0.0;
  } else {
    rotationVector *= 2.0 * asin(norm_v) / norm_v;
  }
  DEBUG_EXPR(siconos::algebra::print(*rotationVector);)
  DEBUG_END(
      "rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> rotationVector )\n");
}

void siconos::geometry::rotationVectorFromConfiguration(
    siconos::algebra::SiconosVector &q, siconos::algebra::SiconosVector &rotationVector) {
  siconos::geometry::rotationVectorFromQuaternion(q(3), q(4), q(5), q(6), rotationVector);
}

void siconos::geometry::quaternionFromAxisAngle(const siconos::algebra::SiconosVector3 &axis,
                                                double angle,
                                                siconos::algebra::SiconosVector7 &q) {
  q(3) = cos(angle / 2.0);
  q(4) = axis(0) * sin(angle * 0.5);
  q(5) = axis(1) * sin(angle * 0.5);
  q(6) = axis(2) * sin(angle * 0.5);
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
    const siconos::algebra::SiconosVector &rotationVector,
    siconos::algebra::SiconosVector &q) {
  double angle =
      sqrt(rotationVector(0) * rotationVector(0) + rotationVector(1) * rotationVector(1) +
           rotationVector(2) * rotationVector(2));

  double f = 0.5 * sin_x(angle * 0.5);

  q(3) = cos(angle / 2.0);
  q(4) = rotationVector(0) * f;
  q(5) = rotationVector(1) * f;
  q(6) = rotationVector(2) * f;
}

double siconos::geometry::quaternionNorm(const siconos::algebra::SiconosVector7 &q) {
  double normq = sqrt(q(3) * q(3) + q(4) * q(4) + q(5) * q(5) + q(6) * q(6));
  return normq;
}

void siconos::geometry::normalizeq(Eigen::Ref<siconos::algebra::SiconosVector7> q) {
  double normq = sqrt(q(3) * q(3) + q(4) * q(4) + q(5) * q(5) + q(6) * q(6));
  assert(normq > 0);
  normq = 1.0 / normq;
  q(3) = q(3) * normq;
  q(4) = q(4) * normq;
  q(5) = q(5) * normq;
  q(6) = q(6) * normq;
}

void siconos::geometry::quaternionFromTwistVector(
    const siconos::algebra::SiconosVector6 &twist,
    Eigen::Ref<siconos::algebra::SiconosVector7> q) {
  double angle = sqrt(twist(3) * twist(3) + twist(4) * twist(4) + twist(5) * twist(5));

  double f = 0.5 * sin_x(angle * 0.5);

  q(3) = cos(angle / 2.0);
  q(4) = twist(3) * f;
  q(5) = twist(4) * f;
  q(6) = twist(5) * f;
}
void siconos::geometry::compositionLawLieGroup(const siconos::algebra::SiconosVector7 &a,
                                               siconos::algebra::SiconosVector7 &b,
                                               siconos::algebra::SiconosVector7 &ab) {
  // For the translational component, the composition law is the addition
  ab(0) = a(0) + b(0);
  ab(1) = a(1) + b(1);
  ab(2) = a(2) + b(2);

  // For the quaternion that encodes rotation, the composition law is the quaternion product.
  boost::math::quaternion<double> quat_a(a(3), a(4), a(5), a(6));
  boost::math::quaternion<double> quat_b(b(3), b(4), b(5), b(6));
  boost::math::quaternion<double> quat_ab = quat_a * quat_b;
  ab(3) = quat_ab.R_component_1();
  ab(4) = quat_ab.R_component_2();
  ab(5) = quat_ab.R_component_3();
  ab(6) = quat_ab.R_component_4();
}

void siconos::geometry::compositionLawLieGroup(
    const siconos::algebra::SiconosVector7 &a,
    Eigen::Ref<siconos::algebra::SiconosVector7> b) {
  // For the translational component, the composition law is the addition
  b(0) = a(0) + b(0);
  b(1) = a(1) + b(1);
  b(2) = a(2) + b(2);

  // For the quaternion that encodes rotation, the composition law is the quaternion product.
  boost::math::quaternion<double> quat_a(a(3), a(4), a(5), a(6));
  boost::math::quaternion<double> quat_b(b(3), b(4), b(5), b(6));
  boost::math::quaternion<double> quat_ab = quat_a * quat_b;
  b(3) = quat_ab.R_component_1();
  b(4) = quat_ab.R_component_2();
  b(5) = quat_ab.R_component_3();
  b(6) = quat_ab.R_component_4();
  // normalizeq(b);
}

void siconos::geometry::copyQuatRot(const siconos::algebra::SiconosVector7 &from,
                                    boost::math::quaternion<double> &to) {
  to = boost::math::quaternion<double>{from(3), from(4), from(5), from(6)};
}

void siconos::geometry::copyQuatPos(const boost::math::quaternion<double> &from,
                                    siconos::algebra::SiconosVector &to) {
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
  to(2) = from.R_component_4();
}

void siconos::geometry::copyQuatPos(const siconos::algebra::SiconosVector &from,
                                    boost::math::quaternion<double> &to) {
  to = boost::math::quaternion<double>{0, from(0), from(1), from(2)};
}

void siconos::geometry::copyQuatRot2d(const siconos::algebra::SiconosVector &from,
                                      boost::math::quaternion<double> &to) {
  double half_angle = from(2) / 2.0;

  to = boost::math::quaternion<double>{cos(half_angle), 0.0, 0.0, sin(half_angle)};
}

void siconos::geometry::copyQuatPos2d(const boost::math::quaternion<double> &from,
                                      siconos::algebra::SiconosVector &to) {
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
}

void siconos::geometry::copyQuatPos2d(const siconos::algebra::SiconosVector &from,
                                      boost::math::quaternion<double> &to) {
  to = boost::math::quaternion<double>{0, from(0), from(1), 0.0};
}

boost::math::quaternion<double> siconos::geometry::rotquat(
    const siconos::algebra::SiconosVector7 &v) {
  if (!v.isZero())
    return boost::math::quaternion<double>(v(3), v(4), v(5), v(6));
  else
    return boost::math::quaternion<double>(1, 0, 0, 0);
}

boost::math::quaternion<double> siconos::geometry::posquat(
    const siconos::algebra::SiconosVector &v) {
  return boost::math::quaternion<double>{0, v(0), v(1), v(2)};
}

void siconos::geometry::computeOrthonormalBaseFromAxis(
    siconos::algebra::SiconosVector3 &axis0, siconos::algebra::SiconosVector3 &axis1,
    siconos::algebra::SiconosVector3 &axis2) {
  if (axis0.norm() < 1e-10)
    throw std::invalid_argument(
        "input vector has a norm equal to zero, can't compute a base.");

  axis0.normalize();

  siconos::algebra::SiconosVector3 arbitrary(1.0, 0.0, 0.0);
  if (std::abs(axis0.dot(arbitrary)) > 0.99)
    arbitrary = siconos::algebra::SiconosVector3(0.0, 1.0, 0.0);

  axis1 = axis0.cross(arbitrary).normalized();
  axis2 = axis0.cross(axis1);
}

bool siconos::geometry::orthoBaseFromVector(siconos::algebra::SiconosVector3 &A,
                                            siconos::algebra::SiconosVector3 &A1,
                                            siconos::algebra::SiconosVector3 &A2) {
  double normA = A.norm();
  if (normA == 0.0) {
    // If A is null, we assign Nan to outputs  and return an error code
    A1 = siconos::algebra::SiconosVector3(std::numeric_limits<double>::quiet_NaN(),
                                          std::numeric_limits<double>::quiet_NaN(),
                                          std::numeric_limits<double>::quiet_NaN());
    A2 = siconos::algebra::SiconosVector3(std::numeric_limits<double>::quiet_NaN(),
                                          std::numeric_limits<double>::quiet_NaN(),
                                          std::numeric_limits<double>::quiet_NaN());
    return false;
  }

  // Normalize A
  A.normalize();

  double sign = std::copysign(1.0, A.z());
  const double a = -1.0 / (sign + A.z());
  const double b = A.x() * A.y() * a;

  // Build the orthonormal basis using a and b
  A1 << 1.0 + sign * A.x() * A.x() * a, sign * b, -sign * A.x();
  A2 << b, sign + A.y() * A.y() * a, -A.y();

  // and normalize them
  // A1.normalize();
  // A2.normalize();

  // Check norms ...
  assert(std::fabs(A1.norm() - 1.0) < 1e-14);
  assert(std::fabs(A.dot(A1)) < 1e-14);
  assert(std::fabs(A2.norm() - 1.0) < 1e-14);
  assert(std::fabs(A.dot(A2)) < 1e-14);
  assert(std::fabs(A1.dot(A2)) < 1e-14);

  return true;
}
