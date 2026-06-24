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
// #include <numbers>  // pi

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "siconos_debug.h"

namespace siconos::geometry {

void computeRotationMatrix(double q0, double q1, double q2, double q3,
                           siconos::algebra::SiconosMatrix33 &rotationMatrix) {
  /* Direct computation using the quaternion-to-rotation-matrix formula
   * See: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
   *
   * For a unit quaternion q = [w, x, y, z], the rotation matrix R is:
   * R = | 1-2y²-2z²   2xy-2zw     2xz+2yw   |
   *     | 2xy+2zw     1-2x²-2z²   2yz-2xw   |
   *     | 2xz-2yw     2yz+2xw     1-2x²-2y² |
   *
   * The implementation below uses an equivalent expanded form.
   */
  rotationMatrix(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
  rotationMatrix(0, 1) = 2.0 * (q1 * q2 - q0 * q3);
  rotationMatrix(0, 2) = 2.0 * (q1 * q3 + q0 * q2);
  rotationMatrix(1, 0) = 2.0 * (q1 * q2 + q0 * q3);
  rotationMatrix(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
  rotationMatrix(1, 2) = 2.0 * (q2 * q3 - q0 * q1);
  rotationMatrix(2, 0) = 2.0 * (q1 * q3 - q0 * q2);
  rotationMatrix(2, 1) = 2.0 * (q2 * q3 + q0 * q1);
  rotationMatrix(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
}

void rotateVector(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                  Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  DEBUG_BEGIN("::rotateVector(q,v)\n");
  DEBUG_EXPR(siconos::algebra::print(v););
  DEBUG_EXPR(std::cout << std::scientific << std::setprecision(12) << std::setw(16)
                       << "q[3:6] " << q.tail(4) << "\n";);

  // Efficient rotation using quaternion cross product formula (Rodrigues' rotation formula)
  // For unit quaternion q = [w, v], rotated vector v' = v + 2*v × (v × v_original) + 2*w*(v × v_original)
  // This is more efficient than building the full rotation matrix.
  siconos::algebra::ConstMapVector3Type qvect(q.data() + 4);  // view onto q4, q5, q6 (vector part)
  auto q0 = q(3);  // scalar part
  siconos::algebra::SiconosVector3 t = 2 * qvect.cross(v);
  v += qvect.cross(t);
  v += q0 * t;
  DEBUG_EXPR(std::cout << v << "\n";);
  DEBUG_END("::rotateVector(q,v)\n");
}

void rotateMatrix(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                  Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  DEBUG_BEGIN("::rotateMatrix(q,m)\n");
  DEBUG_EXPR(std::cout << m << "\n";);
  DEBUG_EXPR(std::cout << std::scientific << std::setprecision(12) << std::setw(16)
                       << "q[3:6] " << q.tail(4) << "\n";);

  // Apply quaternion rotation to each column of the matrix
  siconos::algebra::ConstMapVector3Type qvect(q.data() + 4);  // view onto q4, q5, q6
  auto q0 = q(3);
  for (siconos::algebra::Index j = 0; j < m.cols(); j++) {
    Eigen::Map<Eigen::Vector3d> mcol(m.col(j).data());
    Eigen::Vector3d t = 2 * qvect.cross(mcol);
    mcol += qvect.cross(t) + q0 * t;
  }
  DEBUG_EXPR(std::cout << m << "\n";);
  DEBUG_END("::rotateMatrix(q,m)\n");
}

void rotateVectorFromInertialToBodyFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  DEBUG_BEGIN("::rotateVectorToBodyFrame(q,v)\n");
  // To rotate from inertial to body frame, we use the conjugate quaternion
  // q_conj = [w, -x, -y, -z]
  siconos::algebra::SiconosVector7 qbis;
  qbis << q(0), q(1), q(2), q(3), -q(4), -q(5), -q(6);
  siconos::geometry::rotateVector(qbis, v);
  DEBUG_END("::rotateVectorToBodyFrame(q,v)\n");
}

void rotateMatrixFromInertialToBodyFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  DEBUG_BEGIN("::rotateMatrixToBodyFrame(q,m)\n");
  // Use conjugate quaternion for inverse rotation
  siconos::algebra::SiconosVector7 qbis;
  qbis << q(0), q(1), q(2), q(3), -q(4), -q(5), -q(6);
  siconos::geometry::rotateMatrix(qbis, m);
  DEBUG_END("::rotateMatrixToBodyFrame(q,m)\n");
}

void rotateVectorFromBodyToInertialFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  DEBUG_BEGIN("::rotateVectorToInertialFrame(q,v)\n");
  siconos::geometry::rotateVector(q, v);
  DEBUG_END("::rotateVectorToInertialFrame(q,v)\n");
}

void rotateMatrixFromBodyToInertialFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  DEBUG_BEGIN("::rotateMatrixToInertialFrame(q,m )\n");
  siconos::geometry::rotateMatrix(q, m);
  DEBUG_END("::rotateMatrixToInertialFrame(q,m )\n");
}

void computeRotationMatrix(const siconos::algebra::SiconosVector7 &q,
                           siconos::algebra::SiconosMatrix33 &rotationMatrix) {
  siconos::geometry::computeRotationMatrix(q(3), q(4), q(5), q(6), rotationMatrix);
}

void computeRotationMatrixTransposed(const siconos::algebra::SiconosVector7 &q,
                                     siconos::algebra::SiconosMatrix33 &rotationMatrix) {
  // Transposed rotation matrix = inverse rotation = rotation by conjugate quaternion
  siconos::geometry::computeRotationMatrix(q(3), -q(4), -q(5), -q(6), rotationMatrix);
}

double axisAngleFromQuaternion(double q0, double q1, double q2, double q3,
                               Eigen::Ref<siconos::algebra::SiconosVector3> &axis) {
  DEBUG_BEGIN(
      "axisAngleFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> axis )\n");
  // angle = 2 * acos(w), axis = [x, y, z] / sin(angle/2)
  double angle = acos(q0) * 2.0;
  // Using sqrt(1-w²) is cheaper and more numerically stable than sin(angle/2)
  double f = sqrt(1 - q0 * q0);
  if (f > std::numeric_limits<double>::epsilon()) {
    axis(0) = q1 / f;
    axis(1) = q2 / f;
    axis(2) = q3 / f;
  } else {
    // For zero rotation, axis is undefined - set to zero
    axis.setZero();
  }
  DEBUG_PRINTF("angle= %12.8e\n", angle);
  DEBUG_EXPR(siconos::algebra::print(*axis);)
  DEBUG_END(
      "axisAngleFromQuaternion(double q0, double q1, double q2, double q3, "
      "std::shared_ptr<siconos::algebra::SiconosVector> axis )\n");
  return angle;
}

double axisAngleFromConfiguration(const Eigen::Ref<siconos::algebra::SiconosVector7> &q,
                                  Eigen::Ref<siconos::algebra::SiconosVector3> axis) {
  double angle = siconos::geometry::axisAngleFromQuaternion(q(3), q(4), q(5), q(6), axis);
  return angle;
}

siconos::algebra::SiconosVector3 rotationVectorFromQuaternion(double q0, double q1,
                                                               double q2, double q3) {
  DEBUG_BEGIN("rotationVectorFromQuaternion(...)\n");
  siconos::algebra::SiconosVector3 rotationVector;

  // Rotation vector r = 2 * acos(w) * [x, y, z] / |[x, y, z]|
  // For small rotations, this is approximately 2 * [x, y, z]
  rotationVector << q1, q2, q3;
  double norm_v = rotationVector.norm();
  norm_v = std::clamp(norm_v, 0., 1.);  // Clamp to valid range for asin
  if (norm_v < 1e-12) {
    rotationVector.setZero();
  } else {
    rotationVector *= 2.0 * std::asin(norm_v) / norm_v;
  }
  DEBUG_EXPR(siconos::algebra::print(*rotationVector);)
  DEBUG_END("rotationVectorFromQuaternion(...)\n");
  return rotationVector;
}

siconos::algebra::SiconosVector3 rotationVectorFromConfiguration(
    siconos::algebra::SiconosVector7 &q) {
  return siconos::geometry::rotationVectorFromQuaternion(q(3), q(4), q(5), q(6));
}

void quaternionFromAxisAngle(const siconos::algebra::SiconosVector3 &axis, double angle,
                             siconos::algebra::SiconosVector7 &q) {
  // q = [cos(θ/2), axis * sin(θ/2)]
  q(3) = cos(angle / 2.0);
  q(4) = axis(0) * sin(angle * 0.5);
  q(5) = axis(1) * sin(angle * 0.5);
  q(6) = axis(2) * sin(angle * 0.5);
}

double sinc(double x) {
  if (std::abs(x) <= 1e-3) {
    // Taylor expansion for small x to avoid 0/0
    // sin(x)/x ≈ 1 + x²/6 + x⁴/120 + ...
    // Here we use coefficients for a different expansion order
    return 1.0 + x * x / 3.0 + pow(x, 4) * 2.0 / 15.0 + pow(x, 6) * 17.0 / 315.0 +
           pow(x, 8) * 62.0 / 2835.0;
  } else {
    return sin(x) / x;
  }
}

siconos::algebra::SiconosVector7 quaternionFromRotationVector(
    const siconos::algebra::SiconosVector3 &rotationVector) {
  siconos::algebra::SiconosVector7 q;
  q.setZero();
  double angle = rotationVector.norm();
  double half_angle = angle * 0.5;
  double f = 0.5 * sinc(half_angle);

  q(3) = cos(half_angle);
  q(4) = rotationVector(0) * f;
  q(5) = rotationVector(1) * f;
  q(6) = rotationVector(2) * f;
  return q;
}

double quaternionNorm(const siconos::algebra::SiconosVector7 &q) {
  double normq = sqrt(q(3) * q(3) + q(4) * q(4) + q(5) * q(5) + q(6) * q(6));
  return normq;
}

void normalizeQuaternion(Eigen::Ref<siconos::algebra::SiconosVector7> q) {
  auto quat = q.tail<4>();
  double norm = quat.norm();
  if (norm > std::numeric_limits<double>::epsilon())
    quat /= norm;
  else
    THROW_EXCEPTION("normalizeQuaternion: quaternion part has zero norm");
}

void quaternionFromTwistVector(const siconos::algebra::SiconosVector6 &twist,
                               Eigen::Ref<siconos::algebra::SiconosVector7> q) {
  // Twist vector [vx, vy, vz, wx, wy, wz] -> quaternion from angular part
  double angle = sqrt(twist(3) * twist(3) + twist(4) * twist(4) + twist(5) * twist(5));
  double f = 0.5 * sinc(angle * 0.5);

  q(3) = cos(angle / 2.0);
  q(4) = twist(3) * f;
  q(5) = twist(4) * f;
  q(6) = twist(5) * f;
}

void compositionLawLieGroup(const siconos::algebra::SiconosVector7 &a,
                            siconos::algebra::SiconosVector7 &b,
                            siconos::algebra::SiconosVector7 &ab) {
  // SE(3) composition: ab = a ∘ b
  // Position: ab_pos = a_pos + R(a_rot) * b_pos
  ab(0) = a(0) + b(0);
  ab(1) = a(1) + b(1);
  ab(2) = a(2) + b(2);

  // Rotation: ab_rot = a_rot * b_rot (quaternion product)
  boost::math::quaternion<double> quat_a(a(3), a(4), a(5), a(6));
  boost::math::quaternion<double> quat_b(b(3), b(4), b(5), b(6));
  boost::math::quaternion<double> quat_ab = quat_a * quat_b;
  ab(3) = quat_ab.R_component_1();
  ab(4) = quat_ab.R_component_2();
  ab(5) = quat_ab.R_component_3();
  ab(6) = quat_ab.R_component_4();
}

void compositionLawLieGroup(const siconos::algebra::SiconosVector7 &a,
                            Eigen::Ref<siconos::algebra::SiconosVector7> b) {
  // In-place composition: b = a ∘ b
  b(0) = a(0) + b(0);
  b(1) = a(1) + b(1);
  b(2) = a(2) + b(2);

  boost::math::quaternion<double> quat_a(a(3), a(4), a(5), a(6));
  boost::math::quaternion<double> quat_b(b(3), b(4), b(5), b(6));
  boost::math::quaternion<double> quat_ab = quat_a * quat_b;
  b(3) = quat_ab.R_component_1();
  b(4) = quat_ab.R_component_2();
  b(5) = quat_ab.R_component_3();
  b(6) = quat_ab.R_component_4();
}

void extractRotationQuaternion(const siconos::algebra::SiconosVector7 &from,
                               boost::math::quaternion<double> &to) {
  to = boost::math::quaternion<double>{from(3), from(4), from(5), from(6)};
}

void extractVectorFromQuaternion(const boost::math::quaternion<double> &from,
                                 siconos::algebra::SiconosVector3 &to) {
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
  to(2) = from.R_component_4();
}

void extractPositionToQuaternion(const siconos::algebra::SiconosVector7 &from,
                                 boost::math::quaternion<double> &to) {
  to = boost::math::quaternion<double>{0, from(0), from(1), from(2)};
}

void extractRotationQuaternion2d(const siconos::algebra::SiconosVector &from,
                                 boost::math::quaternion<double> &to) {
  // 2D rotation around z-axis: quaternion [cos(θ/2), 0, 0, sin(θ/2)]
  double half_angle = from(2) / 2.0;
  to = boost::math::quaternion<double>{cos(half_angle), 0.0, 0.0, sin(half_angle)};
}

void extractPositionFromQuaternion2d(const boost::math::quaternion<double> &from,
                                     siconos::algebra::SiconosVector &to) {
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
}

void extractPositionToQuaternion2d(const siconos::algebra::SiconosVector &from,
                                   boost::math::quaternion<double> &to) {
  to = boost::math::quaternion<double>{0, from(0), from(1), 0.0};
}

boost::math::quaternion<double> getRotationQuaternion(const siconos::algebra::SiconosVector7 &v) {
  if (!v.isZero())
    return boost::math::quaternion<double>(v(3), v(4), v(5), v(6));
  else
    return boost::math::quaternion<double>(1, 0, 0, 0);
}

boost::math::quaternion<double> getPositionQuaternion(const siconos::algebra::SiconosVector &v) {
  return boost::math::quaternion<double>{0, v(0), v(1), v(2)};
}

void computeOrthonormalBaseFromAxis(siconos::algebra::SiconosVector3 &axis0,
                                    siconos::algebra::SiconosVector3 &axis1,
                                    siconos::algebra::SiconosVector3 &axis2) {
  if (axis0.norm() < 1e-10)
    throw std::invalid_argument(
        "input vector has a norm equal to zero, can't compute a base.");

  axis0.normalize();

  // Choose an arbitrary vector not parallel to axis0
  siconos::algebra::SiconosVector3 arbitrary(1.0, 0.0, 0.0);
  if (std::abs(axis0.dot(arbitrary)) > 0.99)
    arbitrary = siconos::algebra::SiconosVector3(0.0, 1.0, 0.0);

  // Gram-Schmidt orthogonalization
  axis1 = axis0.cross(arbitrary).normalized();
  axis2 = axis0.cross(axis1);
}

bool orthoBaseFromVector(siconos::algebra::SiconosVector3 &A,
                         siconos::algebra::SiconosVector3 &A1,
                         siconos::algebra::SiconosVector3 &A2) {
  double normA = A.norm();
  if (normA == 0.0) {
    // Return NaN for invalid input
    A = siconos::algebra::SiconosVector3::Constant(std::numeric_limits<double>::quiet_NaN());
    A1 = siconos::algebra::SiconosVector3::Constant(std::numeric_limits<double>::quiet_NaN());
    A2 = siconos::algebra::SiconosVector3::Constant(std::numeric_limits<double>::quiet_NaN());
    return false;
  }

  // Normalize A
  A.normalize();

  // Branchless orthonormal basis construction (Duff et al. 2017)
  double sign = std::copysign(1.0, A.z());
  const double a = -1.0 / (sign + A.z());
  const double b = A.x() * A.y() * a;

  // Build orthonormal basis
  A1 << 1.0 + sign * A.x() * A.x() * a, sign * b, -sign * A.x();
  A2 << b, sign + A.y() * A.y() * a, -A.y();

  // Verify orthonormality (debug builds only)
  assert(std::fabs(A1.norm() - 1.0) < 1e-14);
  assert(std::fabs(A.dot(A1)) < 1e-14);
  assert(std::fabs(A2.norm() - 1.0) < 1e-14);
  assert(std::fabs(A.dot(A2)) < 1e-14);
  assert(std::fabs(A1.dot(A2)) < 1e-14);

  return true;
}

}  // namespace siconos::geometry
