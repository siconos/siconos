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

/*! \file RotationQuaternion.hpp
 *  \brief Quaternion-based rotation utilities for rigid body dynamics
 *
 *  This file provides a comprehensive set of utilities for working with rotations
 *  in 3D space using quaternions. The functions support:
 *  - Conversions between quaternions, rotation matrices, axis-angle, and rotation vectors
 *  - Rotating vectors and matrices by quaternions
 *  - Transforming between body (local) and inertial (absolute) reference frames
 *  - Composition of rigid body configurations (SE(3) Lie group operations)
 *
 *  \section sec_conventions Conventions
 *
 *  - Quaternions are stored as 4D vectors [w, x, y, z] where w is the scalar part
 *  - Configurations (poses) are stored as 7D vectors [x, y, z, qw, qx, qy, qz]
 *    where (x,y,z) is the position and (qw,qx,qy,qz) is the unit quaternion
 *  - Active rotations: the vector is rotated, not the coordinate frame
 *  - Right-handed coordinate systems
 *
 *  \section sec_usage Usage Example
 *  \code
 *  // Create a configuration (position + orientation)
 *  siconos::algebra::SiconosVector7 q;
 *  q << 1.0, 2.0, 3.0,  // position
 *       0.707, 0.0, 0.707, 0.0;  // quaternion (90 deg around y-axis)
 *
 *  // Rotate a vector from body frame to inertial frame
 *  siconos::algebra::SiconosVector3 v_body(1.0, 0.0, 0.0);
 *  siconos::geometry::rotateVectorToInertialFrame(q, v_body);
 *
 *  // Get rotation matrix
 *  siconos::algebra::SiconosMatrix33 R;
 *  siconos::geometry::computeRotationMatrix(q, R);
 *  \endcode
 */

#ifndef ROTATIONQUATERNION_H
#define ROTATIONQUATERNION_H

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// Include full quaternion header (provides boost::math::quaternion)
#include <boost/math/quaternion.hpp>

namespace siconos::algebra {

// class SiconosVector;
// class SiconosMatrix;
}  // namespace siconos::algebra

namespace siconos::geometry {

/** @defgroup QuaternionConversions Quaternion Conversions
 *  Functions for converting between different rotation representations
 */

/** @defgroup QuaternionOperations Quaternion Operations
 *  Functions for applying rotations and composing configurations
 */

/** @defgroup FrameTransformations Frame Transformations
 *  Functions for transforming vectors/matrices between reference frames
 */

/** @defgroup QuaternionAccessors Quaternion Accessors
 *  Functions for extracting components from configuration vectors
 */

/** \brief Convert a quaternion to axis-angle representation
 *
 *  Given a unit quaternion q = [w, x, y, z], compute the equivalent
 *  rotation as an axis-angle pair (axis, angle).
 *
 *  The conversion follows:
 *  - angle = 2 * acos(w)
 *  - axis = [x, y, z] / sin(angle/2)  (if angle != 0)
 *
 *  \param[in] q0 scalar part (w) of the quaternion
 *  \param[in] q1 first vector component (x) of the quaternion
 *  \param[in] q2 second vector component (y) of the quaternion
 *  \param[in] q3 third vector component (z) of the quaternion
 *  \param[out] axis the rotation axis (unit vector), modified in-place
 *  \return the rotation angle in radians (range [0, 2π])
 *  \note If the angle is zero, axis is set to zero vector
 *  \ingroup QuaternionConversions
 */
double axisAngleFromQuaternion(double q0, double q1, double q2, double q3,
                               Eigen::Ref<siconos::algebra::SiconosVector3> &axis);

/** \brief Extract axis-angle from a configuration vector
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[out] axis the rotation axis (unit vector), modified in-place
 *  \return the rotation angle in radians
 *  \see axisAngleFromQuaternion
 *  \ingroup QuaternionConversions
 */
double axisAngleFromConfiguration(const Eigen::Ref<siconos::algebra::SiconosVector7> &q,
                                  Eigen::Ref<siconos::algebra::SiconosVector3> axis);

/** \brief Convert a quaternion to rotation vector (exponential map)
 *
 *  The rotation vector r is related to the quaternion by:
 *  - r = 2 * acos(w) * [x, y, z] / |[x, y, z]|
 *
 *  This is the inverse of quaternionFromRotationVector.
 *
 *  \param[in] q0 scalar part (w) of the quaternion
 *  \param[in] q1 first vector component (x) of the quaternion
 *  \param[in] q2 second vector component (y) of the quaternion
 *  \param[in] q3 third vector component (z) of the quaternion
 *  \return the rotation vector (axis * angle)
 *  \ingroup QuaternionConversions
 */
siconos::algebra::SiconosVector3 rotationVectorFromQuaternion(double q0, double q1, double q2,
                                                              double q3);

/** \brief Extract rotation vector from a configuration
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \return the rotation vector
 *  \see rotationVectorFromQuaternion
 *  \ingroup QuaternionConversions
 */
siconos::algebra::SiconosVector3 rotationVectorFromConfiguration(
    siconos::algebra::SiconosVector7 &q);

/** \brief Convert axis-angle to quaternion and store in configuration
 *
 *  \param[in] axis the rotation axis (will be normalized)
 *  \param[in] angle the rotation angle in radians
 *  \param[out] q configuration vector [x, y, z, qw, qx, qy, qz], modified in-place.
 *               The position part (first 3 components) is left unchanged.
 *  \ingroup QuaternionConversions
 */
void quaternionFromAxisAngle(const siconos::algebra::SiconosVector3 &axis, double angle,
                             siconos::algebra::SiconosVector7 &q);

/** \brief Convert rotation vector to quaternion
 *
 *  The rotation vector r = angle * axis is converted to quaternion using:
 *  - w = cos(|r|/2)
 *  - [x, y, z] = sin(|r|/2) * r / |r|
 *
 *  This is the inverse of rotationVectorFromQuaternion.
 *
 *  \param[in] rotationVector the rotation vector (axis * angle)
 *  \return configuration vector with quaternion set, position zeroed
 *  \ingroup QuaternionConversions
 */
siconos::algebra::SiconosVector7 quaternionFromRotationVector(
    const siconos::algebra::SiconosVector3 &rotationVector);

/** \brief Compute sin(x)/x with Taylor expansion for small x
 *
 *  Uses a 4th-order Taylor expansion for |x| <= 1e-3 to avoid numerical issues:
 *  sinc(x) ≈ 1 + x²/3 + 2x⁴/15 + 17x⁶/315 + 62x⁸/2835
 *
 *  \param[in] x the input value
 *  \return sin(x)/x
 *  \note For x=0, returns 1.0 (the limit value)
 *  \ingroup QuaternionConversions
 */
double sinc(double x);

/** \brief Convert twist vector (spatial velocity) to quaternion increment
 *
 *  Used in Lie group integration to convert angular velocity integrated
 *  over a timestep into a quaternion rotation increment.
 *
 *  \param[in] twist 6D twist vector [vx, vy, vz, wx, wy, wz] (linear + angular)
 *  \param[out] q configuration vector with quaternion set from angular part
 *  \ingroup QuaternionConversions
 */
void quaternionFromTwistVector(const siconos::algebra::SiconosVector6 &twist,
                               Eigen::Ref<siconos::algebra::SiconosVector7> q);

/** \brief Compute the norm of the quaternion part of a configuration
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \return the norm of the quaternion part sqrt(qw² + qx² + qy² + qz²)
 *  \ingroup QuaternionOperations
 */
double quaternionNorm(const siconos::algebra::SiconosVector7 &q);

/** \brief Normalize the quaternion part of a configuration
 *
 *  Ensures the quaternion part has unit norm, which is required for
 *  valid rotations. Throws if the norm is too small.
 *
 *  \param[in,out] q configuration vector, quaternion part normalized in-place
 *  \throw SiconosException if quaternion norm is near zero
 *  \ingroup QuaternionOperations
 */
void normalizeQuaternion(Eigen::Ref<siconos::algebra::SiconosVector7> q);

/** \brief Compute the 3x3 rotation matrix from a quaternion
 *
 *  Computes the rotation matrix R using the direct formula:
 *  \f[ R = \begin{bmatrix}
 *  q_0^2+q_1^2-q_2^2-q_3^2 & 2(q_1q_2-q_0q_3) & 2(q_1q_3+q_0q_2) \\
 *  2(q_1q_2+q_0q_3) & q_0^2-q_1^2+q_2^2-q_3^2 & 2(q_2q_3-q_0q_1) \\
 *  2(q_1q_3-q_0q_2) & 2(q_2q_3+q_0q_1) & q_0^2-q_1^2-q_2^2+q_3^2
 *  \end{bmatrix} \f]
 *
 *  \param[in] q0 scalar part (w) of the quaternion
 *  \param[in] q1 first vector component (x) of the quaternion
 *  \param[in] q2 second vector component (y) of the quaternion
 *  \param[in] q3 third vector component (z) of the quaternion
 *  \param[out] rotationMatrix the 3x3 rotation matrix, modified in-place
 *  \ingroup QuaternionConversions
 */
void computeRotationMatrix(double q0, double q1, double q2, double q3,
                           siconos::algebra::SiconosMatrix33 &rotationMatrix);

/** \brief Compute rotation matrix from configuration
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[out] rotationMatrix the 3x3 rotation matrix, modified in-place
 *  \see computeRotationMatrix(double, double, double, double, SiconosMatrix33&)
 *  \ingroup QuaternionConversions
 */
void computeRotationMatrix(const siconos::algebra::SiconosVector7 &q,
                           siconos::algebra::SiconosMatrix33 &rotationMatrix);

/** \brief Compute transposed rotation matrix from configuration
 *
 *  Equivalent to the inverse rotation (R^T = R^-1 for rotation matrices).
 *  Used to transform from inertial to body frame.
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[out] rotationMatrix the transposed 3x3 rotation matrix, modified in-place
 *  \ingroup QuaternionConversions
 */
void computeRotationMatrixTransposed(const siconos::algebra::SiconosVector7 &q,
                                     siconos::algebra::SiconosMatrix33 &rotationMatrix);

/** \brief Rotate a vector using quaternion (inertial frame rotation)
 *
 *  Rotates vector v by the quaternion in q: v' = R(q) * v
 *  where R(q) is the rotation matrix corresponding to q.
 *
 *  Uses an efficient implementation with cross products:
 *  t = 2 * q_vec × v
 *  v' = v + q_vec × t + q_w * t
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[in,out] v the vector to rotate, modified in-place
 *  \ingroup QuaternionOperations
 */
void rotateVector(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                  Eigen::Ref<siconos::algebra::SiconosVector3> v);

/** \brief Rotate a matrix using quaternion (inertial frame rotation)
 *
 *  Rotates each column of matrix m by the quaternion in q.
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[in,out] m the matrix to rotate (3x3), modified in-place
 *  \ingroup QuaternionOperations
 */
void rotateMatrix(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                  Eigen::Ref<siconos::algebra::SiconosMatrix33> m);

/** \brief Transform vector from inertial frame to body frame
 *
 *  Applies the inverse rotation: v_body = R(q)^T * v_inertial
 *  This transforms a vector expressed in the inertial frame
 *  to the body (local) frame.
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[in,out] v the vector to transform, modified in-place
 *  \see rotateVectorFromBodyToInertialFrame for the inverse operation
 *  \ingroup FrameTransformations
 */
void rotateVectorFromInertialToBodyFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosVector3> v);

/** \brief Transform matrix from inertial frame to body frame
 *
 *  Applies the inverse rotation to each column: M_body = R(q)^T * M_inertial
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[in,out] m the matrix to transform, modified in-place
 *  \see rotateMatrixFromBodyToInertialFrame for the inverse operation
 *  \ingroup FrameTransformations
 */
void rotateMatrixFromInertialToBodyFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosMatrix33> m);

/** \brief Transform vector from body frame to inertial frame
 *
 *  Applies the forward rotation: v_inertial = R(q) * v_body
 *  This transforms a vector expressed in the body (local) frame
 *  to the inertial frame.
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[in,out] v the vector to transform, modified in-place
 *  \see rotateVectorFromInertialToBodyFrame for the inverse operation
 *  \ingroup FrameTransformations
 */
void rotateVectorFromBodyToInertialFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosVector3> v);

/** \brief Transform matrix from body frame to inertial frame
 *
 *  Applies the forward rotation to each column: M_inertial = R(q) * M_body
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[in,out] m the matrix to transform, modified in-place
 *  \see rotateMatrixFromInertialToBodyFrame for the inverse operation
 *  \ingroup FrameTransformations
 */
void rotateMatrixFromBodyToInertialFrame(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                                         Eigen::Ref<siconos::algebra::SiconosMatrix33> m);

/** \brief Compose two rigid body configurations (SE(3) group operation)
 *
 *  Computes ab = a ∘ b where ∘ is the composition law for rigid body motions:
 *  - Position: ab_pos = a_pos + R(a_rot) * b_pos
 *  - Rotation: ab_rot = a_rot * b_rot (quaternion product)
 *
 *  \param[in] a first configuration [x, y, z, qw, qx, qy, qz]
 *  \param[in] b second configuration [x, y, z, qw, qx, qy, qz]
 *  \param[out] ab result configuration, modified in-place
 *  \ingroup QuaternionOperations
 */
void compositionLawLieGroup(const siconos::algebra::SiconosVector7 &a,
                            siconos::algebra::SiconosVector7 &b,
                            siconos::algebra::SiconosVector7 &ab);

/** \brief In-place composition of rigid body configurations
 *
 *  Computes b = a ∘ b (modifies b in-place).
 *
 *  \param[in] a first configuration
 *  \param[in,out] b second configuration, replaced by composition result
 *  \see compositionLawLieGroup
 *  \ingroup QuaternionOperations
 */
void compositionLawLieGroup(const siconos::algebra::SiconosVector7 &a,
                            Eigen::Ref<siconos::algebra::SiconosVector7> b);

/** \brief Copy rotation quaternion from configuration
 *
 *  \param[in] from configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[out] to boost quaternion [qw, qx, qy, qz]
 *  \ingroup QuaternionAccessors
 */
void extractRotationQuaternion(const siconos::algebra::SiconosVector7 &from,
                               boost::math::quaternion<double> &to);

/** \brief Copy vector part from pure quaternion to 3D vector
 *
 *  \param[in] from boost quaternion [w, x, y, z]
 *  \param[out] to 3D vector [x, y, z]
 *  \ingroup QuaternionAccessors
 */
void extractVectorFromQuaternion(const boost::math::quaternion<double> &from,
                                 siconos::algebra::SiconosVector3 &to);

/** \brief Copy position from configuration to pure quaternion
 *
 *  \param[in] from configuration vector [x, y, z, qw, qx, qy, qz]
 *  \param[out] to boost quaternion [0, x, y, z]
 *  \ingroup QuaternionAccessors
 */
void extractPositionToQuaternion(const siconos::algebra::SiconosVector7 &from,
                                 boost::math::quaternion<double> &to);

/** \brief Copy 2D rotation to quaternion
 *
 *  For 2D rotations around z-axis, converts angle to quaternion [cos(θ/2), 0, 0, sin(θ/2)]
 *
 *  \param[in] from 2D configuration vector where from(2) is the angle
 *  \param[out] to boost quaternion
 *  \ingroup QuaternionAccessors
 */
void extractRotationQuaternion2d(const siconos::algebra::SiconosVector &from,
                                 boost::math::quaternion<double> &to);

/** \brief Copy 2D position from quaternion
 *
 *  \param[in] from boost quaternion
 *  \param[out] to 2D vector where to(0)=x, to(1)=y
 *  \ingroup QuaternionAccessors
 */
void extractPositionFromQuaternion2d(const boost::math::quaternion<double> &from,
                                     siconos::algebra::SiconosVector &to);

/** \brief Copy 2D position to quaternion
 *
 *  \param[in] from 2D configuration vector
 *  \param[out] to boost quaternion [0, x, y, 0]
 *  \ingroup QuaternionAccessors
 */
void extractPositionToQuaternion2d(const siconos::algebra::SiconosVector &from,
                                   boost::math::quaternion<double> &to);

/** \brief Get rotation quaternion from configuration
 *
 *  \param[in] q configuration vector [x, y, z, qw, qx, qy, qz]
 *  \return boost quaternion [qw, qx, qy, qz]
 *  \ingroup QuaternionAccessors
 */
boost::math::quaternion<double> getRotationQuaternion(const siconos::algebra::SiconosVector7 &q);

/** \brief Get position quaternion from configuration
 *
 *  \param[in] v configuration or position vector
 *  \return boost quaternion [0, x, y, z]
 *  \ingroup QuaternionAccessors
 */
boost::math::quaternion<double> getPositionQuaternion(const siconos::algebra::SiconosVector &v);

/** \brief Compute an orthonormal basis from a given input axis (Gram-Schmidt)
 *
 *  Given a primary axis, computes two orthogonal vectors using Gram-Schmidt
 *  orthogonalization to form a right-handed orthonormal basis.
 *
 *  \param[in,out] axis0 primary axis (will be normalized), becomes first basis vector
 *  \param[out] axis1 second basis vector (orthogonal to axis0)
 *  \param[out] axis2 third basis vector (orthogonal to axis0 and axis1)
 *  \throw std::invalid_argument if axis0 has near-zero norm
 *  \ingroup QuaternionOperations
 */
void computeOrthonormalBaseFromAxis(siconos::algebra::SiconosVector3 &axis0,
                                    siconos::algebra::SiconosVector3 &axis1,
                                    siconos::algebra::SiconosVector3 &axis2);

/** \brief Compute an orthonormal basis from a vector (branchless method)
 *
 *  Implements the branchless orthonormal basis algorithm from:
 *  "Building an Orthonormal Basis, Revisited" (Duff et al., JCGT 2017)
 *
 *  This method is numerically stable and avoids branching based on the
 *  input vector's orientation.
 *
 *  \param[in,out] A reference vector (normalized after call), becomes first basis vector
 *  \param[out] A1 second basis vector (orthogonal to A)
 *  \param[out] A2 third basis vector (orthogonal to A and A1)
 *  \return true if successful, false if A has zero norm (outputs set to NaN)
 *  \ingroup QuaternionOperations
 *
 *  \note Reference:
 *  Tom Duff et al., "Building an Orthonormal Basis, Revisited",
 *  Journal of Computer Graphics Techniques, vol. 6, no. 1, 2017
 *  http://jcgt.org/published/0006/01/01/
 */
bool orthoBaseFromVector(siconos::algebra::SiconosVector3 &A,
                         siconos::algebra::SiconosVector3 &A1,
                         siconos::algebra::SiconosVector3 &A2);

/** @name Deprecated Function Names
 *  Old function names kept for backward compatibility.
 *  \deprecated Use the new function names instead.
 */
//@{
[[deprecated("Use rotateVector instead")]] inline void quaternionRotateVector(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  rotateVector(q, v);
}

[[deprecated("Use rotateMatrix instead")]] inline void quaternionRotateMatrix(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  rotateMatrix(q, m);
}

[[deprecated("Use rotateVectorFromInertialToBodyFrame instead")]] inline void rewriteVectorFromAbsoluteToBodyFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  rotateVectorFromInertialToBodyFrame(q, v);
}

[[deprecated("Use rotateMatrixFromInertialToBodyFrame instead")]] inline void rewriteMatrixFromAbsoluteToBodyFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  rotateMatrixFromInertialToBodyFrame(q, m);
}

[[deprecated("Use rotateVectorFromBodyToInertialFrame instead")]] inline void rewriteVectorFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v) {
  rotateVectorFromBodyToInertialFrame(q, v);
}

[[deprecated("Use rotateMatrixFromBodyToInertialFrame instead")]] inline void rewriteMatrixFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m) {
  rotateMatrixFromBodyToInertialFrame(q, m);
}

[[deprecated("Use normalizeQuaternion instead")]] inline void normalizeq(Eigen::Ref<siconos::algebra::SiconosVector7> q) {
  normalizeQuaternion(q);
}

[[deprecated("Use extractRotationQuaternion instead")]] inline void copyQuatRot(
    const siconos::algebra::SiconosVector7 &from,
    boost::math::quaternion<double> &to) {
  extractRotationQuaternion(from, to);
}

[[deprecated("Use extractPositionToQuaternion instead")]] inline void copyQuatPos(
    const siconos::algebra::SiconosVector7 &from,
    boost::math::quaternion<double> &to) {
  extractPositionToQuaternion(from, to);
}

[[deprecated("Use extractRotationQuaternion2d instead")]] inline void copyQuatRot2d(
    const siconos::algebra::SiconosVector &from,
    boost::math::quaternion<double> &to) {
  extractRotationQuaternion2d(from, to);
}

[[deprecated("Use extractPositionToQuaternion2d instead")]] inline void copyQuatPos2d(
    const siconos::algebra::SiconosVector &from,
    boost::math::quaternion<double> &to) {
  extractPositionToQuaternion2d(from, to);
}

[[deprecated("Use getRotationQuaternion instead")]] inline boost::math::quaternion<double> rotquat(
    const siconos::algebra::SiconosVector7 &v) {
  return getRotationQuaternion(v);
}

[[deprecated("Use getPositionQuaternion instead")]] inline boost::math::quaternion<double> posquat(
    const siconos::algebra::SiconosVector &v) {
  return getPositionQuaternion(v);
}

[[deprecated("Use sinc instead")]] inline double sin_x(double x) {
  return sinc(x);
}
//@}

}  // namespace siconos::geometry
#endif  // ROTATIONQUATERNION_H
