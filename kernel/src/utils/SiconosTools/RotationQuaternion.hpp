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

/** \file Rotationquaternion.hpp
 */

#ifndef ROTATIONQUATERNION_H
#define ROTATIONQUATERNION_H

#include <boost/math_fwd.hpp>  // for quaternion
#include <memory>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

// class SiconosVector;
// class SiconosMatrix;
}  // namespace siconos::algebra

namespace siconos::geometry {

/* For a given quaternion q, compute the angle/axis representation
 */
double axisAngleFromQuaternion(double q0, double q1, double q2, double q3,
                               Eigen::Ref<siconos::algebra::SiconosVector3> &axis);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the angle/axis representation
 */
double axisAngleFromConfiguration(const Eigen::Ref<siconos::algebra::SiconosVector7> &q,
                                  Eigen::Ref<siconos::algebra::SiconosVector3> axis);

/* For a given quaternion, compute the rotation vector representation
 */
siconos::algebra::SiconosVector3 rotationVectorFromQuaternion(double q0, double q1, double q2,
                                                              double q3);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the rotation vector representation
 */
siconos::algebra::SiconosVector3 rotationVectorFromConfiguration(
    siconos::algebra::SiconosVector7 &q);

/* For a given angle and rotation vector, compute the unit quaternion
 */
void quaternionFromAxisAngle(const siconos::algebra::SiconosVector3 &axis, double angle,
                             siconos::algebra::SiconosVector7 &q);

/* For a given  rotation vector, compute the quaternion
 */
siconos::algebra::SiconosVector7 quaternionFromRotationVector(
    const siconos::algebra::SiconosVector3 &rotationVector);

double sin_x(double x);

void quaternionFromTwistVector(const siconos::algebra::SiconosVector6 &twist,
                               Eigen::Ref<siconos::algebra::SiconosVector7> q);

/* For a given quaternion q, compute the norm
 */
double quaternionNorm(const siconos::algebra::SiconosVector7 &q);

/* For a given quaternion q, compute the unit quaternion by normalization
 */
void normalizeq(Eigen::Ref<siconos::algebra::SiconosVector7> q);

/* For a given quaternion q, compute the associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */

void computeRotationMatrix(double q0, double q1, double q2, double q3,
                           siconos::algebra::SiconosMatrix33 &rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */

void computeRotationMatrix(const siconos::algebra::SiconosVector7 &q,
                           siconos::algebra::SiconosMatrix33 &rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the transposed associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void computeRotationMatrixTransposed(const siconos::algebra::SiconosVector7 &q,
                                     siconos::algebra::SiconosMatrix33 &rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 *  performs the rotation of the vector v
 *  w.r.t the quaternion that parametrize the rotation in q
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void quaternionRotateVector(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                            Eigen::Ref<siconos::algebra::SiconosVector3> v);

/* For a given configuration vector q composed of a position and a quaternion,
 * performs the rotation of the matrix m
 * w.r.t the quaternion that parametrize the rotation in q
 * \param[in] q the position vector
 * \param[in,out] m the vector to be rotated
 */
void quaternionRotateMatrix(const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
                            Eigen::Ref<siconos::algebra::SiconosMatrix33> m);

/* For a given  configuration vector q composed of a position and a quaternion,
 * express the vector v given in
 * the inertial frame into to the body frame
 * w.r.t the quaternion that parametrize the rotation in q.
 * The operation amounts to multiplying by the transposed rotation matrix.
 * the result is return in v
 * \param[in] q the position vector
 * \param[in,out] v the vector to be reexpressed
 */
void rewriteVectorFromAbsoluteToBodyFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v);

/* For a given  configuration vector q composed of a position and a quaternion,
 * express the matrix m given in
 * the inertial frame into to the body frame
 * w.r.t the quaternion that parametrize the rotation in q.
 * The operation amounts to multiplying by the transposed rotation matrix.
 * the result is return in v
 * \param[in] q the position vector
 * \param[in,out] m the matrix to be reexpressed
 */
void rewriteMatrixFromAbsoluteToBodyFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m);

void rewriteVectorFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosVector3> v);

void rewriteMatrixFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector7> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> m);

void compositionLawLieGroup(const siconos::algebra::SiconosVector7 &a,
                            siconos::algebra::SiconosVector7 &b,
                            siconos::algebra::SiconosVector7 &ab);

void compositionLawLieGroup(const siconos::algebra::SiconosVector7 &a,
                            Eigen::Ref<siconos::algebra::SiconosVector7> b);

void copyQuatRot(const siconos::algebra::SiconosVector7 &from,
                 boost::math::quaternion<double> &to);

void copyQuatPos(const boost::math::quaternion<double> &from,
                 siconos::algebra::SiconosVector &to);

void copyQuatPos(const siconos::algebra::SiconosVector &from,
                 boost::math::quaternion<double> &to);

void copyQuatRot2d(const siconos::algebra::SiconosVector &from,
                   boost::math::quaternion<double> &to);

void copyQuatPos2d(const boost::math::quaternion<double> &from,
                   siconos::algebra::SiconosVector &to);

void copyQuatPos2d(const siconos::algebra::SiconosVector &from,
                   boost::math::quaternion<double> &to);

boost::math::quaternion<double> rotquat(const siconos::algebra::SiconosVector7 &q);

boost::math::quaternion<double> posquat(const siconos::algebra::SiconosVector &q);

/** Compute an orthonormal basis from a given input axis

    \param[in,out] axis0 reference axis (will be normalized)
    \param[inout] axis1 second axis of the base
    \param[inout] axis2 third axis of the base

*/
void computeOrthonormalBaseFromAxis(siconos::algebra::SiconosVector3 &axis0,
                                    siconos::algebra::SiconosVector3 &axis1,
                                    siconos::algebra::SiconosVector3 &axis2);

/* function to compute an orthonormal basis form a given vector
 * adapted with eigen vectors from
 * Tom Duff, James Burgess, Per Christensen, Christophe Hery, Andrew Kensler, Max Liani, and
 * Ryusuke Villemin, Building an Orthonormal Basis, Revisited, Journal of Computer Graphics
 * Techniques (JCGT), vol. 6, no. 1, 1-8, 2017 Available online
 * http://jcgt.org/published/0006/01/01/ void branchlessONB(const Vec3f &n, Vec3f &b1, Vec3f
 * &b2)
 * {
 *   float sign = copysignf(1.0f, n.z);
 *   const float a = -1.0f / (sign + n.z);
 *   const float b = n.x * n.y * a;
 *   b1 = Vec3f(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
 *   b2 = Vec3f(b, sign + n.y * n.y * a, -n.y);
 * }
 */
/**  compute an orthonormal basis
 *   \param[in,out] A a reference vector (normalized after call)
 *   \param[out] A1 second base vector
 *   \param[out] A2 third base vector
 *   \return true if all went right else false (e.g. if A.norm =0)
 */
bool orthoBaseFromVector(siconos::algebra::SiconosVector3 &A,
                         siconos::algebra::SiconosVector3 &A1,
                         siconos::algebra::SiconosVector3 &A2);

}  // namespace siconos::geometry
#endif  // ROTATIONQUATERNION_H
