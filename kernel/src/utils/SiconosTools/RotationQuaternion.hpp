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
                               Eigen::Ref<siconos::algebra::SiconosVector> &axis);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the angle/axis representation
 */
double axisAngleFromConfiguration(const Eigen::Ref<siconos::algebra::SiconosVector> &q,
                                  Eigen::Ref<siconos::algebra::SiconosVector> axis);

/* For a given quaternion, compute the rotation vector representation
 */
void rotationVectorFromQuaternion(double q0, double q1, double q2, double q3,
                                  siconos::algebra::SiconosVector &rotationVector);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the rotation vector representation
 */
void rotationVectorFromConfiguration(siconos::algebra::SiconosVector &q,
                                     siconos::algebra::SiconosVector &rotationVector);

/* For a given angle and rotation vector, compute the unit quaternion
 */
void quaternionFromAxisAngle(const siconos::algebra::SiconosVector &axis, double angle,
                             siconos::algebra::SiconosVector &q);

/* For a given  rotation vector, compute the quaternion
 */
void quaternionFromRotationVector(const siconos::algebra::SiconosVector &rotationVector,
                                  siconos::algebra::SiconosVector &q);

double sin_x(double x);

void quaternionFromTwistVector(const siconos::algebra::SiconosVector &twist,
                               siconos::algebra::SiconosVector &q);

/* For a given quaternion q, compute the norm
 */
double quaternionNorm(const siconos::algebra::SiconosVector &q);

/* For a given quaternion q, compute the unit quaternion by normalization
 */
void normalizeq(siconos::algebra::SiconosVector &q);

/* For a given quaternion q, compute the associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */

void computeRotationMatrix(double q0, double q1, double q2, double q3,
                           siconos::algebra::SiconosMatrix &rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */

void computeRotationMatrix(const siconos::algebra::SiconosVector &q,
                           siconos::algebra::SiconosMatrix &rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the transposed associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void computeRotationMatrixTransposed(const siconos::algebra::SiconosVector &q,
                                     siconos::algebra::SiconosMatrix &rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 *  performs the rotation of the vector v
 *  w.r.t the quaternion that parametrize the rotation in q
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void quaternionRotateVector(const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
                            Eigen::Ref<siconos::algebra::SiconosVector3> v);

/* For a given configuration vector q composed of a position and a quaternion,
 * performs the rotation of the matrix m
 * w.r.t the quaternion that parametrize the rotation in q
 * \param[in] q the position vector
 * \param[in,out] m the vector to be rotated
 */
void quaternionRotateMatrix(const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
                            Eigen::Ref<siconos::algebra::SiconosMatrix> m);

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
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
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
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix> m);

void rewriteVectorFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
    Eigen::Ref<siconos::algebra::SiconosVector> v);

void rewriteMatrixFromBodyToAbsoluteFrame(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
    Eigen::Ref<siconos::algebra::SiconosMatrix> m);

void compositionLawLieGroup(const siconos::algebra::SiconosVector &a,
                            siconos::algebra::SiconosVector &b,
                            siconos::algebra::SiconosVector &ab);

void compositionLawLieGroup(const siconos::algebra::SiconosVector &a,
                            siconos::algebra::SiconosVector &b);

void copyQuatRot(const siconos::algebra::SiconosVector &from,
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

boost::math::quaternion<double> rotquat(const siconos::algebra::SiconosVector &v);

boost::math::quaternion<double> posquat(const siconos::algebra::SiconosVector &v);

}  // namespace siconos::geometry
#endif  // ROTATIONQUATERNION_H
