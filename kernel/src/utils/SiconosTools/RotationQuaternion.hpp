/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/* For a given quaternion q, compute the angle/axis representation
 */
double axisAngleFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector axis);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the angle/axis representation
 */
double axisAngleFromConfiguration(SP::SiconosVector q, SP::SiconosVector axis);

/* For a given quaternion, compute the rotation vector representation
 */
void rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector rotationVector);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the rotation vector representation
 */
void rotationVectorFromConfiguration(SP::SiconosVector q, SP::SiconosVector rotationVector);

/* For a given angle and rotation vector, compute the unit quaternion
 */
void quaternionFromAxisAngle(SP::SiconosVector axis, double angle, SP::SiconosVector q);

/* For a given  rotation vector, compute the quaternion
 */
void quaternionFromRotationVector(SP::SiconosVector rotationVector, SP::SiconosVector q);
/* For a given quaternion q, compute the unit quaternion by normalization
 */
void normalizeq(SP::SiconosVector q);

/* For a given quaternion q, compute the associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */

void computeRotationMatrix(double q0, double q1, double q2, double q3, SP::SimpleMatrix rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */

void computeRotationMatrix(SP::SiconosVector q,  SP::SimpleMatrix rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 * compute the transposed associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void computeRotationMatrixTransposed(SP::SiconosVector q, SP::SimpleMatrix rotationMatrix);

/* For a given configuration vector q composed of a position and a quaternion,
 *  performs the rotation of the vector v
 * w.r.t the quaternion that parametrize the rotation in q
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void quaternionRotate(SP::SiconosVector q, SP::SiconosVector v);

/* For a given quaternion q, compute the associated rotation matrix
 * w.r.t the quaternion that parametrize the rotation in q,
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void quaternionRotate(double q0, double q1, double q2, double q3, SiconosVector& v);
void quaternionRotate(double q0, double q1, double q2, double q3, SP::SiconosVector v);
void quaternionRotate(double q0, double q1, double q2, double q3, SP::SimpleMatrix m);

/* For a given configuration vector q composed of a position and a quaternion,
 *  performs the rotation of the vector v
 * w.r.t the quaternion that parametrize the rotation in q
 * \param[in] q the position vector
 * \param[in,out] v the vector to be rotated
 */
void quaternionRotate(SP::SiconosVector q, SP::SiconosVector v);

/* For a given configuration vector q composed of a position and a quaternion,
 * performs the rotation of the matrix m
 * w.r.t the quaternion that parametrize the rotation in q
 * \param[in] q the position vector
 * \param[in,out] m the vector to be rotated
 */
void quaternionRotate(SP::SiconosVector q, SP::SimpleMatrix m);



/* For a given  configuration vector q composed of a position and a quaternion,
 * express the vector v given in
 * the inertial frame into to the body frame
 * w.r.t the quaternion that parametrize the rotation in q.
 * The operation amounts to multiplying by the transposed rotation matrix.
 * the result is return in v
 * \param[in] q the position vector
 * \param[in,out] v the vector to be reexpressed
 */
void changeFrameAbsToBody(const SiconosVector& q, SiconosVector& v);
void changeFrameAbsToBody(SP::SiconosVector q, SP::SiconosVector v);
void changeFrameAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m);

void changeFrameBodyToAbs(const SiconosVector& q, SiconosVector& v);
void changeFrameBodyToAbs(SP::SiconosVector q, SP::SiconosVector v);
void changeFrameBodyToAbs(SP::SiconosVector q, SP::SimpleMatrix m);





#endif // ROTATIONQUATERNION_H
