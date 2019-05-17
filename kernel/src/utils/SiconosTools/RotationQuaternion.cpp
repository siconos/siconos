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
//#define DEBUG_NOCOLOR
//#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <debug.h>
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include <boost/math/quaternion.hpp>


void computeRotationMatrix(double q0, double q1, double q2, double q3,
                           SP::SimpleMatrix rotationMatrix)
{

  /* Brute force version by multiplication of quaternion
   */
  // ::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  // ::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  // ::boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  // ::boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  // ::boost::math::quaternion<double>    quatz(0, 0, 0, 1);
  // ::boost::math::quaternion<double>    quatBuff;
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
  rotationMatrix->setValue(0, 0,     q0*q0 +q1*q1 -q2*q2 -q3*q3);
  rotationMatrix->setValue(0, 1, 2.0*(q1*q2        - q0*q3));
  rotationMatrix->setValue(0, 2, 2.0*(q1*q3        + q0*q2));

  rotationMatrix->setValue(1, 0, 2.0*(q1*q2        + q0*q3));
  rotationMatrix->setValue(1, 1,     q0*q0 -q1*q1 +q2*q2 -q3*q3);
  rotationMatrix->setValue(1, 2, 2.0*(q2*q3        - q0*q1));

  rotationMatrix->setValue(2, 0, 2.0*(q1*q3        - q0*q2));
  rotationMatrix->setValue(2, 1, 2.0*(q2*q3         + q0*q1));
  rotationMatrix->setValue(2, 2,     q0*q0 -q1*q1 -q2*q2 +q3*q3);
}



void quaternionRotate(double q0, double q1, double q2, double q3, SiconosVector& v)
{
  DEBUG_BEGIN("::quaternionRotate(double q0, double q1, double q2, double q3, SiconosVector& v )\n");
  DEBUG_EXPR(v.display(););
  DEBUG_PRINTF("( q0 = %16.12e,  q1 = %16.12e,  q2= %16.12e,  q3= %16.12e )\n", q0,q1,q2,q3);
  assert(v.size()==3);

  // First way. Using the rotation matrix
  // SP::SimpleMatrix rotationMatrix(new SimpleMatrix(3,3));
  // SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,q1,q2,q3, rotationMatrix);
  // prod(*rotationMatrix, v, tmp);
  // v = tmp;
  // return;

  // Second way. Using the transpose of the rotation matrix
  // SP::SimpleMatrix rotationMatrix(new SimpleMatrix(3,3));
  // SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,-q1,-q2,-q3, rotationMatrix);
  // prod(v, *rotationMatrix, tmp);
  // v = tmp;

  // Third way. cross product and axis angle
  // see http://www.geometrictools.com/Documentation/RotationIssues.pdf
  // SP::SiconosVector axis(new SiconosVector(3));
  // double angle = ::axisAngleFromQuaternion(q0,q1,q2,q3, axis);
  // SiconosVector t(3), tmp(3);
  // cross_product(*axis,v,t);
  // cross_product(*axis,t,tmp);
  // v += sin(angle)*t + (1.0-cos(angle))*tmp;

  // Direct computation with cross product
  // Works only with unit quaternion
  SiconosVector t(3), tmp(3);
  SiconosVector qvect(3);
  qvect(0)=q1;
  qvect(1)=q2;
  qvect(2)=q3;
  cross_product(qvect,v,t);
  t *= 2.0;
  cross_product(qvect,t,tmp);
  v += tmp;
  v += q0*t;
  DEBUG_EXPR(v.display(););
  DEBUG_END("::quaternionRotate(double q0, double q1, double q2, double q3, SP::SiconosVector v )\n");
}

void quaternionRotate(double q0, double q1, double q2, double q3, SP::SiconosVector v)
{
  ::quaternionRotate(q0, q1, q2, q3, *v);
}

void quaternionRotate(double q0, double q1, double q2, double q3, SP::SimpleMatrix m)
{
  DEBUG_BEGIN("::quaternionRotate(double q0, double q1, double q2, double q3, SP::SimpleMatrix m )\n");
  DEBUG_EXPR(m->display(););
  DEBUG_PRINTF("( q0 = %16.12e,  q1 = %16.12e,  q2= %16.12e,  q3= %16.12e )\n", q0,q1,q2,q3);

  // Direct computation with cross product for each column
  assert(m->size(0) == 3 && "::quaternionRotate(double q0, double q1, double q2, double q3, SP::SimpleMatrix m ) m must have 3 rows");
  SiconosVector v(3);
  SiconosVector t(3), tmp(3);
  SiconosVector qvect(3);
  qvect(0)=q1;
  qvect(1)=q2;
  qvect(2)=q3;
  for(unsigned int j = 0; j < m->size(1); j++)
  {
    v(0) = m->getValue(0,j);
    v(1) = m->getValue(1,j);
    v(2) = m->getValue(2,j);
    cross_product(qvect,v,t);
    t *= 2.0;
    cross_product(qvect,t,tmp);
    v += tmp;
    v += q0*t;
    m->setValue(0,j,v(0));
    m->setValue(1,j,v(1));
    m->setValue(2,j,v(2));
  }
  DEBUG_EXPR(m->display(););
  DEBUG_END("::quaternionRotate(double q0, double q1, double q2, double q3, SP::SimpleMatrix m )\n");
}


void quaternionRotate(SP::SiconosVector q, SP::SiconosVector v)
{
  DEBUG_BEGIN("::quaternionRotate(SP::SiconosVector q, SP::SiconosVector v )\n");
  ::quaternionRotate(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6), v);
  DEBUG_END("::quaternionRotate(SP::SiconosVector q, SP::SiconosVector v )\n");
}

void quaternionRotate(SP::SiconosVector q, SP::SimpleMatrix m)
{
  DEBUG_BEGIN("::quaternionRotate(SP::SiconosVector q, SP::SimpleMatrix m )\n");
  ::quaternionRotate(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6),m);
  DEBUG_END("::quaternionRotate(SP::SiconosVector q, SP::SimpleMatrix m)\n");
}

void changeFrameAbsToBody(const SiconosVector& q, SiconosVector& v)
{
  DEBUG_BEGIN("::changeFrameAbsToBody(const SiconosVector& q, SiconosVector& v )\n");
  ::quaternionRotate(q.getValue(3),-q.getValue(4),-q.getValue(5),-q.getValue(6), v);
  DEBUG_END("::changeFrameAbsToBody(const SiconosVector& q, SiconosVector& v )\n");
}
void changeFrameAbsToBody(SP::SiconosVector q, SP::SiconosVector v)
{
  DEBUG_BEGIN("::changeFrameAbsToBody(SP::SiconosVector q, SP::SiconosVector v )\n");
  ::quaternionRotate(q->getValue(3),-q->getValue(4),-q->getValue(5),-q->getValue(6), v);
  DEBUG_END("::changeFrameAbsToBody(SP::SiconosVector q, SP::SiconosVector v )\n");
}
void changeFrameAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m)
{
  DEBUG_BEGIN("::changeFrameAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m )\n");
  ::quaternionRotate(q->getValue(3),-q->getValue(4),-q->getValue(5),-q->getValue(6), m);
  DEBUG_END("::changeFrameAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m )\n");
}

void changeFrameBodyToAbs(const SiconosVector& q, SiconosVector& v)
{
  DEBUG_BEGIN("::changeFrameBodyToAbs(const SiconosVector& q, SiconosVector& v )\n");
  ::quaternionRotate(q.getValue(3),q.getValue(4),q.getValue(5),q.getValue(6), v);
  DEBUG_END("::changeFrameBodyToAbs(const SiconosVector& q, SiconosVector& v )\n");
}
void changeFrameBodyToAbs(SP::SiconosVector q, SP::SiconosVector v)
{
  DEBUG_BEGIN("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SiconosVector v )\n");
  ::quaternionRotate(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6), *v);
  DEBUG_END("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SiconosVector v )\n");
}
void changeFrameBodyToAbs(SP::SiconosVector q, SP::SimpleMatrix m)
{
  DEBUG_BEGIN("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SimpleMatrix m )\n");
  ::quaternionRotate(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6), m);
  DEBUG_END("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SimpleMatrix m )\n");
}



void computeRotationMatrix(SP::SiconosVector q, SP::SimpleMatrix rotationMatrix)
{
  ::computeRotationMatrix(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6),
                          rotationMatrix);
}
void computeRotationMatrixTransposed(SP::SiconosVector q, SP::SimpleMatrix rotationMatrix)
{
  ::computeRotationMatrix(q->getValue(3),-q->getValue(4),-q->getValue(5),-q->getValue(6),
                          rotationMatrix);
}

double axisAngleFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector axis)
{
  DEBUG_BEGIN("axisAngleFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector axis )\n");
  double angle = acos(q0) *2.0;
  //double f = sin( angle *0.5);
  double f = sqrt(1-q0*q0); // cheaper than sin ?
  if(f !=0.0)
  {
    axis->setValue(0, q1/f);
    axis->setValue(1, q2/f);
    axis->setValue(2, q3/f);
  }
  else
  {
    axis->zero();
  }
  DEBUG_PRINTF("angle= %12.8e\n", angle);
  DEBUG_EXPR(axis->display(););
  DEBUG_END("axisAngleFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector axis )\n");
  return angle;
}

double axisAngleFromConfiguration(SP::SiconosVector q, SP::SiconosVector axis)
{
  double angle = ::axisAngleFromQuaternion(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6),axis);
  return angle;
}

void rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector rotationVector)
{
  DEBUG_BEGIN("rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector rotationVector )\n");

  rotationVector->setValue(0, q1);
  rotationVector->setValue(1, q2);
  rotationVector->setValue(2, q3);

  double norm_v = sqrt(q1*q1+q2*q2+q3*q3);
  assert(norm_v <= M_PI);  /* it should be called for a unit quaternion */
  if(norm_v < 1e-12)
  {
    rotationVector->setValue(0, 0.0);
    rotationVector->setValue(1, 0.0);
    rotationVector->setValue(2, 0.0);
  }
  else
  {
    *rotationVector *=  2.0 * asin(norm_v)/norm_v;
  }
  DEBUG_EXPR(rotationVector->display(););
  DEBUG_END("rotationVectorFromQuaternion(double q0, double q1, double q2, double q3, SP::SiconosVector rotationVector )\n");
}

void rotationVectorFromConfiguration(SP::SiconosVector q, SP::SiconosVector rotationVector)
{
  ::rotationVectorFromQuaternion(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6), rotationVector);
}


void quaternionFromAxisAngle(SP::SiconosVector axis, double angle, SP::SiconosVector q)
{
  q->setValue(3,cos(angle/2.0));
  q->setValue(4,axis->getValue(0)* sin(angle *0.5));
  q->setValue(5,axis->getValue(1)* sin(angle *0.5));
  q->setValue(6,axis->getValue(2)* sin(angle *0.5));
}

static
double sin_x(double x)
{
  if(std::abs(x) <= 1e-3)
  {
    return 1.0 + x*x / 3.0 + pow(x,4) * 2.0 / 15.0 + pow(x,6) * 17.0 / 315.0 + pow(x,8) * 62.0 / 2835.0;
  }
  else
  {
    return sin(x)/x;
  }
}

void quaternionFromRotationVector(SP::SiconosVector rotationVector, SP::SiconosVector q)
{
  double angle = sqrt(rotationVector->getValue(0)*rotationVector->getValue(0)+
               rotationVector->getValue(1)*rotationVector->getValue(1)+
               rotationVector->getValue(2)*rotationVector->getValue(2));

  double f = 0.5 * sin_x(angle *0.5);

  q->setValue(3,cos(angle/2.0));
  q->setValue(4,rotationVector->getValue(0)* f);
  q->setValue(5,rotationVector->getValue(1)* f);
  q->setValue(6,rotationVector->getValue(2)* f);
}


void normalizeq(SP::SiconosVector q)
{
  double normq = sqrt(q->getValue(3) * q->getValue(3) +
                      q->getValue(4) * q->getValue(4) +
                      q->getValue(5) * q->getValue(5) +
                      q->getValue(6) * q->getValue(6));
  assert(normq > 0);
  normq = 1.0 / normq;
  q->setValue(3, q->getValue(3) * normq);
  q->setValue(4, q->getValue(4) * normq);
  q->setValue(5, q->getValue(5) * normq);
  q->setValue(6, q->getValue(6) * normq);
}
