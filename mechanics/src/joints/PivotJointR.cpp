/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <Interaction.hpp>
#include <NewtonEulerDS.hpp>
#include <boost/math/quaternion.hpp>

#include <BlockVector.hpp>

#include <iostream>

#include <op3x3.h>

/*
 * This file contains some code generated using sympy.  The following
 * is the necessary predule:
 *
 * from sympy import Symbol
 * import numpy as np
 *
 * A1 = np.array([0, Symbol('_A1x'), Symbol('_A1y'), Symbol('_A1z')])
 * A2 = np.array([0, Symbol('_A2x'), Symbol('_A2y'), Symbol('_A2z')])
 * q1 = np.array([Symbol('q10'), Symbol('q11'), Symbol('q12'), Symbol('q13')])
 * q2 = np.array([Symbol('q20'), Symbol('q21'), Symbol('q22'), Symbol('q23')])
 * cq2q10 = np.array([Symbol('_cq2q101'),Symbol('_cq2q102'),
 *                    Symbol('_cq2q103'),Symbol('_cq2q104')])
 *
 * qinv = lambda q: np.array([q[0],-q[1],-q[2],-q[3]])
 * qmul = lambda a,b: np.array([
 *          a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
 *          a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
 *          a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
 *          a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]])
 */

// Wrap value in interval [-pi,pi]
static double piwrap(double x)
{
  return fmod(x + 3*M_PI, 2*M_PI) - M_PI;
}

PivotJointR::PivotJointR(SP::SiconosVector P, SP::SiconosVector A, bool absoluteRef,
                         SP::NewtonEulerDS d1, SP::NewtonEulerDS d2)
  : KneeJointR(P, absoluteRef, d1, d2)
  , _A(std11::make_shared<SiconosVector>(3))
{
  _axes.resize(1);

  setAxis(0, A);
  if (d1)
    setInitialConditions(d1->q(), d2 ? d2->q() : SP::SiconosVector());
}

static ::boost::math::quaternion<double> quat(const SP::SiconosVector& v)
{
  if (v && v->size()==7)
    return ::boost::math::quaternion<double>((*v)(3),(*v)(4),(*v)(5),(*v)(6));
  else if (v && v->size()==3)
    return ::boost::math::quaternion<double>(0, (*v)(0),(*v)(1),(*v)(2));
  else
    return ::boost::math::quaternion<double>(1, 0, 0, 0);
}

void PivotJointR::setInitialConditions(SP::SiconosVector q1,
                                       SP::SiconosVector q2)
{
  *_A = *_axes[0];
  // TODO: add support for absolute frame here

  ::boost::math::quaternion<double> quat1(quat(q1));
  ::boost::math::quaternion<double> quat2(quat(q2));
  ::boost::math::quaternion<double> cq2q10(1.0 / quat2 * quat1);

  _cq2q101 = cq2q10.R_component_1();
  _cq2q102 = cq2q10.R_component_2();
  _cq2q103 = cq2q10.R_component_3();
  _cq2q104 = cq2q10.R_component_4();

  buildA1A2();

  double rot2to1w, rot2to1x, rot2to1y, rot2to1z;
  if (q2)
    rot2to1((*q1)(3), (*q1)(4), (*q1)(5), (*q1)(6),
            (*q2)(3), (*q2)(4), (*q2)(5), (*q2)(6),
            &rot2to1w, &rot2to1x, &rot2to1y, &rot2to1z);
  else
    rot2to1((*q1)(3), (*q1)(4), (*q1)(5), (*q1)(6),
            1, 0, 0, 0,
            &rot2to1w, &rot2to1x, &rot2to1y, &rot2to1z);

  _initial_AscalA1 = AscalA1(rot2to1x, rot2to1y, rot2to1z);
  _initial_AscalA2 = AscalA2(rot2to1x, rot2to1y, rot2to1z);

  // In case of joint constraints, it's okay to use dot product=0, but
  // in the case of the free axis we must measure the actual angle
  // using atan2 so that stops can be placed correctly.
  double Adot2to1 = AscalA(rot2to1x, rot2to1y, rot2to1z);
  _initial_AscalA = 2*atan2(rot2to1w, Adot2to1);

  _twistCount = 0;
  _previousAngle = 0;
}

void PivotJointR::buildA1A2()
{
  double Ax = (*_A)(0);
  double Ay = (*_A)(1);
  double Az = (*_A)(2);
  if (orthoBaseFromVector(&Ax, &Ay, &Az,
                          &_A1x, &_A1y, &_A1z,
                          &_A2x, &_A2y, &_A2z))
    RuntimeException::selfThrow("PivotJointR::initComponents. Problem in calling orthoBaseFromVector");

  assert(fabs(_A1x * Ax + _A1y * Ay + _A1z * Az) < 1e-9 && "PivotJoint, _A1 wrong\n");
  assert(fabs(_A2x * Ax + _A2y * Ay + _A2z * Az) < 1e-9 && "PivotJoint, _A2 wrong\n");
  assert(fabs(_A1x * _A2x + _A1y * _A2y + _A1z * _A2z) < 1e-9 && "PivotJoint, _A12 wrong\n");
  // std::cout << "JointPivot: _A1x _A1y _A1z :" << _A1x << " " << _A1y << " " << _A1z << std::endl;
  // std::cout << "JointPivot: _A2x _A2y _A2z :" << _A2x << " " << _A2y << " " << _A2z << std::endl;
}
void PivotJointR::Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  KneeJointR::Jd1d2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);

  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);

  // sympy expression: [AscalA1.diff(x) for x in q1]
  _jachq->setValue(3, 3,
                   _A1x*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
                   + _A1y*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
                   + _A1z*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20));
  _jachq->setValue(3, 4,
                   _A1x*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
                   + _A1y*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
                   + _A1z*(_cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21));
  _jachq->setValue(3, 5,
                   _A1x*(_cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20)
                   + _A1y*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
                   + _A1z*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22));
  _jachq->setValue(3, 6,
                   _A1x*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
                   + _A1y*(_cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22)
                   + _A1z*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23));

  _jachq->setValue(3, 7, 0);
  _jachq->setValue(3, 8, 0);
  _jachq->setValue(3, 9, 0);

  // sympy expression: [AscalA1.diff(x) for x in q2]
  _jachq->setValue(3, 10,
                   _A1x*(_cq2q101*q11 - _cq2q102*q10 - _cq2q103*q13 + _cq2q104*q12)
                   + _A1y*(_cq2q101*q12 + _cq2q102*q13 - _cq2q103*q10 - _cq2q104*q11)
                   + _A1z*(_cq2q101*q13 - _cq2q102*q12 + _cq2q103*q11 - _cq2q104*q10));
  _jachq->setValue(3, 11,
                   _A1x*(-_cq2q101*q10 - _cq2q102*q11 + _cq2q103*q12 + _cq2q104*q13)
                   + _A1y*(_cq2q101*q13 - _cq2q102*q12 - _cq2q103*q11 + _cq2q104*q10)
                   + _A1z*(-_cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 - _cq2q104*q11));
  _jachq->setValue(3, 12,
                   _A1x*(-_cq2q101*q13 - _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10)
                   + _A1y*(-_cq2q101*q10 + _cq2q102*q11 - _cq2q103*q12 + _cq2q104*q13)
                   + _A1z*(_cq2q101*q11 + _cq2q102*q10 - _cq2q103*q13 - _cq2q104*q12));
  _jachq->setValue(3, 13,
                   _A1x*(_cq2q101*q12 - _cq2q102*q13 + _cq2q103*q10 - _cq2q104*q11)
                   + _A1y*(-_cq2q101*q11 - _cq2q102*q10 - _cq2q103*q13 - _cq2q104*q12)
                   + _A1z*(-_cq2q101*q10 + _cq2q102*q11 + _cq2q103*q12 - _cq2q104*q13));

  _jachq->setValue(4, 0, 0);
  _jachq->setValue(4, 1, 0);
  _jachq->setValue(4, 2, 0);

  // sympy expression: [AscalA2.diff(x) for x in q1]
  _jachq->setValue(4, 3,
                   _A2x*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
                   + _A2y*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
                   + _A2z*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20));
  _jachq->setValue(4, 4,
                   _A2x*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
                   + _A2y*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
                   + _A2z*(_cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21));
  _jachq->setValue(4, 5,
                   _A2x*(_cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20)
                   + _A2y*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
                   + _A2z*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22));
  _jachq->setValue(4, 6,
                   _A2x*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
                   + _A2y*(_cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22)
                   + _A2z*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23));

  _jachq->setValue(4, 7, 0);
  _jachq->setValue(4, 8, 0);
  _jachq->setValue(4, 9, 0);

  // sympy expression: [AscalA2.diff(x) for x in q1]
  _jachq->setValue(4, 10,
                   _A2x*(_cq2q101*q11 - _cq2q102*q10 - _cq2q103*q13 + _cq2q104*q12)
                   + _A2y*(_cq2q101*q12 + _cq2q102*q13 - _cq2q103*q10 - _cq2q104*q11)
                   + _A2z*(_cq2q101*q13 - _cq2q102*q12 + _cq2q103*q11 - _cq2q104*q10));
  _jachq->setValue(4, 11,
                   _A2x*(-_cq2q101*q10 - _cq2q102*q11 + _cq2q103*q12 + _cq2q104*q13)
                   + _A2y*(_cq2q101*q13 - _cq2q102*q12 - _cq2q103*q11 + _cq2q104*q10)
                   + _A2z*(-_cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 - _cq2q104*q11));
  _jachq->setValue(4, 12,
                   _A2x*(-_cq2q101*q13 - _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10)
                   + _A2y*(-_cq2q101*q10 + _cq2q102*q11 - _cq2q103*q12 + _cq2q104*q13)
                   + _A2z*(_cq2q101*q11 + _cq2q102*q10 - _cq2q103*q13 - _cq2q104*q12));
  _jachq->setValue(4, 13,
                   _A2x*(_cq2q101*q12 - _cq2q102*q13 + _cq2q103*q10 - _cq2q104*q11)
                   + _A2y*(-_cq2q101*q11 - _cq2q102*q10 - _cq2q103*q13 - _cq2q104*q12)
                   + _A2z*(-_cq2q101*q10 + _cq2q102*q11 + _cq2q103*q12 - _cq2q104*q13));

  /*proj_with_q
  for (unsigned int ii=0; ii <_jachq->size(0); ii++)
    for (unsigned int jj=0; jj <_jachq->size(1); jj++)
  _jachqProj->setValue(ii,jj,_jachq->getValue(ii,jj));

  _jachqProj->setValue(5,0,0);
  _jachqProj->setValue(5,1,0);
  _jachqProj->setValue(5,2,0);
  _jachqProj->setValue(5,3,2.0*q10);
  _jachqProj->setValue(5,4,2.0*q11);
  _jachqProj->setValue(5,5,2.0*q12);
  _jachqProj->setValue(5,6,2.0*q13);
  _jachqProj->setValue(6,0+7,0);
  _jachqProj->setValue(6,1+7,0);
  _jachqProj->setValue(6,2+7,0);
  _jachqProj->setValue(6,3+7,2.0*q20);
  _jachqProj->setValue(6,4+7,2.0*q21);
  _jachqProj->setValue(6,5+7,2.0*q22);
  _jachqProj->setValue(6,6+7,2.0*q23);
  */

  //_jachq->display();
}

void PivotJointR::Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13)
{

  KneeJointR::Jd1(X1, Y1, Z1, q10, q11, q12, q13);


  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);

  // sympy expression: [AscalA1.diff(x) for x in q1]
  _jachq->setValue(3, 3,
                   _A1x*(- _cq2q102)
                   + _A1y*(- _cq2q103)
                   + _A1z*(- _cq2q104));
  _jachq->setValue(3, 4,
                   _A1x*(_cq2q101)
                   + _A1y*(- _cq2q104)
                   + _A1z*(_cq2q103));
  _jachq->setValue(3, 5,
                   _A1x*(_cq2q104)
                   + _A1y*(_cq2q101)
                   + _A1z*(- _cq2q102));
  _jachq->setValue(3, 6,
                   _A1x*(- _cq2q103)
                   + _A1y*(_cq2q102)
                   + _A1z*(_cq2q101));


  _jachq->setValue(4, 0, 0);
  _jachq->setValue(4, 1, 0);
  _jachq->setValue(4, 2, 0);

  // sympy expression: [AscalA2.diff(x) for x in q1]
  _jachq->setValue(4, 3,
                   _A2x*(- _cq2q102)
                   + _A2y*(- _cq2q103)
                   + _A2z*(- _cq2q104));
  _jachq->setValue(4, 4,
                   _A2x*(_cq2q101)
                   + _A2y*(- _cq2q104)
                   + _A2z*(_cq2q103));
  _jachq->setValue(4, 5,
                   _A2x*(_cq2q104)
                   + _A2y*(_cq2q101)
                   + _A2z*(- _cq2q102));
  _jachq->setValue(4, 6,
                   _A2x*(- _cq2q103)
                   + _A2y*(_cq2q102)
                   + _A2z*(_cq2q101));

  /*proj_with_q
      for (unsigned int ii=0; ii <_jachq->size(0); ii++)
        for (unsigned int jj=0; jj <_jachq->size(1); jj++)
    _jachqProj->setValue(ii,jj,_jachq->getValue(ii,jj));

      _jachqProj->setValue(5,0,0);
      _jachqProj->setValue(5,1,0);
      _jachqProj->setValue(5,2,0);
      _jachqProj->setValue(5,3,2.0*q10);
      _jachqProj->setValue(5,4,2.0*q11);
      _jachqProj->setValue(5,5,2.0*q12);
      _jachqProj->setValue(5,6,2.0*q13);
  */

}

void PivotJointR::rot2to1(double q10, double q11, double q12, double q13,
                          double q20, double q21, double q22, double q23,
                          double *rot2to1w, double *rot2to1x,
                          double *rot2to1y, double *rot2to1z)
{
  /*
   * The current rotation vector taking into account initial rotation
   * difference.
   *
   * sympy expression:
   * rot2to1 = qmul(qinv(qmul(q2,cq2q10)),q1)
   */

  if (rot2to1w)
  *rot2to1w = (q10*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
               - q11*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
               - q12*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
               - q13*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20));
  *rot2to1x = (q10*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
               + q11*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
               - q12*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
               + q13*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21));
  *rot2to1y = (q10*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
               + q11*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
               + q12*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
               - q13*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22));
  *rot2to1z = (q10*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
               - q11*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
               + q12*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
               + q13*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23));
}

double PivotJointR::AscalA1(double rot2to1x, double rot2to1y, double rot2to1z)
{
  /*
   * The angle between A1 and rotation q2-to-q1 must be zero,
   * (taking into account original difference in orientation rot2to1).
   *
   * sympy expression:
   * AscalA1 = np.dot(A1,rot2to1) - initial_AscalA1
   */

  return _A1x*rot2to1x + _A1y*rot2to1y + _A1z*rot2to1z;
}

double PivotJointR::AscalA2(double rot2to1x, double rot2to1y, double rot2to1z)
{
  /*
   * The angle between A2 and rotation q2-to-q1 must be zero,
   * (taking into account original difference in orientation rot2to1).
   *
   * sympy expression:
   * AscalA2 = np.dot(A2,rot2to1) - initial_AscalA2
   */

  return _A2x*rot2to1x + _A2y*rot2to1y + _A2z*rot2to1z;
}

double PivotJointR::AscalA(double rot2to1x, double rot2to1y, double rot2to1z)
{
  /*
   * The angle between A and rotation q2-to-q1 must be zero,
   * (taking into account original difference in orientation rot2to1).
   *
   * sympy expression:
   * AscalA = np.dot(A,rot2to1) - initial_AscalA
   */

  return _A->getValue(0)*rot2to1x
    + _A->getValue(1)*rot2to1y
    + _A->getValue(2)*rot2to1z;
}

void PivotJointR::computeh(double time, BlockVector& q0, SiconosVector& y)
{

  KneeJointR::computeh(time, q0,  y);

  double q10 = q0.getValue(3);
  double q11 = q0.getValue(4);
  double q12 = q0.getValue(5);
  double q13 = q0.getValue(6);
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;
  if(q0.numberOfBlocks()>1)
  {
    q20 = q0.getValue(10);
    q21 = q0.getValue(11);
    q22 = q0.getValue(12);
    q23 = q0.getValue(13);
  }

  double rot2to1x, rot2to1y, rot2to1z;
  rot2to1(q10, q11, q12, q13, q20, q21, q22, q23,
          NULL, &rot2to1x, &rot2to1y, &rot2to1z);

  y.setValue(3, AscalA1(rot2to1x, rot2to1y, rot2to1z) - _initial_AscalA1);
  y.setValue(4, AscalA2(rot2to1x, rot2to1y, rot2to1z) - _initial_AscalA2);
}

/** Compute the vector of linear and angular positions of the free axes */
void PivotJointR::computehDoF(double time, BlockVector& q0, SiconosVector& y,
                              unsigned int axis)
{
  // Normally we fill y starting at axis up to the number of columns,
  // but in this case there is only one, so just don't do anything if
  // it doesn't match.
  if (axis != 0)
    return;

  SP::SiconosVector q1 = (q0.getAllVect())[0];
  double q10 = q1->getValue(3);
  double q11 = q1->getValue(4);
  double q12 = q1->getValue(5);
  double q13 = q1->getValue(6);

  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;

  if (q0.numberOfBlocks()>1)
  {
    SP::SiconosVector q2 = (q0.getAllVect())[1];
    q20 = q2->getValue(3);
    q21 = q2->getValue(4);
    q22 = q2->getValue(5);
    q23 = q2->getValue(6);
  }

  double rot2to1w, rot2to1x, rot2to1y, rot2to1z;
  rot2to1(q10, q11, q12, q13, q20, q21, q22, q23,
          &rot2to1w, &rot2to1x, &rot2to1y, &rot2to1z);

  // In case of joint constraints, it's okay to use dot product=0, but
  // in the case of the free axis we must measure the actual angle
  // using atan2 so that stops can be placed correctly.
  double Adot2to1 = AscalA(rot2to1x, rot2to1y, rot2to1z);
  double wrappedAngle = piwrap(2*atan2(rot2to1w, Adot2to1) - _initial_AscalA);

  // Count the number of twists around the angle, and report the
  // unwrapped angle.  Needed to implement joint stops near pi.
  if (wrappedAngle < -M_PI*3/4 && _previousAngle > M_PI*3/4)
    _twistCount ++;
  else if (wrappedAngle > M_PI*3/4 && _previousAngle < -M_PI*3/4)
    _twistCount --;
  _previousAngle = wrappedAngle;
  double unwrappedAngle = wrappedAngle + 2*M_PI*_twistCount;

  y.setValue(0, unwrappedAngle);
}

/** Compute the jacobian of linear and angular DoF with respect to some q */
void PivotJointR::computeJachqDoF(double time, Interaction& inter,
                                  SP::BlockVector q0, SimpleMatrix& jachq,
                                  unsigned int axis)
{
  // Normally we fill jachq starting at axis up to the number of rows,
  // but in this case there is only one, so just don't do anything if
  // it doesn't match.
  if (axis != 0)
    return;

  SP::SiconosVector q1 = (q0->getAllVect())[0];
  double q10 = q1->getValue(3);
  double q11 = q1->getValue(4);
  double q12 = q1->getValue(5);
  double q13 = q1->getValue(6);

  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;

  if (q0->numberOfBlocks()>1)
  {
    SP::SiconosVector q2 = (q0->getAllVect())[1];
    q20 = q2->getValue(3);
    q21 = q2->getValue(4);
    q22 = q2->getValue(5);
    q23 = q2->getValue(6);
  }

  jachq.setValue(0, 0, 0);
  jachq.setValue(0, 1, 0);
  jachq.setValue(0, 2, 0);

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

  double x0 = _cq2q103*q23;
  double x1 = _cq2q101*q21;
  double x2 = _cq2q102*q20;
  double x3 = _cq2q104*q22;
  double x4 = x0 - x1 - x2 - x3;
  double x5 = _cq2q104*q21;
  double x6 = _cq2q101*q22;
  double x7 = _cq2q102*q23;
  double x8 = _cq2q103*q20;
  double x9 = x5 - x6 - x7 - x8;
  double x10 = _cq2q102*q22;
  double x11 = _cq2q101*q23;
  double x12 = _cq2q103*q21;
  double x13 = _cq2q104*q20;
  double x14 = x10 - x11 - x12 - x13;
  double x15 = q11*x4;
  double x16 = q12*x9;
  double x17 = q13*x14;
  double x18 = _cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23;
  double x19 = q10*x18;
  double x20 = _A->getValue(0)*(q10*x4 + q11*x18 - q12*x14 + q13*x9) + _A->getValue(1)*(q10*x9 + q11*x14 + q12*x18 - q13*x4) + _A->getValue(2)*(q10*x14 - q11*x9 + q12*x4 + q13*x18);
  double x21 = 1.0/(pow(x20, 2) + pow(-x15 - x16 - x17 + x19, 2));
  double x22 = 2*x21*(x15 + x16 + x17 - x19);
  double x23 = 2*x20*x21;
  double x24 = -x5 + x6 + x7 + x8;
  double x25 = -x0 + x1 + x2 + x3;
  double x26 = -x10 + x11 + x12 + x13;
  double x27 = _cq2q101*q11;
  double x28 = _cq2q102*q10;
  double x29 = -x28;
  double x30 = x27 + x29;
  double x31 = _cq2q104*q12;
  double x32 = _cq2q103*q13;
  double x33 = -x32;
  double x34 = _cq2q101*q12;
  double x35 = _cq2q102*q13;
  double x36 = _cq2q103*q10;
  double x37 = -x36;
  double x38 = _cq2q104*q11;
  double x39 = -x38;
  double x40 = x37 + x39;
  double x41 = _cq2q101*q13;
  double x42 = _cq2q102*q12;
  double x43 = -x42;
  double x44 = x41 + x43;
  double x45 = _cq2q103*q11;
  double x46 = _cq2q104*q10;
  double x47 = -x46;
  double x48 = _cq2q101*q10;
  double x49 = _cq2q102*q11;
  double x50 = _cq2q103*q12;
  double x51 = _cq2q104*q13;
  double x52 = x50 + x51;
  double x53 = -x48;
  double x54 = -x45;
  double x55 = -x35;
  double x56 = -x31;
  double x57 = x49 + x53;
  double x58 = x33 + x56;
  double x59 = x47 + x54;
  double x60 = x34 + x55;

  jachq.setValue(0, 3, x18*x23 + x22*(_A->getValue(0)*x4 + _A->getValue(1)*x9 + _A->getValue(2)*x14));
  jachq.setValue(0, 4, x22*(_A->getValue(0)*x18 + _A->getValue(1)*x14 + _A->getValue(2)*x24) + x23*x25);
  jachq.setValue(0, 5, x22*(_A->getValue(0)*x26 + _A->getValue(1)*x18 + _A->getValue(2)*x4) + x23*x24);
  jachq.setValue(0, 6, x22*(_A->getValue(0)*x9 + _A->getValue(1)*x25 + _A->getValue(2)*x18) + x23*x26);

  if (q0->numberOfBlocks()<2)
    return;

  jachq.setValue(0, 7, 0);
  jachq.setValue(0, 8, 0);
  jachq.setValue(0, 9, 0);

  /*
   * sympy expression:
   *
   * for i in range(4): print('jachq.setValue(0, {}, {});'.format(i+10,e[i+4]))
  */

  jachq.setValue(0, 10, x22*(_A->getValue(0)*(x30 + x31 + x33) + _A->getValue(1)*(x34 + x35 + x40) + _A->getValue(2)*(x44 + x45 + x47)) + x23*(x48 + x49 + x52));
  jachq.setValue(0, 11, x22*(_A->getValue(0)*(-x49 + x52 + x53) + _A->getValue(1)*(x44 + x46 + x54) + _A->getValue(2)*(-x34 + x40 + x55)) + x23*(x30 + x32 + x56));
  jachq.setValue(0, 12, x22*(_A->getValue(0)*(-x41 + x43 + x59) + _A->getValue(1)*(-x50 + x51 + x57) + _A->getValue(2)*(x27 + x28 + x58)) + x23*(x37 + x38 + x60));
  jachq.setValue(0, 13, x22*(_A->getValue(0)*(x36 + x39 + x60) + _A->getValue(1)*(-x27 + x29 + x58) + _A->getValue(2)*(x50 - x51 + x57)) + x23*(x41 + x42 + x59));
}


/** Return the normal of the angular DoF axis of rotation.
 * \param axis must be 0 */
void PivotJointR::_normalDoF(const BlockVector& q0, SiconosVector& ans, int axis,
                             bool absoluteRef)
{
  assert(axis == 0);
  if (axis != 0) return;

  // We assume that A is normalized.
  ans = *_A;

  if (absoluteRef)
  {
    SP::SiconosVector tmp2(std11::make_shared<SiconosVector>(ans));
    changeFrameAbsToBody(q0.getAllVect()[0], tmp2);
    ans = *tmp2;
    return;
    ::boost::math::quaternion<double> q1(q0.getValue(3), q0.getValue(4),
                                         q0.getValue(5), q0.getValue(6));

    // _A is in the q1 frame, so change it to the inertial frame.
    ::boost::math::quaternion<double> aq(0, (*_A)(0), (*_A)(1), (*_A)(2));
    //::boost::math::quaternion<double> tmp( (1.0/q1) * aq * q1 );
    ::boost::math::quaternion<double> tmp( q1 * aq / q1 );
    ans(0) = tmp.R_component_2();
    ans(1) = tmp.R_component_3();
    ans(2) = tmp.R_component_4();
    printf("PivotJointR::_normalDoF: _A = (%0.02f, %0.02f, %0.02f)\n",
           (*_A)(0), (*_A)(1), (*_A)(2));
    printf("PivotJointR::_normalDoF: ans = (%0.02f, %0.02f, %0.02f)\n",
           ans(0), ans(1), ans(2));
  }
}
