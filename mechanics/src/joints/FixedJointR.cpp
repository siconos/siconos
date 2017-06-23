/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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

/*! \file FixedJointR.cpp */

#include "FixedJointR.hpp"
#include <NewtonEulerDS.hpp>
#include <Interaction.hpp>
#include <boost/math/quaternion.hpp>
#include <BlockVector.hpp>
#include <cfloat>
#include <iostream>

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

/*
 * This file contains some code generated using sympy.  The following
 * is the necessary predule:
 *
 * from sympy import Symbol
 * from sympy import ccode
 * import numpy as np
 *
 * G10G20d1 = np.array([0,Symbol('_G10G20d1x'),Symbol('_G10G20d1y'),Symbol('_G10G20d1z')])
 * cq2q10 = np.array([Symbol('_cq2q101'),Symbol('_cq2q102'),
 *                    Symbol('_cq2q103'),Symbol('_cq2q104')])
 * q1 = np.array([Symbol('q10'), Symbol('q11'), Symbol('q12'), Symbol('q13')])
 * q2 = np.array([Symbol('q20'), Symbol('q21'), Symbol('q22'), Symbol('q23')])
 * G1 = np.array([0, Symbol('X1'), Symbol('Y1'), Symbol('Z1')])
 * G2 = np.array([0, Symbol('X2'), Symbol('Y2'), Symbol('Z2')])
 *
 * qinv = lambda q: np.array([q[0],-q[1],-q[2],-q[3]])
 * qmul = lambda a,b: np.array([
 *          a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
 *          a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
 *          a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
 *          a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]])
 *
 * unrot = lambda V,q: qmul(qinv(q), qmul(V, q))
 * rot = lambda V,q: qmul(q, qmul(V, qinv(q)))
 */

FixedJointR::FixedJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2)
  : NewtonEulerJointR()
{
  setInitialConditions(d1->q(), d2 ? d2->q() : SP::SiconosVector());
}

void FixedJointR::setInitialConditions(SP::SiconosVector q1, SP::SiconosVector q2)
{
  ::boost::math::quaternion<double> quat1((*q1)(3), (*q1)(4), (*q1)(5), (*q1)(6));
  ::boost::math::quaternion<double> quat2(q2 ? (*q2)(3) : 1, q2 ? (*q2)(4) : 0,
                                          q2 ? (*q2)(5) : 0, q2 ? (*q2)(6) : 0);
  ::boost::math::quaternion<double> cq2q10(1.0 / quat2 * quat1);

  _cq2q101 = cq2q10.R_component_1();
  _cq2q102 = cq2q10.R_component_2();
  _cq2q103 = cq2q10.R_component_3();
  _cq2q104 = cq2q10.R_component_4();

  ::boost::math::quaternion<double>    quatG10G20_abs(
    0, (q2 ? q2->getValue(0) : 0) - q1->getValue(0),
    (q2 ? q2->getValue(1) : 0) - q1->getValue(1),
    (q2 ? q2->getValue(2) : 0) - q1->getValue(2));
  ::boost::math::quaternion<double>    quatBuff(0, 0, 0, 0);

  quatBuff = 1.0/quat1 * quatG10G20_abs * quat1;
  _G10G20d1x = quatBuff.R_component_2();
  _G10G20d1y = quatBuff.R_component_3();
  _G10G20d1z = quatBuff.R_component_4();
}

void FixedJointR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  double X1 = q0.getValue(0);
  double Y1 = q0.getValue(1);
  double Z1 = q0.getValue(2);
  double q10 = q0.getValue(3);
  double q11 = q0.getValue(4);
  double q12 = q0.getValue(5);
  double q13 = q0.getValue(6);

  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;

  if (q0.numberOfBlocks()>1)
  {
    X2 = q0.getValue(7);
    Y2 = q0.getValue(8);
    Z2 = q0.getValue(9);
    q20 = q0.getValue(10);
    q21 = q0.getValue(11);
    q22 = q0.getValue(12);
    q23 = q0.getValue(13);
  }

  /* sympy expression:
   *
   * G1G2d1 = unrot(G2-G1, q1) - G10G20d1
   * q2to1 = qmul(q1,qinv(qmul(q2,cq2q10)))
   *
   * H = list(G1G2d1[1:]) + list(q2to1[1:])
   */

  y.setValue(0, -_G10G20d1x + q10*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2))
             - q11*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2))
             - q12*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
             + q13*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2)));
  y.setValue(1, -_G10G20d1y + q10*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2))
             + q11*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
             - q12*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2))
             - q13*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2)));
  y.setValue(2, -_G10G20d1z + q10*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
             - q11*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2))
             + q12*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2))
             - q13*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2)));
  y.setValue(3, q10*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
             + q11*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
             + q12*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
             - q13*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21));
  y.setValue(4, q10*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
             - q11*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
             + q12*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
             + q13*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22));
  y.setValue(5, q10*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
             + q11*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
             - q12*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
             + q13*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23));
}

void FixedJointR::computeJachq(double time, Interaction& inter, SP::BlockVector q0)
{
  _jachq->zero();

  SP::SiconosVector q1 = (q0->getAllVect())[0];
  double X1 = q1->getValue(0);
  double Y1 = q1->getValue(1);
  double Z1 = q1->getValue(2);
  double q10 = q1->getValue(3);
  double q11 = q1->getValue(4);
  double q12 = q1->getValue(5);
  double q13 = q1->getValue(6);

  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;

  if(q0->numberOfBlocks()>1)
  {
    SP::SiconosVector q2 = (q0->getAllVect())[1];
    X2 = q2->getValue(0);
    Y2 = q2->getValue(1);
    Z2 = q2->getValue(2);
    q20 = q2->getValue(3);
    q21 = q2->getValue(4);
    q22 = q2->getValue(5);
    q23 = q2->getValue(6);
    Jd1d2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);
  }
  else
    Jd1(X1, Y1, Z1, q10, q11, q12, q13);
}

void FixedJointR::Jd1d2(double X1, double Y1, double Z1,
                        double q10, double q11, double q12, double q13,
                        double X2, double Y2, double Z2,
                        double q20, double q21, double q22, double q23)
{
  /* sympy expression:
   *
   * H = list(G1G2d1[1:]) + list(q2to1[1:])
   * dq = list(G1[1:])+list(q1)+list(G2[1:])+list(q2)
   * jachq = [[h.diff(d) for d in dq] for h in H]
   */

  _jachq->setValue(0, 0, -pow(q10, 2) - pow(q11, 2) + pow(q12, 2) + pow(q13, 2));
  _jachq->setValue(0, 1, -2*q10*q13 - 2*q11*q12);
  _jachq->setValue(0, 2, 2*q10*q12 - 2*q11*q13);
  _jachq->setValue(0, 3, 2*q10*(-X1 + X2) - 2*q12*(-Z1 + Z2) + 2*q13*(-Y1 + Y2));
  _jachq->setValue(0, 4, q11*(-X1 + X2) - q11*(X1 - X2) + q12*(-Y1 + Y2)
                   - q12*(Y1 - Y2) + 2*q13*(-Z1 + Z2));
  _jachq->setValue(0, 5, -q10*(-Z1 + Z2) + q10*(Z1 - Z2) + q11*(-Y1 + Y2)
                   - q11*(Y1 - Y2) - 2*q12*(-X1 + X2));
  _jachq->setValue(0, 6, 2*q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q11*(Z1 - Z2)
                   - q13*(-X1 + X2) + q13*(X1 - X2));
  _jachq->setValue(0, 7, pow(q10, 2) + pow(q11, 2) - pow(q12, 2) - pow(q13, 2));
  _jachq->setValue(0, 8, 2*q10*q13 + 2*q11*q12);
  _jachq->setValue(0, 9, -2*q10*q12 + 2*q11*q13);
  _jachq->setValue(0, 10, 0);
  _jachq->setValue(0, 11, 0);
  _jachq->setValue(0, 12, 0);
  _jachq->setValue(0, 13, 0);
  _jachq->setValue(1, 0, 2*q10*q13 - 2*q11*q12);
  _jachq->setValue(1, 1, -pow(q10, 2) + pow(q11, 2) - pow(q12, 2) + pow(q13, 2));
  _jachq->setValue(1, 2, -2*q10*q11 - 2*q12*q13);
  _jachq->setValue(1, 3, 2*q10*(-Y1 + Y2) + 2*q11*(-Z1 + Z2) - 2*q13*(-X1 + X2));
  _jachq->setValue(1, 4, 2*q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q11*(Y1 - Y2)
                   + q12*(-X1 + X2) - q12*(X1 - X2));
  _jachq->setValue(1, 5, 2*q11*(-X1 + X2) + q12*(-Y1 + Y2) - q12*(Y1 - Y2)
                   + q13*(-Z1 + Z2) - q13*(Z1 - Z2));
  _jachq->setValue(1, 6, -q10*(-X1 + X2) + q10*(X1 - X2) + q12*(-Z1 + Z2)
                   - q12*(Z1 - Z2) - 2*q13*(-Y1 + Y2));
  _jachq->setValue(1, 7, -2*q10*q13 + 2*q11*q12);
  _jachq->setValue(1, 8, pow(q10, 2) - pow(q11, 2) + pow(q12, 2) - pow(q13, 2));
  _jachq->setValue(1, 9, 2*q10*q11 + 2*q12*q13);
  _jachq->setValue(1, 10, 0);
  _jachq->setValue(1, 11, 0);
  _jachq->setValue(1, 12, 0);
  _jachq->setValue(1, 13, 0);
  _jachq->setValue(2, 0, -2*q10*q12 - 2*q11*q13);
  _jachq->setValue(2, 1, 2*q10*q11 - 2*q12*q13);
  _jachq->setValue(2, 2, -pow(q10, 2) + pow(q11, 2) + pow(q12, 2) - pow(q13, 2));
  _jachq->setValue(2, 3, 2*q10*(-Z1 + Z2) - 2*q11*(-Y1 + Y2) + 2*q12*(-X1 + X2));
  _jachq->setValue(2, 4, -q10*(-Y1 + Y2) + q10*(Y1 - Y2) - 2*q11*(-Z1 + Z2)
                   + q13*(-X1 + X2) - q13*(X1 - X2));
  _jachq->setValue(2, 5, 2*q10*(-X1 + X2) - q12*(-Z1 + Z2) + q12*(Z1 - Z2)
                   + q13*(-Y1 + Y2) - q13*(Y1 - Y2));
  _jachq->setValue(2, 6, q11*(-X1 + X2) - q11*(X1 - X2) + 2*q12*(-Y1 + Y2)
                   + q13*(-Z1 + Z2) - q13*(Z1 - Z2));
  _jachq->setValue(2, 7, 2*q10*q12 + 2*q11*q13);
  _jachq->setValue(2, 8, -2*q10*q11 + 2*q12*q13);
  _jachq->setValue(2, 9, pow(q10, 2) - pow(q11, 2) - pow(q12, 2) + pow(q13, 2));
  _jachq->setValue(2, 10, 0);
  _jachq->setValue(2, 11, 0);
  _jachq->setValue(2, 12, 0);
  _jachq->setValue(2, 13, 0);
  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);
  _jachq->setValue(3, 3, -_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22);
  _jachq->setValue(3, 4, _cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23);
  _jachq->setValue(3, 5, -_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20);
  _jachq->setValue(3, 6, _cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21);
  _jachq->setValue(3, 7, 0);
  _jachq->setValue(3, 8, 0);
  _jachq->setValue(3, 9, 0);
  _jachq->setValue(3, 10, _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12);
  _jachq->setValue(3, 11, -_cq2q101*q10 - _cq2q102*q11 - _cq2q103*q12 - _cq2q104*q13);
  _jachq->setValue(3, 12, _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10);
  _jachq->setValue(3, 13, -_cq2q101*q12 + _cq2q102*q13 + _cq2q103*q10 - _cq2q104*q11);
  _jachq->setValue(4, 0, 0);
  _jachq->setValue(4, 1, 0);
  _jachq->setValue(4, 2, 0);
  _jachq->setValue(4, 3, -_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21);
  _jachq->setValue(4, 4, _cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20);
  _jachq->setValue(4, 5, _cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23);
  _jachq->setValue(4, 6, -_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22);
  _jachq->setValue(4, 7, 0);
  _jachq->setValue(4, 8, 0);
  _jachq->setValue(4, 9, 0);
  _jachq->setValue(4, 10, _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11);
  _jachq->setValue(4, 11, -_cq2q101*q13 - _cq2q102*q12 + _cq2q103*q11 + _cq2q104*q10);
  _jachq->setValue(4, 12, -_cq2q101*q10 - _cq2q102*q11 - _cq2q103*q12 - _cq2q104*q13);
  _jachq->setValue(4, 13, _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12);
  _jachq->setValue(5, 0, 0);
  _jachq->setValue(5, 1, 0);
  _jachq->setValue(5, 2, 0);
  _jachq->setValue(5, 3, -_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20);
  _jachq->setValue(5, 4, -_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21);
  _jachq->setValue(5, 5, _cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22);
  _jachq->setValue(5, 6, _cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23);
  _jachq->setValue(5, 7, 0);
  _jachq->setValue(5, 8, 0);
  _jachq->setValue(5, 9, 0);
  _jachq->setValue(5, 10, _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10);
  _jachq->setValue(5, 11, _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11);
  _jachq->setValue(5, 12, -_cq2q101*q11 + _cq2q102*q10 - _cq2q103*q13 + _cq2q104*q12);
  _jachq->setValue(5, 13, -_cq2q101*q10 - _cq2q102*q11 - _cq2q103*q12 - _cq2q104*q13);
}

void FixedJointR::Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13)
{
  /* sympy expression:
   *
   * q2 = np.array([1,0,0,0])
   * G2 = np.array([0, 0,0,0])
   * G1G2d1 = unrot(G2-G1, q1) - G10G20d1
   * q2to1 = qmul(q1,qinv(qmul(q2,cq2q10)))
   *
   * H = list(G1G2d1[1:]) + list(q2to1[1:])
   * dq = list(G1[1:])+list(q1)
   * jachq = [[h.diff(d) for d in dq] for h in H]
   */

  _jachq->setValue(0, 0, -pow(q10, 2) - pow(q11, 2) + pow(q12, 2) + pow(q13, 2));
  _jachq->setValue(0, 1, -2*q10*q13 - 2*q11*q12);
  _jachq->setValue(0, 2, 2*q10*q12 - 2*q11*q13);
  _jachq->setValue(0, 3, -2*X1*q10 - 2*Y1*q13 + 2*Z1*q12);
  _jachq->setValue(0, 4, -2*X1*q11 - 2*Y1*q12 - 2*Z1*q13);
  _jachq->setValue(0, 5, 2*X1*q12 - 2*Y1*q11 + 2*Z1*q10);
  _jachq->setValue(0, 6, 2*X1*q13 - 2*Y1*q10 - 2*Z1*q11);
  _jachq->setValue(1, 0, 2*q10*q13 - 2*q11*q12);
  _jachq->setValue(1, 1, -pow(q10, 2) + pow(q11, 2) - pow(q12, 2) + pow(q13, 2));
  _jachq->setValue(1, 2, -2*q10*q11 - 2*q12*q13);
  _jachq->setValue(1, 3, 2*X1*q13 - 2*Y1*q10 - 2*Z1*q11);
  _jachq->setValue(1, 4, -2*X1*q12 + 2*Y1*q11 - 2*Z1*q10);
  _jachq->setValue(1, 5, -2*X1*q11 - 2*Y1*q12 - 2*Z1*q13);
  _jachq->setValue(1, 6, 2*X1*q10 + 2*Y1*q13 - 2*Z1*q12);
  _jachq->setValue(2, 0, -2*q10*q12 - 2*q11*q13);
  _jachq->setValue(2, 1, 2*q10*q11 - 2*q12*q13);
  _jachq->setValue(2, 2, -pow(q10, 2) + pow(q11, 2) + pow(q12, 2) - pow(q13, 2));
  _jachq->setValue(2, 3, -2*X1*q12 + 2*Y1*q11 - 2*Z1*q10);
  _jachq->setValue(2, 4, -2*X1*q13 + 2*Y1*q10 + 2*Z1*q11);
  _jachq->setValue(2, 5, -2*X1*q10 - 2*Y1*q13 + 2*Z1*q12);
  _jachq->setValue(2, 6, -2*X1*q11 - 2*Y1*q12 - 2*Z1*q13);
  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);
  _jachq->setValue(3, 3, -_cq2q102);
  _jachq->setValue(3, 4, _cq2q101);
  _jachq->setValue(3, 5, -_cq2q104);
  _jachq->setValue(3, 6, _cq2q103);
  _jachq->setValue(4, 0, 0);
  _jachq->setValue(4, 1, 0);
  _jachq->setValue(4, 2, 0);
  _jachq->setValue(4, 3, -_cq2q103);
  _jachq->setValue(4, 4, _cq2q104);
  _jachq->setValue(4, 5, _cq2q101);
  _jachq->setValue(4, 6, -_cq2q102);
  _jachq->setValue(5, 0, 0);
  _jachq->setValue(5, 1, 0);
  _jachq->setValue(5, 2, 0);
  _jachq->setValue(5, 3, -_cq2q104);
  _jachq->setValue(5, 4, -_cq2q103);
  _jachq->setValue(5, 5, _cq2q102);
  _jachq->setValue(5, 6, _cq2q101);
}
