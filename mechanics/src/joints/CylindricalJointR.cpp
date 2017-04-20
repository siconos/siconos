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

/*! \file CylindricalJointR.cpp
*/

#include "CylindricalJointR.hpp"
#include <NewtonEulerDS.hpp>
#include <Interaction.hpp>
#include <boost/math/quaternion.hpp>
#include <BlockVector.hpp>

#include <iostream>

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

/*
 * This file contains some code generated using sympy.  The following
 * is the necessary prelude:
 *
 * from sympy import Symbol
 * import numpy as np
 *
 * q1 = np.array([Symbol('q10'), Symbol('q11'), Symbol('q12'), Symbol('q13')])
 * q2 = np.array([Symbol('q20'), Symbol('q21'), Symbol('q22'), Symbol('q23')])
 * cq2q10 = np.array([Symbol('_cq2q101'),Symbol('_cq2q102'),
 *                    Symbol('_cq2q103'),Symbol('_cq2q104')])
 * G10G20d1 = np.array([0, Symbol('_G10G20d1x'),
 *                      Symbol('_G10G20d1y'), Symbol('_G10G20d1z')])
 * G1 = np.array([0, Symbol('X1'), Symbol('Y1'), Symbol('Z1')])
 * G2 = np.array([0, Symbol('X2'), Symbol('Y2'), Symbol('Z2')])
 * V1 = np.array([0, Symbol('_V1x'), Symbol('_V1y'), Symbol('_V1z')])
 * V2 = np.array([0, Symbol('_V2x'), Symbol('_V2y'), Symbol('_V2z')])
 * A = np.array([0, Symbol('_axis0->getValue(0)'),
 *               Symbol('_axis0->getValue(1)'), Symbol('_axis0->getValue(2)')])
 * P = np.array([0, Symbol('_P->getValue(0)'),
 *               Symbol('_P->getValue(1)'), Symbol('_P->getValue(2)')])
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

/**axe is the axis of the cylindrical joint, in the frame of the first DS, d1.*/
CylindricalJointR::CylindricalJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2,
                                     SP::SiconosVector A)
  : NewtonEulerR()
  , _P(std11::make_shared<SiconosVector>(3))
{
  // in the two-DS case, _P is unused.
  _P->zero();

  _axis0 = A;
  computeFromInitialPosition(d1->q(),d2->q());
}

/*axis is the axis of the cylindrical joint, in the absolute frame.*/
CylindricalJointR::CylindricalJointR(SP::NewtonEulerDS d1,
                                     SP::SiconosVector P, SP::SiconosVector A,
                                     bool absoluteRef)
  : NewtonEulerR()
  , _P(std11::make_shared<SiconosVector>(3))
{
  if (P) *_P = *P; else _P->zero();
  if (!absoluteRef) {
    // add d1's frame
    _P->setValue(0, _P->getValue(0) + d1->q()->getValue(0));
    _P->setValue(1, _P->getValue(1) + d1->q()->getValue(1));
    _P->setValue(2, _P->getValue(2) + d1->q()->getValue(2));
  }

  _axis0 = A;
  computeFromInitialPosition(d1->q());
}

void CylindricalJointR::computeFromInitialPosition(SP::SiconosVector q1, SP::SiconosVector q2)
{
  computeV1V2FromAxis();
  SP::SiconosVector q2i(new SiconosVector(7));
  q2i->zero();
  q2i->setValue(0, _P->getValue(0));
  q2i->setValue(1, _P->getValue(1));
  q2i->setValue(2, _P->getValue(2));
  q2i->setValue(3, 1);

  if(q2)
    *q2i = *q2;

  ::boost::math::quaternion<double>    quat1(q1->getValue(3), q1->getValue(4), q1->getValue(5), q1->getValue(6));
  ::boost::math::quaternion<double>    quat2(q2i->getValue(3), q2i->getValue(4), q2i->getValue(5), q2i->getValue(6));
  ::boost::math::quaternion<double>    quat1_inv(q1->getValue(3), -q1->getValue(4), -q1->getValue(5), -q1->getValue(6));
  ::boost::math::quaternion<double>    quatG10G20_abs(
    0, q2i->getValue(0) - q1->getValue(0),
    q2i->getValue(1) - q1->getValue(1),
    q2i->getValue(2) - q1->getValue(2));
  ::boost::math::quaternion<double>    quatBuff(0, 0, 0, 0);

  quatBuff = quat1_inv * (quatG10G20_abs * quat1);
  _G10G20d1x = quatBuff.R_component_2();
  _G10G20d1y = quatBuff.R_component_3();
  _G10G20d1z = quatBuff.R_component_4();
  std::cout << "G10G20d1: " << quatBuff << std::endl;

  quatBuff = 1.0/quat2 * quat1;
  _cq2q101 = quatBuff.R_component_1();
  _cq2q102 = quatBuff.R_component_2();
  _cq2q103 = quatBuff.R_component_3();
  _cq2q104 = quatBuff.R_component_4();
}

void CylindricalJointR::computeV1V2FromAxis()
{
  _V1.reset(new SiconosVector(3));
  _V2.reset(new SiconosVector(3));
  _V1->zero();
  _V2->zero();
  //build _V1
  if(_axis0->getValue(0) > _axis0->getValue(1))
    if(_axis0->getValue(0) > _axis0->getValue(2))
    {
      _V1->setValue(1, -_axis0->getValue(0));
      _V1->setValue(0, _axis0->getValue(1));
    }
    else
    {
      _V1->setValue(1, -_axis0->getValue(2));
      _V1->setValue(2, _axis0->getValue(1));
    }
  else if(_axis0->getValue(2) > _axis0->getValue(1))
  {
    _V1->setValue(1, -_axis0->getValue(2));
    _V1->setValue(2, _axis0->getValue(1));
  }
  else
  {
    _V1->setValue(1, -_axis0->getValue(0));
    _V1->setValue(0, _axis0->getValue(1));
  }
  double aux = 1 / _V1->norm2();
  scal(aux, *_V1, *_V1);
  cross_product(*_axis0, *_V1, *_V2);
  _V1x = _V1->getValue(0);
  _V1y = _V1->getValue(1);
  _V1z = _V1->getValue(2);
  _V2x = _V2->getValue(0);
  _V2y = _V2->getValue(1);
  _V2z = _V2->getValue(2);
}

void CylindricalJointR::computeJachq(double time, Interaction& inter,  SP::BlockVector q0)
{
  SP::SiconosVector q1 = (q0->getAllVect())[0];
  double X1 = q1->getValue(0);
  double Y1 = q1->getValue(1);
  double Z1 = q1->getValue(2);
  double q10 = q1->getValue(3);
  double q11 = q1->getValue(4);
  double q12 = q1->getValue(5);
  double q13 = q1->getValue(6);

  if(q0->getNumberOfBlocks()>1)
  {
    SP::SiconosVector q2 = (q0->getAllVect())[1];
    double X2 = q2->getValue(0);
    double Y2 = q2->getValue(1);
    double Z2 = q2->getValue(2);
    double q20 = q2->getValue(3);
    double q21 = q2->getValue(4);
    double q22 = q2->getValue(5);
    double q23 = q2->getValue(6);
    Jd1d2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);
  }
  else
    Jd1(X1, Y1, Z1, q10, q11, q12, q13);
}

void CylindricalJointR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  SP::SiconosVector q1 = (q0.getAllVect())[0];
  double X1 = q1->getValue(0);
  double Y1 = q1->getValue(1);
  double Z1 = q1->getValue(2);
  double q10 = q1->getValue(3);
  double q11 = q1->getValue(4);
  double q12 = q1->getValue(5);
  double q13 = q1->getValue(6);
  double X2 = _P->getValue(0);
  double Y2 = _P->getValue(1);
  double Z2 = _P->getValue(2);
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;

  if (q0.getNumberOfBlocks()>1)
  {
    SP::SiconosVector q2 = (q0.getAllVect())[1];
    X2 = q2->getValue(0);
    Y2 = q2->getValue(1);
    Z2 = q2->getValue(2);
    q20 = q2->getValue(3);
    q21 = q2->getValue(4);
    q22 = q2->getValue(5);
    q23 = q2->getValue(6);
  }

  y.setValue(0, H1(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(1, H2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(2, H3(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(3, H4(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
}

// The rest of the code is generated.
// we can disable some warning
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

/* sympy expression:
 *
 * G1G2d1 = unrot(G2-G1, q1)
 *
 * H  = lambda V: np.dot(V, G1G2d1)
 * H0 = lambda V: np.dot(V, G10G20d1)
 *
 * H1 = H(V1) - H0(V1)
 * H2 = H(V2) - H0(V2)
 *
 */

double CylindricalJointR::H1(
  double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
  double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  return -_G10G20d1x*_V1x - _G10G20d1y*_V1y - _G10G20d1z*_V1z
    + _V1x*(q10*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2))
            - q11*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2))
            - q12*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
            + q13*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2)))
    + _V1y*(q10*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2))
            + q11*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
            - q12*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2))
            - q13*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2)))
    + _V1z*(q10*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
            - q11*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2))
            + q12*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2))
            - q13*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2)));
}

double CylindricalJointR::H2(
  double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
  double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  return -_G10G20d1x*_V2x - _G10G20d1y*_V2y - _G10G20d1z*_V2z
    + _V2x*(q10*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2))
            - q11*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2))
            - q12*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
            + q13*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2)))
    + _V2y*(q10*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2))
            + q11*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
            - q12*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2))
            - q13*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2)))
    + _V2z*(q10*(q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q12*(-X1 + X2))
            - q11*(q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q13*(-X1 + X2))
            + q12*(q10*(-X1 + X2) - q12*(-Z1 + Z2) + q13*(-Y1 + Y2))
            - q13*(-q11*(-X1 + X2) - q12*(-Y1 + Y2) - q13*(-Z1 + Z2)));
}

/* sympy expression:
 *
 * q2to1 = qmul(q1,qinv(qmul(q2,cq2q10)))
 *
 * H3 = np.dot(rot(V1,q1), q2to1)
 * H4 = np.dot(rot(V2,q1), q2to1)
 */

double CylindricalJointR::H3(
  double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
  double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  return
    (q10*(_V1x*q10 - _V1y*q13 + _V1z*q12) + q11*(_V1x*q11 + _V1y*q12 + _V1z*q13)
     + q12*(-_V1x*q12 + _V1y*q11 + _V1z*q10) - q13*(_V1x*q13 + _V1y*q10 - _V1z*q11))
    *(q10*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
      + q11*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
      + q12*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
      - q13*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21))
    + (q10*(_V1x*q11 + _V1y*q12 + _V1z*q13) - q11*(_V1x*q10 - _V1y*q13 + _V1z*q12)
       - q12*(_V1x*q13 + _V1y*q10 - _V1z*q11) - q13*(-_V1x*q12 + _V1y*q11 + _V1z*q10))
    *(q10*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
      - q11*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
      - q12*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
      - q13*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20))
    + (q10*(-_V1x*q12 + _V1y*q11 + _V1z*q10) + q11*(_V1x*q13 + _V1y*q10 - _V1z*q11)
       - q12*(_V1x*q10 - _V1y*q13 + _V1z*q12) + q13*(_V1x*q11 + _V1y*q12 + _V1z*q13))
    *(q10*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
      + q11*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
      - q12*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
      + q13*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23))
    + (q10*(_V1x*q13 + _V1y*q10 - _V1z*q11) - q11*(-_V1x*q12 + _V1y*q11 + _V1z*q10)
       + q12*(_V1x*q11 + _V1y*q12 + _V1z*q13) + q13*(_V1x*q10 - _V1y*q13 + _V1z*q12))
    *(q10*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
      - q11*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
      + q12*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
      + q13*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22));
}

double CylindricalJointR::H4(
  double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
  double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  return
    (q10*(_V2x*q10 - _V2y*q13 + _V2z*q12) + q11*(_V2x*q11 + _V2y*q12 + _V2z*q13)
     + q12*(-_V2x*q12 + _V2y*q11 + _V2z*q10) - q13*(_V2x*q13 + _V2y*q10 - _V2z*q11))
    *(q10*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
      + q11*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
      + q12*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
      - q13*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21))
    + (q10*(_V2x*q11 + _V2y*q12 + _V2z*q13) - q11*(_V2x*q10 - _V2y*q13 + _V2z*q12)
       - q12*(_V2x*q13 + _V2y*q10 - _V2z*q11) - q13*(-_V2x*q12 + _V2y*q11 + _V2z*q10))
    *(q10*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
      - q11*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
      - q12*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
      - q13*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20))
    + (q10*(-_V2x*q12 + _V2y*q11 + _V2z*q10) + q11*(_V2x*q13 + _V2y*q10 - _V2z*q11)
       - q12*(_V2x*q10 - _V2y*q13 + _V2z*q12) + q13*(_V2x*q11 + _V2y*q12 + _V2z*q13))
    *(q10*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
      + q11*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
      - q12*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
      + q13*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23))
    + (q10*(_V2x*q13 + _V2y*q10 - _V2z*q11) - q11*(-_V2x*q12 + _V2y*q11 + _V2z*q10)
       + q12*(_V2x*q11 + _V2y*q12 + _V2z*q13) + q13*(_V2x*q10 - _V2y*q13 + _V2z*q12))
    *(q10*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
      - q11*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
      + q12*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
      + q13*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22));
}

void CylindricalJointR::Jd1d2(
  double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
  double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  /*
   * sympy expression:
   *
   * H = [H1,H2,H3,H4]
   * dq = list(G1[1:])+list(q1)+list(G2[1:])+list(q2)
   *
   * jachq = [[h.diff(d) for d in dq] for h in H]
   *
   */

  /* Prismatic constraints (H1, H2)
   */
  const double x0 = 2*q10*q12;
  const double x1 = 2*q11*q13 + x0;
  const double x2 = _V1z*x1;
  const double x3 = 2*q10*q13;
  const double x4 = -2*q11*q12 + x3;
  const double x5 = _V1y*x4;
  const double x6 = pow(q10, 2);
  const double x7 = pow(q11, 2);
  const double x8 = pow(q12, 2);
  const double x9 = -x8;
  const double x10 = pow(q13, 2);
  const double x11 = -x10;
  const double x12 = x11 + x6 + x7 + x9;
  const double x13 = _V1x*x12;
  const double x14 = 2*q11*q12 + x3;
  const double x15 = _V1x*x14;
  const double x16 = 2*q10*q11;
  const double x17 = -2*q12*q13 + x16;
  const double x18 = _V1z*x17;
  const double x19 = x6 - x7;
  const double x20 = x11 + x19 + x8;
  const double x21 = _V1y*x20;
  const double x22 = 2*q12*q13 + x16;
  const double x23 = _V1y*x22;
  const double x24 = -2*q11*q13 + x0;
  const double x25 = _V1x*x24;
  const double x26 = x10 + x19 + x9;
  const double x27 = _V1z*x26;
  const double x28 = X1 - X2;
  const double x29 = q10*x28;
  const double x30 = Y1 - Y2;
  const double x31 = Z1 - Z2;
  const double x32 = q12*x31;
  const double x33 = q13*x30 + x29 - x32;
  const double x34 = q10*x30;
  const double x35 = q13*x28;
  const double x36 = q11*x31 + x34 - x35;
  const double x37 = q10*x31;
  const double x38 = q11*x30;
  const double x39 = q12*x28 + x37 - x38;
  const double x40 = 2*q11*x28 + 2*q12*x30 + 2*q13*x31;
  const double x41 = -2*q10*x31 + 2*q11*x30 - 2*q12*x28;
  const double x42 = -Y1 + Y2;
  const double x43 = -Z1 + Z2;
  const double x44 = -X1 + X2;
  const double x45 = q10*x42 + 2*q11*x43 - q13*x44 - x34 + x35;
  const double x46 = -2*q10*x28 + 2*q12*x31 - 2*q13*x30;
  const double x47 = q10*x43 - q11*x42 + 2*q12*x44 - x37 + x38;
  const double x48 = -2*q10*x30 - 2*q11*x31 + 2*q13*x28;
  const double x49 = q10*x44 - q12*x43 + 2*q13*x42 - x29 + x32;
  const double x50 = _V2z*x1;
  const double x51 = _V2y*x4;
  const double x52 = _V2x*x12;
  const double x53 = _V2x*x14;
  const double x54 = _V2z*x17;
  const double x55 = _V2y*x20;
  const double x56 = _V2y*x22;
  const double x57 = _V2x*x24;
  const double x58 = _V2z*x26;
  const double x59 = -_cq2q101*q20 + _cq2q102*q21 + _cq2q103*q22 + _cq2q104*q23;
  const double x60 = _V1x*q11 + _V1y*q12 + _V1z*q13;
  const double x61 = _V1x*q10 - _V1y*q13 + _V1z*q12;
  const double x62 = _V1x*q13 + _V1y*q10 - _V1z*q11;
  const double x63 = -_V1x*q12 + _V1y*q11 + _V1z*q10;
  const double x64 = -q10*x60 + q11*x61 + q12*x62 + q13*x63;
  const double x65 = _cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22;
  const double x66 = q10*x61 + q11*x60 + q12*x63 - q13*x62;
  const double x67 = _cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21;
  const double x68 = q10*x62 - q11*x63 + q12*x60 + q13*x61;
  const double x69 = _cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20;
  const double x70 = q10*x63 + q11*x62 - q12*x61 + q13*x60;
  const double x71 = 2*q10*x65 + 2*q11*x59 + 2*q12*x69 - 2*q13*x67;
  const double x72 = 2*q10*x69 + 2*q11*x67 - 2*q12*x65 + 2*q13*x59;
  const double x73 = 2*q10*x67 - 2*q11*x69 + 2*q12*x59 + 2*q13*x65;
  const double x74 = 2*_V1x*q11 + 2*_V1y*q12 + 2*_V1z*q13;
  const double x75 = -q10*x65 - q11*x59 - q12*x69 + q13*x67;
  const double x76 = -2*_V1x*q12 + 2*_V1y*q11 + 2*_V1z*q10;
  const double x77 = -q10*x67 + q11*x69 - q12*x59 - q13*x65;
  const double x78 = 2*_V1x*q13 + 2*_V1y*q10 - 2*_V1z*q11;
  const double x79 = -q10*x69 - q11*x67 + q12*x65 - q13*x59;
  const double x80 = 2*_V1x*q10 - 2*_V1y*q13 + 2*_V1z*q12;
  const double x81 = _cq2q101*q10 + _cq2q102*q11 + _cq2q103*q12 + _cq2q104*q13;
  const double x82 = _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12;
  const double x83 = _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11;
  const double x84 = _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10;
  const double x85 = _V2x*q11 + _V2y*q12 + _V2z*q13;
  const double x86 = _V2x*q10 - _V2y*q13 + _V2z*q12;
  const double x87 = _V2x*q13 + _V2y*q10 - _V2z*q11;
  const double x88 = -_V2x*q12 + _V2y*q11 + _V2z*q10;
  const double x89 = -q10*x85 + q11*x86 + q12*x87 + q13*x88;
  const double x90 = q10*x86 + q11*x85 + q12*x88 - q13*x87;
  const double x91 = q10*x87 - q11*x88 + q12*x85 + q13*x86;
  const double x92 = q10*x88 + q11*x87 - q12*x86 + q13*x85;
  const double x93 = 2*_V2x*q11 + 2*_V2y*q12 + 2*_V2z*q13;
  const double x94 = -2*_V2x*q12 + 2*_V2y*q11 + 2*_V2z*q10;
  const double x95 = 2*_V2x*q13 + 2*_V2y*q10 - 2*_V2z*q11;
  const double x96 = 2*_V2x*q10 - 2*_V2y*q13 + 2*_V2z*q12;
  _jachq->setValue(0, 0, -x13 - x2 + x5);
  _jachq->setValue(0, 1, -x15 + x18 - x21);
  _jachq->setValue(0, 2, -x23 + x25 - x27);
  _jachq->setValue(0, 3, -2*_V1x*x33 - 2*_V1y*x36 - 2*_V1z*x39);
  _jachq->setValue(0, 4, -_V1x*x40 + _V1y*x41 - _V1z*x45);
  _jachq->setValue(0, 5, -_V1x*x47 - _V1y*x40 + _V1z*x46);
  _jachq->setValue(0, 6, _V1x*x48 - _V1y*x49 - _V1z*x40);
  _jachq->setValue(0, 7, x13 + x2 - x5);
  _jachq->setValue(0, 8, x15 - x18 + x21);
  _jachq->setValue(0, 9, x23 - x25 + x27);
  _jachq->setValue(0, 10, 0);
  _jachq->setValue(0, 11, 0);
  _jachq->setValue(0, 12, 0);
  _jachq->setValue(0, 13, 0);
  _jachq->setValue(1, 0, -x50 + x51 - x52);
  _jachq->setValue(1, 1, -x53 + x54 - x55);
  _jachq->setValue(1, 2, -x56 + x57 - x58);
  _jachq->setValue(1, 3, -2*_V2x*x33 - 2*_V2y*x36 - 2*_V2z*x39);
  _jachq->setValue(1, 4, -_V2x*x40 + _V2y*x41 - _V2z*x45);
  _jachq->setValue(1, 5, -_V2x*x47 - _V2y*x40 + _V2z*x46);
  _jachq->setValue(1, 6, _V2x*x48 - _V2y*x49 - _V2z*x40);
  _jachq->setValue(1, 7, x50 - x51 + x52);
  _jachq->setValue(1, 8, x53 - x54 + x55);
  _jachq->setValue(1, 9, x56 - x57 + x58);
  _jachq->setValue(1, 10, 0);
  _jachq->setValue(1, 11, 0);
  _jachq->setValue(1, 12, 0);
  _jachq->setValue(1, 13, 0);

  /* Orientation constraints (H3, H4, H5)
   */
  _jachq->setValue(2, 0, 0);
  _jachq->setValue(2, 1, 0);
  _jachq->setValue(2, 2, 0);
  _jachq->setValue(2, 3, x59*x64 - x61*x71 - x62*x73 - x63*x72
                   - x65*x66 - x67*x68 - x69*x70);
  _jachq->setValue(2, 4, -x59*x66 - x64*x65 - x67*x70 + x68*x69
                   + x74*x75 - x76*x77 + x78*x79);
  _jachq->setValue(2, 5, -x59*x68 - x64*x67 + x65*x70 - x66*x69
                   + x74*x77 + x75*x76 - x79*x80);
  _jachq->setValue(2, 6, -x59*x70 - x64*x69 - x65*x68 + x66*x67
                   + x74*x79 - x75*x78 + x77*x80);
  _jachq->setValue(2, 7, 0);
  _jachq->setValue(2, 8, 0);
  _jachq->setValue(2, 9, 0);
  _jachq->setValue(2, 10, -x64*x81 + x66*x82 + x68*x83 + x70*x84);
  _jachq->setValue(2, 11, -x64*x82 - x66*x81 - x68*x84 + x70*x83);
  _jachq->setValue(2, 12, -x64*x83 + x66*x84 - x68*x81 - x70*x82);
  _jachq->setValue(2, 13, -x64*x84 - x66*x83 + x68*x82 - x70*x81);
  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);
  _jachq->setValue(3, 3, x59*x89 - x65*x90 - x67*x91 - x69*x92
                   - x71*x86 - x72*x88 - x73*x87);
  _jachq->setValue(3, 4, -x59*x90 - x65*x89 - x67*x92 + x69*x91
                   + x75*x93 - x77*x94 + x79*x95);
  _jachq->setValue(3, 5, -x59*x91 + x65*x92 - x67*x89 - x69*x90
                   + x75*x94 + x77*x93 - x79*x96);
  _jachq->setValue(3, 6, -x59*x92 - x65*x91 + x67*x90 - x69*x89
                   - x75*x95 + x77*x96 + x79*x93);
  _jachq->setValue(3, 7, 0);
  _jachq->setValue(3, 8, 0);
  _jachq->setValue(3, 9, 0);
  _jachq->setValue(3, 10, -x81*x89 + x82*x90 + x83*x91 + x84*x92);
  _jachq->setValue(3, 11, -x81*x90 - x82*x89 + x83*x92 - x84*x91);
  _jachq->setValue(3, 12, -x81*x91 - x82*x92 - x83*x89 + x84*x90);
  _jachq->setValue(3, 13, -x81*x92 + x82*x91 - x83*x90 - x84*x89);
}

void CylindricalJointR::Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13)
{
  /*
   * sympy expression:
   *
   * (same as Jd1d2 case but with..)
   * q2 = np.array([1,0,0,0])
   * G2 = P
   *
   * H = [H1,H2,H3,H4]
   * dq = list(G1[1:])+list(q1)
   *
   * jachq = [[h.diff(d) for d in dq] for h in H]
   */

  /* Cylindrical constraints (H1, H2)
   */
  const double x0 = 2*q10*q12;
  const double x1 = 2*q11*q13 + x0;
  const double x2 = 2*q10*q13;
  const double x3 = -2*q11*q12 + x2;
  const double x4 = pow(q10, 2);
  const double x5 = pow(q11, 2);
  const double x6 = pow(q12, 2);
  const double x7 = -x6;
  const double x8 = pow(q13, 2);
  const double x9 = -x8;
  const double x10 = x4 + x5 + x7 + x9;
  const double x11 = 2*q11*q12 + x2;
  const double x12 = 2*q10*q11;
  const double x13 = -2*q12*q13 + x12;
  const double x14 = x4 - x5;
  const double x15 = x14 + x6 + x9;
  const double x16 = 2*q12*q13 + x12;
  const double x17 = -2*q11*q13 + x0;
  const double x18 = x14 + x7 + x8;
  const double x19 = X1 - _P->getValue(0);
  const double x20 = q10*x19;
  const double x21 = Y1 - _P->getValue(1);
  const double x22 = Z1 - _P->getValue(2);
  const double x23 = q12*x22;
  const double x24 = q13*x21 + x20 - x23;
  const double x25 = q10*x21;
  const double x26 = q13*x19;
  const double x27 = q11*x22 + x25 - x26;
  const double x28 = q10*x22;
  const double x29 = q11*x21;
  const double x30 = q12*x19 + x28 - x29;
  const double x31 = 2*q11*x19 + 2*q12*x21 + 2*q13*x22;
  const double x32 = -2*q10*x22 + 2*q11*x21 - 2*q12*x19;
  const double x33 = -Y1 + _P->getValue(1);
  const double x34 = -Z1 + _P->getValue(2);
  const double x35 = -X1 + _P->getValue(0);
  const double x36 = q10*x33 + 2*q11*x34 - q13*x35 - x25 + x26;
  const double x37 = -2*q10*x19 + 2*q12*x22 - 2*q13*x21;
  const double x38 = q10*x34 - q11*x33 + 2*q12*x35 - x28 + x29;
  const double x39 = -2*q10*x21 - 2*q11*x22 + 2*q13*x19;
  const double x40 = q10*x35 - q12*x34 + 2*q13*x33 - x20 + x23;
  const double x41 = 2*_V1x*q10 - 2*_V1y*q13 + 2*_V1z*q12;
  const double x42 = _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12;
  const double x43 = -2*_V1x*q12 + 2*_V1y*q11 + 2*_V1z*q10;
  const double x44 = _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10;
  const double x45 = 2*_V1x*q13 + 2*_V1y*q10 - 2*_V1z*q11;
  const double x46 = _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11;
  const double x47 = _V1x*q11 + _V1y*q12 + _V1z*q13;
  const double x48 = _V1x*q10 - _V1y*q13 + _V1z*q12;
  const double x49 = _V1x*q13 + _V1y*q10 - _V1z*q11;
  const double x50 = -_V1x*q12 + _V1y*q11 + _V1z*q10;
  const double x51 = -q10*x47 + q11*x48 + q12*x49 + q13*x50;
  const double x52 = q10*x48 + q11*x47 + q12*x50 - q13*x49;
  const double x53 = q10*x49 - q11*x50 + q12*x47 + q13*x48;
  const double x54 = q10*x50 + q11*x49 - q12*x48 + q13*x47;
  const double x55 = 2*_V1x*q11 + 2*_V1y*q12 + 2*_V1z*q13;
  const double x56 = 2*_V2x*q10 - 2*_V2y*q13 + 2*_V2z*q12;
  const double x57 = -2*_V2x*q12 + 2*_V2y*q11 + 2*_V2z*q10;
  const double x58 = 2*_V2x*q13 + 2*_V2y*q10 - 2*_V2z*q11;
  const double x59 = _V2x*q11 + _V2y*q12 + _V2z*q13;
  const double x60 = _V2x*q10 - _V2y*q13 + _V2z*q12;
  const double x61 = _V2x*q13 + _V2y*q10 - _V2z*q11;
  const double x62 = -_V2x*q12 + _V2y*q11 + _V2z*q10;
  const double x63 = -q10*x59 + q11*x60 + q12*x61 + q13*x62;
  const double x64 = q10*x60 + q11*x59 + q12*x62 - q13*x61;
  const double x65 = q10*x61 - q11*x62 + q12*x59 + q13*x60;
  const double x66 = q10*x62 + q11*x61 - q12*x60 + q13*x59;
  const double x67 = 2*_V2x*q11 + 2*_V2y*q12 + 2*_V2z*q13;
  _jachq->setValue(0, 0, -_V1x*x10 + _V1y*x3 - _V1z*x1);
  _jachq->setValue(0, 1, -_V1x*x11 - _V1y*x15 + _V1z*x13);
  _jachq->setValue(0, 2, _V1x*x17 - _V1y*x16 - _V1z*x18);
  _jachq->setValue(0, 3, -2*_V1x*x24 - 2*_V1y*x27 - 2*_V1z*x30);
  _jachq->setValue(0, 4, -_V1x*x31 + _V1y*x32 - _V1z*x36);
  _jachq->setValue(0, 5, -_V1x*x38 - _V1y*x31 + _V1z*x37);
  _jachq->setValue(0, 6, _V1x*x39 - _V1y*x40 - _V1z*x31);
  _jachq->setValue(1, 0, -_V2x*x10 + _V2y*x3 - _V2z*x1);
  _jachq->setValue(1, 1, -_V2x*x11 - _V2y*x15 + _V2z*x13);
  _jachq->setValue(1, 2, _V2x*x17 - _V2y*x16 - _V2z*x18);
  _jachq->setValue(1, 3, -2*_V2x*x24 - 2*_V2y*x27 - 2*_V2z*x30);
  _jachq->setValue(1, 4, -_V2x*x31 + _V2y*x32 - _V2z*x36);
  _jachq->setValue(1, 5, -_V2x*x38 - _V2y*x31 + _V2z*x37);
  _jachq->setValue(1, 6, _V2x*x39 - _V2y*x40 - _V2z*x31);

  /* Orientation constraints (H3, H4, H5)
   */
  _jachq->setValue(2, 0, 0);
  _jachq->setValue(2, 1, 0);
  _jachq->setValue(2, 2, 0);
  _jachq->setValue(2, 3, -_cq2q101*x51 - _cq2q102*x52 - _cq2q103*x53 - _cq2q104*x54
                   + x41*x42 + x43*x44 + x45*x46);
  _jachq->setValue(2, 4, _cq2q101*x52 - _cq2q102*x51 - _cq2q103*x54 + _cq2q104*x53
                   + x42*x55 - x43*x46 + x44*x45);
  _jachq->setValue(2, 5, _cq2q101*x53 + _cq2q102*x54 - _cq2q103*x51 - _cq2q104*x52
                   - x41*x44 + x42*x43 + x46*x55);
  _jachq->setValue(2, 6, _cq2q101*x54 - _cq2q102*x53 + _cq2q103*x52 - _cq2q104*x51
                   + x41*x46 - x42*x45 + x44*x55);
  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);
  _jachq->setValue(3, 3, -_cq2q101*x63 - _cq2q102*x64 - _cq2q103*x65 - _cq2q104*x66
                   + x42*x56 + x44*x57 + x46*x58);
  _jachq->setValue(3, 4, _cq2q101*x64 - _cq2q102*x63 - _cq2q103*x66 + _cq2q104*x65
                   + x42*x67 + x44*x58 - x46*x57);
  _jachq->setValue(3, 5, _cq2q101*x65 + _cq2q102*x66 - _cq2q103*x63 - _cq2q104*x64
                   + x42*x57 - x44*x56 + x46*x67);
  _jachq->setValue(3, 6, _cq2q101*x66 - _cq2q102*x65 + _cq2q103*x64 - _cq2q104*x63
                   - x42*x58 + x44*x67 + x46*x56);
}
