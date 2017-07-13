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
 *                   Symbol('_cq2q103'),Symbol('_cq2q104')])
 * G1P0 = np.array([0, Symbol('_G1P0->getValue(0)'), Symbol('_G1P0->getValue(1)'),
 *                  Symbol('_G1P0->getValue(2)')])
 * G2P0 = np.array([0, Symbol('_G2P0->getValue(0)'), Symbol('_G2P0->getValue(1)'),
 *                  Symbol('_G2P0->getValue(2)')])
 * G1 = np.array([0, Symbol('X1'), Symbol('Y1'), Symbol('Z1')])
 * G2 = np.array([0, Symbol('X2'), Symbol('Y2'), Symbol('Z2')])
 * V1 = np.array([0, Symbol('_V1->getValue(0)'), Symbol('_V1->getValue(1)'),
 *                Symbol('_V1->getValue(2)')])
 * V2 = np.array([0, Symbol('_V2->getValue(0)'), Symbol('_V2->getValue(1)'),
 *                Symbol('_V2->getValue(2)')])
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

class quatPos : public ::boost::math::quaternion<double>
{
public: quatPos(const SiconosVector& v)
  : ::boost::math::quaternion<double>
  (0,v.getValue(0),v.getValue(1),v.getValue(2)) {}
};

class quatOrientation : public ::boost::math::quaternion<double>
{
public: quatOrientation(const SiconosVector& v)
  : ::boost::math::quaternion<double>
  (v.getValue(3),v.getValue(4),v.getValue(5),v.getValue(6)) {}
};

// Wrap value in interval [-pi,pi]
static double piwrap(double x)
{
  return fmod(x + 3*M_PI, 2*M_PI) - M_PI;
}

CylindricalJointR::CylindricalJointR(SP::SiconosVector P, SP::SiconosVector A,
                                     bool absoluteRef,
                                     SP::NewtonEulerDS d1, SP::NewtonEulerDS d2)
  : NewtonEulerJointR()
  , _axis0(std11::make_shared<SiconosVector>(3))
  , _G1P0(std11::make_shared<SiconosVector>(3))
{
  _points.resize(1);
  _axes.resize(1);
  setAbsolute(absoluteRef);
  setPoint(0, P);
  setAxis(0, A);
  if (d1)
    setInitialConditions(d1->q(), d2 ? d2->q() : SP::SiconosVector());
}

void CylindricalJointR::setInitialConditions(SP::SiconosVector q1,
                                             SP::SiconosVector q2)
{
  // in the two-DS case, _P is unused.
  _G1P0->zero();
  *_G1P0 = *_points[0];
  _axis0->zero();
  *_axis0 = *_axes[0];

  if (_absoluteRef) {
    _G1P0->setValue(0, _G1P0->getValue(0) - q1->getValue(0));
    _G1P0->setValue(1, _G1P0->getValue(1) - q1->getValue(1));
    _G1P0->setValue(2, _G1P0->getValue(2) - q1->getValue(2));
  }

  computeV1V2FromAxis();
  SP::SiconosVector q2i(new SiconosVector(7));
  q2i->zero();
  q2i->setValue(3, 1);

  if(q2)
    *q2i = *q2;

  ::boost::math::quaternion<double> tmp;
  quatOrientation quat1(*q1);
  quatOrientation quat2(*q2i);
  quatPos quatG10G20_abs(*q1 - *q1);

  // Initial orientation offset
  tmp = 1.0/quat2 * quat1;
  _cq2q101 = tmp.R_component_1();
  _cq2q102 = tmp.R_component_2();
  _cq2q103 = tmp.R_component_3();
  _cq2q104 = tmp.R_component_4();

  // Initial G1-P vector in G1 frame is just P.  Initial G2-P vector
  // in G1 frame is calculated by subtracting (P - G2) in the absolute
  // frame and un-rotating from the q2 frame.
  quatPos quatG1P0(*_G1P0);
  tmp = quat1 * quatG1P0 / quat1;

  SiconosVector P0_abs(3);
  P0_abs.setValue(0, tmp.R_component_2() + q1->getValue(0));
  P0_abs.setValue(1, tmp.R_component_3() + q1->getValue(1));
  P0_abs.setValue(2, tmp.R_component_4() + q1->getValue(2));

  SiconosVector G2_abs(3);
  G2_abs.setValue(0, q2i->getValue(0));
  G2_abs.setValue(1, q2i->getValue(1));
  G2_abs.setValue(2, q2i->getValue(2));

  SiconosVector G2P0_abs(3);
  G2P0_abs = P0_abs - G2_abs;
  quatPos quatG2P0_abs(G2P0_abs);
  tmp = 1.0/quat2 * quatG2P0_abs * quat2;

  _G2P0 = std11::make_shared<SiconosVector>(3);
  _G2P0->setValue(0, tmp.R_component_2());
  _G2P0->setValue(1, tmp.R_component_3());
  _G2P0->setValue(2, tmp.R_component_4());

  /* Compute initial positions/angles of degrees of freedom */
  _previousAngle = 0.0;
  _twistCount = 0;
  _initialAngle = 0.0;
  BlockVector q0(q1, q2);
  SiconosVector tmpy(1);
  this->computehDoF(0, q0, tmpy, 1);
  _initialAngle = tmpy(0);

  _twistCount = 0;
  _previousAngle = 0;
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

  if(q0->numberOfBlocks()>1)
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
  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;

  if (q0.numberOfBlocks()>1)
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

  /* sympy expression:
   *
   * HP = (G2-G1) + rot(G2P0, q2) - rot(G1P0, q1)
   * H1 = np.dot(rot(V1,q1), HP)
   * H2 = np.dot(rot(V2,q1), HP)
   *
   * i.e., the path formed by the 3 vectors, projected onto the plane
   * V1,V2 orthogonal to V, should sum to zero.
   *
   * q2to1 = qmul(q1,qinv(qmul(q2,cq2q10)))
   *
   * H3 = np.dot(rot(V1,q1), q2to1)
   * H4 = np.dot(rot(V2,q1), q2to1)
   *
   * i.e., the rotation from q2 to q1, corrected for initial offset,
   * should be zero when projected onto the plane V1,V2 orthogonal to V
   * and rotated into q1 frame.
   */

  const double x0 = _V1->getValue(0)*q11 + _V1->getValue(1)*q12 + _V1->getValue(2)*q13;
  const double x1 = _V1->getValue(0)*q10 - _V1->getValue(1)*q13 + _V1->getValue(2)*q12;
  const double x2 = _V1->getValue(0)*q13 + _V1->getValue(1)*q10 - _V1->getValue(2)*q11;
  const double x3 = -_V1->getValue(0)*q12 + _V1->getValue(1)*q11 + _V1->getValue(2)*q10;
  const double x4 = -q10*x0 + q11*x1 + q12*x2 + q13*x3;
  const double x5 = _G1P0->getValue(0)*q11 + _G1P0->getValue(1)*q12 + _G1P0->getValue(2)*q13;
  const double x6 = _G2P0->getValue(0)*q21 + _G2P0->getValue(1)*q22 + _G2P0->getValue(2)*q23;
  const double x7 = _G2P0->getValue(0)*q20 - _G2P0->getValue(1)*q23 + _G2P0->getValue(2)*q22;
  const double x8 = _G2P0->getValue(0)*q23 + _G2P0->getValue(1)*q20 - _G2P0->getValue(2)*q21;
  const double x9 = -_G2P0->getValue(0)*q22 + _G2P0->getValue(1)*q21 + _G2P0->getValue(2)*q20;
  const double x10 = _G1P0->getValue(0)*q10 - _G1P0->getValue(1)*q13 + _G1P0->getValue(2)*q12;
  const double x11 = _G1P0->getValue(0)*q13 + _G1P0->getValue(1)*q10 - _G1P0->getValue(2)*q11;
  const double x12 = -_G1P0->getValue(0)*q12 + _G1P0->getValue(1)*q11 + _G1P0->getValue(2)*q10;
  const double x13 = q10*x5 - q11*x10 - q12*x11 - q13*x12 - q20*x6 + q21*x7 + q22*x8 + q23*x9;
  const double x14 = q10*x1 + q11*x0 + q12*x3 - q13*x2;
  const double x15 = X1 - X2 + q10*x10 + q11*x5 + q12*x12 - q13*x11 - q20*x7 - q21*x6 - q22*x9 + q23*x8;
  const double x16 = q10*x3 + q11*x2 - q12*x1 + q13*x0;
  const double x17 = Z1 - Z2 + q10*x12 + q11*x11 - q12*x10 + q13*x5 - q20*x9 - q21*x8 + q22*x7 - q23*x6;
  const double x18 = q10*x2 - q11*x3 + q12*x0 + q13*x1;
  const double x19 = Y1 - Y2 + q10*x11 - q11*x12 + q12*x5 + q13*x10 - q20*x8 + q21*x9 - q22*x6 - q23*x7;
  const double x20 = _V2->getValue(0)*q11 + _V2->getValue(1)*q12 + _V2->getValue(2)*q13;
  const double x21 = _V2->getValue(0)*q10 - _V2->getValue(1)*q13 + _V2->getValue(2)*q12;
  const double x22 = _V2->getValue(0)*q13 + _V2->getValue(1)*q10 - _V2->getValue(2)*q11;
  const double x23 = -_V2->getValue(0)*q12 + _V2->getValue(1)*q11 + _V2->getValue(2)*q10;
  const double x24 = -q10*x20 + q11*x21 + q12*x22 + q13*x23;
  const double x25 = q10*x21 + q11*x20 + q12*x23 - q13*x22;
  const double x26 = q10*x23 + q11*x22 - q12*x21 + q13*x20;
  const double x27 = q10*x22 - q11*x23 + q12*x20 + q13*x21;
  const double x28 = _cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22;
  const double x29 = -_cq2q101*q20 + _cq2q102*q21 + _cq2q103*q22 + _cq2q104*q23;
  const double x30 = _cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20;
  const double x31 = _cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21;
  const double x32 = q10*x28 + q11*x29 + q12*x30 - q13*x31;
  const double x33 = q10*x30 + q11*x31 - q12*x28 + q13*x29;
  const double x34 = q10*x31 - q11*x30 + q12*x29 + q13*x28;
  const double x35 = q10*x29 - q11*x28 - q12*x31 - q13*x30;
  y.setValue(0, x13*x4 - x14*x15 - x16*x17 - x18*x19);
  y.setValue(1, x13*x24 - x15*x25 - x17*x26 - x19*x27);
  y.setValue(2, -x14*x32 - x16*x33 - x18*x34 + x35*x4);
  y.setValue(3, x24*x35 - x25*x32 - x26*x33 - x27*x34);
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
  const double x0 = _V1->getValue(0)*q11 + _V1->getValue(1)*q12 + _V1->getValue(2)*q13;
  const double x1 = q11*x0;
  const double x2 = _V1->getValue(0)*q13 + _V1->getValue(1)*q10 - _V1->getValue(2)*q11;
  const double x3 = q13*x2;
  const double x4 = _V1->getValue(0)*q10 - _V1->getValue(1)*q13 + _V1->getValue(2)*q12;
  const double x5 = q10*x4;
  const double x6 = -_V1->getValue(0)*q12 + _V1->getValue(1)*q11 + _V1->getValue(2)*q10;
  const double x7 = q12*x6;
  const double x8 = q11*x6;
  const double x9 = q12*x0;
  const double x10 = q10*x2;
  const double x11 = q13*x4;
  const double x12 = q12*x4;
  const double x13 = q13*x0;
  const double x14 = q10*x6;
  const double x15 = q11*x2;
  const double x16 = _G1P0->getValue(0)*q10 - _G1P0->getValue(1)*q13 + _G1P0->getValue(2)*q12;
  const double x17 = x1 - x3 + x5 + x7;
  const double x18 = -_G1P0->getValue(0)*q12 + _G1P0->getValue(1)*q11 + _G1P0->getValue(2)*q10;
  const double x19 = -x12 + x13 + x14 + x15;
  const double x20 = _G1P0->getValue(0)*q13 + _G1P0->getValue(1)*q10 - _G1P0->getValue(2)*q11;
  const double x21 = x10 + x11 - x8 + x9;
  const double x22 = _G1P0->getValue(0)*q11 + _G1P0->getValue(1)*q12 + _G1P0->getValue(2)*q13;
  const double x23 = _G2P0->getValue(0)*q21 + _G2P0->getValue(1)*q22 + _G2P0->getValue(2)*q23;
  const double x24 = _G2P0->getValue(0)*q23 + _G2P0->getValue(1)*q20 - _G2P0->getValue(2)*q21;
  const double x25 = _G2P0->getValue(0)*q20 - _G2P0->getValue(1)*q23 + _G2P0->getValue(2)*q22;
  const double x26 = -_G2P0->getValue(0)*q22 + _G2P0->getValue(1)*q21 + _G2P0->getValue(2)*q20;
  const double x27 = X1 - X2 + q10*x16 + q11*x22 + q12*x18 - q13*x20 - q20*x25 - q21*x23 - q22*x26 + q23*x24;
  const double x28 = Z1 - Z2 + q10*x18 + q11*x20 - q12*x16 + q13*x22 - q20*x26 - q21*x24 + q22*x25 - q23*x23;
  const double x29 = Y1 - Y2 + q10*x20 - q11*x18 + q12*x22 + q13*x16 - q20*x24 + q21*x26 - q22*x23 - q23*x25;
  const double x30 = _V2->getValue(0)*q11 + _V2->getValue(1)*q12 + _V2->getValue(2)*q13;
  const double x31 = q11*x30;
  const double x32 = _V2->getValue(0)*q13 + _V2->getValue(1)*q10 - _V2->getValue(2)*q11;
  const double x33 = q13*x32;
  const double x34 = _V2->getValue(0)*q10 - _V2->getValue(1)*q13 + _V2->getValue(2)*q12;
  const double x35 = q10*x34;
  const double x36 = -_V2->getValue(0)*q12 + _V2->getValue(1)*q11 + _V2->getValue(2)*q10;
  const double x37 = q12*x36;
  const double x38 = q11*x36;
  const double x39 = q12*x30;
  const double x40 = q10*x32;
  const double x41 = q13*x34;
  const double x42 = q12*x34;
  const double x43 = q13*x30;
  const double x44 = q10*x36;
  const double x45 = q11*x32;
  const double x46 = x31 - x33 + x35 + x37;
  const double x47 = -x42 + x43 + x44 + x45;
  const double x48 = -x38 + x39 + x40 + x41;
  const double x49 = -_cq2q101*q20 + _cq2q102*q21 + _cq2q103*q22 + _cq2q104*q23;
  const double x50 = -q10*x0 + q11*x4 + q12*x2 + q13*x6;
  const double x51 = _cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22;
  const double x52 = _cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21;
  const double x53 = _cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20;
  const double x54 = 2*q10*x51 + 2*q11*x49 + 2*q12*x53 - 2*q13*x52;
  const double x55 = 2*q10*x53 + 2*q11*x52 - 2*q12*x51 + 2*q13*x49;
  const double x56 = 2*q10*x52 - 2*q11*x53 + 2*q12*x49 + 2*q13*x51;
  const double x57 = 2*_V1->getValue(0)*q11 + 2*_V1->getValue(1)*q12 + 2*_V1->getValue(2)*q13;
  const double x58 = -q10*x51 - q11*x49 - q12*x53 + q13*x52;
  const double x59 = -2*_V1->getValue(0)*q12 + 2*_V1->getValue(1)*q11 + 2*_V1->getValue(2)*q10;
  const double x60 = -q10*x52 + q11*x53 - q12*x49 - q13*x51;
  const double x61 = 2*_V1->getValue(0)*q13 + 2*_V1->getValue(1)*q10 - 2*_V1->getValue(2)*q11;
  const double x62 = -q10*x53 - q11*x52 + q12*x51 - q13*x49;
  const double x63 = 2*_V1->getValue(0)*q10 - 2*_V1->getValue(1)*q13 + 2*_V1->getValue(2)*q12;
  const double x64 = _cq2q101*q10 + _cq2q102*q11 + _cq2q103*q12 + _cq2q104*q13;
  const double x65 = _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12;
  const double x66 = _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11;
  const double x67 = _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10;
  const double x68 = -q10*x30 + q11*x34 + q12*x32 + q13*x36;
  const double x69 = 2*_V2->getValue(0)*q11 + 2*_V2->getValue(1)*q12 + 2*_V2->getValue(2)*q13;
  const double x70 = -2*_V2->getValue(0)*q12 + 2*_V2->getValue(1)*q11 + 2*_V2->getValue(2)*q10;
  const double x71 = 2*_V2->getValue(0)*q13 + 2*_V2->getValue(1)*q10 - 2*_V2->getValue(2)*q11;
  const double x72 = 2*_V2->getValue(0)*q10 - 2*_V2->getValue(1)*q13 + 2*_V2->getValue(2)*q12;
  _jachq->setValue(0, 0, -x1 + x3 - x5 - x7);
  _jachq->setValue(0, 1, -x10 - x11 + x8 - x9);
  _jachq->setValue(0, 2, x12 - x13 - x14 - x15);
  _jachq->setValue(0, 3, -2*x16*x17 - 2*x18*x19 - 2*x2*x29 - 2*x20*x21 - 2*x27*x4 - 2*x28*x6);
  _jachq->setValue(0, 4, -2*x0*x27 - 2*x17*x22 + 2*x18*x21 - 2*x19*x20 - 2*x2*x28 + 2*x29*x6);
  _jachq->setValue(0, 5, -2*x0*x29 + 2*x16*x19 - 2*x17*x18 - 2*x21*x22 - 2*x27*x6 + 2*x28*x4);
  _jachq->setValue(0, 6, -2*x0*x28 - 2*x16*x21 + 2*x17*x20 - 2*x19*x22 + 2*x2*x27 - 2*x29*x4);
  _jachq->setValue(0, 7, x17);
  _jachq->setValue(0, 8, x21);
  _jachq->setValue(0, 9, x19);
  _jachq->setValue(0, 10, 2*x17*x25 + 2*x19*x26 + 2*x21*x24);
  _jachq->setValue(0, 11, 2*x17*x23 + 2*x19*x24 - 2*x21*x26);
  _jachq->setValue(0, 12, 2*x17*x26 - 2*x19*x25 + 2*x21*x23);
  _jachq->setValue(0, 13, -2*x17*x24 + 2*x19*x23 + 2*x21*x25);
  _jachq->setValue(1, 0, -x31 + x33 - x35 - x37);
  _jachq->setValue(1, 1, x38 - x39 - x40 - x41);
  _jachq->setValue(1, 2, x42 - x43 - x44 - x45);
  _jachq->setValue(1, 3, -2*x16*x46 - 2*x18*x47 - 2*x20*x48 - 2*x27*x34 - 2*x28*x36 - 2*x29*x32);
  _jachq->setValue(1, 4, 2*x18*x48 - 2*x20*x47 - 2*x22*x46 - 2*x27*x30 - 2*x28*x32 + 2*x29*x36);
  _jachq->setValue(1, 5, 2*x16*x47 - 2*x18*x46 - 2*x22*x48 - 2*x27*x36 + 2*x28*x34 - 2*x29*x30);
  _jachq->setValue(1, 6, -2*x16*x48 + 2*x20*x46 - 2*x22*x47 + 2*x27*x32 - 2*x28*x30 - 2*x29*x34);
  _jachq->setValue(1, 7, x46);
  _jachq->setValue(1, 8, x48);
  _jachq->setValue(1, 9, x47);
  _jachq->setValue(1, 10, 2*x24*x48 + 2*x25*x46 + 2*x26*x47);
  _jachq->setValue(1, 11, 2*x23*x46 + 2*x24*x47 - 2*x26*x48);
  _jachq->setValue(1, 12, 2*x23*x48 - 2*x25*x47 + 2*x26*x46);
  _jachq->setValue(1, 13, 2*x23*x47 - 2*x24*x46 + 2*x25*x48);

  /* Orientation constraints (H3, H4, H5)
   */
  _jachq->setValue(2, 0, 0);
  _jachq->setValue(2, 1, 0);
  _jachq->setValue(2, 2, 0);
  _jachq->setValue(2, 3, -x17*x51 - x19*x53 - x2*x56 - x21*x52 - x4*x54 + x49*x50 - x55*x6);
  _jachq->setValue(2, 4, -x17*x49 - x19*x52 + x21*x53 - x50*x51 + x57*x58 - x59*x60 + x61*x62);
  _jachq->setValue(2, 5, -x17*x53 + x19*x51 - x21*x49 - x50*x52 + x57*x60 + x58*x59 - x62*x63);
  _jachq->setValue(2, 6, x17*x52 - x19*x49 - x21*x51 - x50*x53 + x57*x62 - x58*x61 + x60*x63);
  _jachq->setValue(2, 7, 0);
  _jachq->setValue(2, 8, 0);
  _jachq->setValue(2, 9, 0);
  _jachq->setValue(2, 10, x17*x65 + x19*x67 + x21*x66 - x50*x64);
  _jachq->setValue(2, 11, -x17*x64 + x19*x66 - x21*x67 - x50*x65);
  _jachq->setValue(2, 12, x17*x67 - x19*x65 - x21*x64 - x50*x66);
  _jachq->setValue(2, 13, -x17*x66 - x19*x64 + x21*x65 - x50*x67);
  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);
  _jachq->setValue(3, 3, -x32*x56 - x34*x54 - x36*x55 - x46*x51 - x47*x53 - x48*x52 + x49*x68);
  _jachq->setValue(3, 4, -x46*x49 - x47*x52 + x48*x53 - x51*x68 + x58*x69 - x60*x70 + x62*x71);
  _jachq->setValue(3, 5, -x46*x53 + x47*x51 - x48*x49 - x52*x68 + x58*x70 + x60*x69 - x62*x72);
  _jachq->setValue(3, 6, x46*x52 - x47*x49 - x48*x51 - x53*x68 - x58*x71 + x60*x72 + x62*x69);
  _jachq->setValue(3, 7, 0);
  _jachq->setValue(3, 8, 0);
  _jachq->setValue(3, 9, 0);
  _jachq->setValue(3, 10, x46*x65 + x47*x67 + x48*x66 - x64*x68);
  _jachq->setValue(3, 11, -x46*x64 + x47*x66 - x48*x67 - x65*x68);
  _jachq->setValue(3, 12, x46*x67 - x47*x65 - x48*x64 - x66*x68);
  _jachq->setValue(3, 13, -x46*x66 - x47*x64 + x48*x65 - x67*x68);
}

void CylindricalJointR::Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13)
{
  /*
   * sympy expression:
   *
   * (same as Jd1d2 case but with..)
   * q2 = np.array([1,0,0,0])
   * G2 = np.array([0,0,0,0])
   *
   * H = [H1,H2,H3,H4]
   * dq = list(G1[1:])+list(q1)
   *
   * jachq = [[h.diff(d) for d in dq] for h in H]
   */

  /* Cylindrical constraints (H1, H2)
   */
  const double x0 = _V1->getValue(0)*q11 + _V1->getValue(1)*q12 + _V1->getValue(2)*q13;
  const double x1 = q11*x0;
  const double x2 = _V1->getValue(0)*q13 + _V1->getValue(1)*q10 - _V1->getValue(2)*q11;
  const double x3 = q13*x2;
  const double x4 = _V1->getValue(0)*q10 - _V1->getValue(1)*q13 + _V1->getValue(2)*q12;
  const double x5 = q10*x4;
  const double x6 = -_V1->getValue(0)*q12 + _V1->getValue(1)*q11 + _V1->getValue(2)*q10;
  const double x7 = q12*x6;
  const double x8 = q11*x6;
  const double x9 = q12*x0;
  const double x10 = q10*x2;
  const double x11 = q13*x4;
  const double x12 = q12*x4;
  const double x13 = q13*x0;
  const double x14 = q10*x6;
  const double x15 = q11*x2;
  const double x16 = _G1P0->getValue(0)*q10 - _G1P0->getValue(1)*q13 + _G1P0->getValue(2)*q12;
  const double x17 = x1 - x3 + x5 + x7;
  const double x18 = -_G1P0->getValue(0)*q12 + _G1P0->getValue(1)*q11 + _G1P0->getValue(2)*q10;
  const double x19 = -x12 + x13 + x14 + x15;
  const double x20 = _G1P0->getValue(0)*q13 + _G1P0->getValue(1)*q10 - _G1P0->getValue(2)*q11;
  const double x21 = x10 + x11 - x8 + x9;
  const double x22 = _G1P0->getValue(0)*q11 + _G1P0->getValue(1)*q12 + _G1P0->getValue(2)*q13;
  const double x23 = X1 - 1.0*_G2P0->getValue(0) + q10*x16 + q11*x22 + q12*x18 - q13*x20;
  const double x24 = Z1 - 1.0*_G2P0->getValue(2) + q10*x18 + q11*x20 - q12*x16 + q13*x22;
  const double x25 = Y1 - 1.0*_G2P0->getValue(1) + q10*x20 - q11*x18 + q12*x22 + q13*x16;
  const double x26 = _V2->getValue(0)*q11 + _V2->getValue(1)*q12 + _V2->getValue(2)*q13;
  const double x27 = q11*x26;
  const double x28 = _V2->getValue(0)*q13 + _V2->getValue(1)*q10 - _V2->getValue(2)*q11;
  const double x29 = q13*x28;
  const double x30 = _V2->getValue(0)*q10 - _V2->getValue(1)*q13 + _V2->getValue(2)*q12;
  const double x31 = q10*x30;
  const double x32 = -_V2->getValue(0)*q12 + _V2->getValue(1)*q11 + _V2->getValue(2)*q10;
  const double x33 = q12*x32;
  const double x34 = q11*x32;
  const double x35 = q12*x26;
  const double x36 = q10*x28;
  const double x37 = q13*x30;
  const double x38 = q12*x30;
  const double x39 = q13*x26;
  const double x40 = q10*x32;
  const double x41 = q11*x28;
  const double x42 = x27 - x29 + x31 + x33;
  const double x43 = -x38 + x39 + x40 + x41;
  const double x44 = -x34 + x35 + x36 + x37;
  const double x45 = 2*_V1->getValue(0)*q10 - 2*_V1->getValue(1)*q13 + 2*_V1->getValue(2)*q12;
  const double x46 = _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12;
  const double x47 = -2*_V1->getValue(0)*q12 + 2*_V1->getValue(1)*q11 + 2*_V1->getValue(2)*q10;
  const double x48 = _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10;
  const double x49 = 2*_V1->getValue(0)*q13 + 2*_V1->getValue(1)*q10 - 2*_V1->getValue(2)*q11;
  const double x50 = _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11;
  const double x51 = -q10*x0 + q11*x4 + q12*x2 + q13*x6;
  const double x52 = 2*_V1->getValue(0)*q11 + 2*_V1->getValue(1)*q12 + 2*_V1->getValue(2)*q13;
  const double x53 = 2*_V2->getValue(0)*q10 - 2*_V2->getValue(1)*q13 + 2*_V2->getValue(2)*q12;
  const double x54 = -2*_V2->getValue(0)*q12 + 2*_V2->getValue(1)*q11 + 2*_V2->getValue(2)*q10;
  const double x55 = 2*_V2->getValue(0)*q13 + 2*_V2->getValue(1)*q10 - 2*_V2->getValue(2)*q11;
  const double x56 = -q10*x26 + q11*x30 + q12*x28 + q13*x32;
  const double x57 = 2*_V2->getValue(0)*q11 + 2*_V2->getValue(1)*q12 + 2*_V2->getValue(2)*q13;
  _jachq->setValue(0, 0, -x1 + x3 - x5 - x7);
  _jachq->setValue(0, 1, -x10 - x11 + x8 - x9);
  _jachq->setValue(0, 2, x12 - x13 - x14 - x15);
  _jachq->setValue(0, 3, -2*x16*x17 - 2*x18*x19 - 2*x2*x25 - 2*x20*x21 - 2*x23*x4 - 2*x24*x6);
  _jachq->setValue(0, 4, -2*x0*x23 - 2*x17*x22 + 2*x18*x21 - 2*x19*x20 - 2*x2*x24 + 2*x25*x6);
  _jachq->setValue(0, 5, -2*x0*x25 + 2*x16*x19 - 2*x17*x18 - 2*x21*x22 - 2*x23*x6 + 2*x24*x4);
  _jachq->setValue(0, 6, -2*x0*x24 - 2*x16*x21 + 2*x17*x20 - 2*x19*x22 + 2*x2*x23 - 2*x25*x4);
  _jachq->setValue(1, 0, -x27 + x29 - x31 - x33);
  _jachq->setValue(1, 1, x34 - x35 - x36 - x37);
  _jachq->setValue(1, 2, x38 - x39 - x40 - x41);
  _jachq->setValue(1, 3, -2*x16*x42 - 2*x18*x43 - 2*x20*x44 - 2*x23*x30 - 2*x24*x32 - 2*x25*x28);
  _jachq->setValue(1, 4, 2*x18*x44 - 2*x20*x43 - 2*x22*x42 - 2*x23*x26 - 2*x24*x28 + 2*x25*x32);
  _jachq->setValue(1, 5, 2*x16*x43 - 2*x18*x42 - 2*x22*x44 - 2*x23*x32 + 2*x24*x30 - 2*x25*x26);
  _jachq->setValue(1, 6, -2*x16*x44 + 2*x20*x42 - 2*x22*x43 + 2*x23*x28 - 2*x24*x26 - 2*x25*x30);

  /* Orientation constraints (H3, H4, H5)
   */
  _jachq->setValue(2, 0, 0);
  _jachq->setValue(2, 1, 0);
  _jachq->setValue(2, 2, 0);
  _jachq->setValue(2, 3, -_cq2q101*x51 - _cq2q102*x17 - _cq2q103*x21 - _cq2q104*x19 + x45*x46 + x47*x48 + x49*x50);
  _jachq->setValue(2, 4, _cq2q101*x17 - _cq2q102*x51 - _cq2q103*x19 + _cq2q104*x21 + x46*x52 - x47*x50 + x48*x49);
  _jachq->setValue(2, 5, _cq2q101*x21 + _cq2q102*x19 - _cq2q103*x51 - _cq2q104*x17 - x45*x48 + x46*x47 + x50*x52);
  _jachq->setValue(2, 6, _cq2q101*x19 - _cq2q102*x21 + _cq2q103*x17 - _cq2q104*x51 + x45*x50 - x46*x49 + x48*x52);
  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);
  _jachq->setValue(3, 3, -_cq2q101*x56 - _cq2q102*x42 - _cq2q103*x44 - _cq2q104*x43 + x46*x53 + x48*x54 + x50*x55);
  _jachq->setValue(3, 4, _cq2q101*x42 - _cq2q102*x56 - _cq2q103*x43 + _cq2q104*x44 + x46*x57 + x48*x55 - x50*x54);
  _jachq->setValue(3, 5, _cq2q101*x44 + _cq2q102*x43 - _cq2q103*x56 - _cq2q104*x42 + x46*x54 - x48*x53 + x50*x57);
  _jachq->setValue(3, 6, _cq2q101*x43 - _cq2q102*x44 + _cq2q103*x42 - _cq2q104*x56 - x46*x55 + x48*x57 + x50*x53);
}

/** Compute the vector of linear and angular positions of the free axes */
void CylindricalJointR::computehDoF(double time, BlockVector& q0, SiconosVector& y,
                                    unsigned int axis)
{
  if (axis > 1)
    return;

  SP::SiconosVector q1 = (q0.getAllVect())[0];
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

  if (q0.numberOfBlocks()>1)
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

  /*
   * sympy expression:
   *
   * HPos = np.dot(rot(A,q1), HP)
   * HAngDot = np.dot(rot(A,q1), q2to1)
   * HAng = 2*sympy.atan2(q2to1[0], HAngDot)
   *
   * # We call atan2 ourselves
   * HDoF = [HPos,HAngDot,q2to1[0]]
   *
   * # But calculate the full derivative symbolically
   * jachqDoF = [[h.diff(d) for d in dq] for h in [HPos,HAng]]
   */

  /* Pre-calculations */
  const double x0 = _axis0->getValue(0)*q11 + _axis0->getValue(1)*q12 + _axis0->getValue(2)*q13;
  const double x1 = _axis0->getValue(0)*q10 - _axis0->getValue(1)*q13 + _axis0->getValue(2)*q12;
  const double x2 = _axis0->getValue(0)*q13 + _axis0->getValue(1)*q10 - _axis0->getValue(2)*q11;
  const double x3 = -_axis0->getValue(0)*q12 + _axis0->getValue(1)*q11 + _axis0->getValue(2)*q10;
  const double x4 = -q10*x0 + q11*x1 + q12*x2 + q13*x3;
  const double x5 = _G1P0->getValue(0)*q11 + _G1P0->getValue(1)*q12 + _G1P0->getValue(2)*q13;
  const double x6 = _G2P0->getValue(0)*q21 + _G2P0->getValue(1)*q22 + _G2P0->getValue(2)*q23;
  const double x7 = _G2P0->getValue(0)*q20 - _G2P0->getValue(1)*q23 + _G2P0->getValue(2)*q22;
  const double x8 = _G2P0->getValue(0)*q23 + _G2P0->getValue(1)*q20 - _G2P0->getValue(2)*q21;
  const double x9 = -_G2P0->getValue(0)*q22 + _G2P0->getValue(1)*q21 + _G2P0->getValue(2)*q20;
  const double x10 = _G1P0->getValue(0)*q10 - _G1P0->getValue(1)*q13 + _G1P0->getValue(2)*q12;
  const double x11 = _G1P0->getValue(0)*q13 + _G1P0->getValue(1)*q10 - _G1P0->getValue(2)*q11;
  const double x12 = -_G1P0->getValue(0)*q12 + _G1P0->getValue(1)*q11 + _G1P0->getValue(2)*q10;
  const double x13 = q10*x1 + q11*x0 + q12*x3 - q13*x2;
  const double x14 = q10*x3 + q11*x2 - q12*x1 + q13*x0;
  const double x15 = q10*x2 - q11*x3 + q12*x0 + q13*x1;
  const double x16 = _cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22;
  const double x17 = -_cq2q101*q20 + _cq2q102*q21 + _cq2q103*q22 + _cq2q104*q23;
  const double x18 = _cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20;
  const double x19 = _cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21;
  const double x20 = q10*x17;
  const double x21 = q11*x16;
  const double x22 = q12*x19;
  const double x23 = q13*x18;

  /* Linear axis */
  unsigned int i = 0;
  if (axis+i == 0 && y.size() > i)
  {
    /* axis 0 = Linear position */
    y.setValue(i, -x13*(X1 - X2 + q10*x10 + q11*x5 + q12*x12 - q13*x11
                        - q20*x7 - q21*x6 - q22*x9 + q23*x8)
                  - x14*(Z1 - Z2 + q10*x12 + q11*x11 - q12*x10 + q13*x5
                         - q20*x9 - q21*x8 + q22*x7 - q23*x6)
                  - x15*(Y1 - Y2 + q10*x11 - q11*x12 + q12*x5 + q13*x10
                         - q20*x8 + q21*x9 - q22*x6 - q23*x7)
                  + x4*(q10*x5 - q11*x10 - q12*x11 - q13*x12 - q20*x6
                        + q21*x7 + q22*x8 + q23*x9));

    i ++;
  }

  /* Rotational axis */
  if (axis+i == 1 && y.size() > i)
  {
    /* The dot product with the rotational free axis taking into
     * account original angular difference. */
    double Adot2to1 = -x13*(q10*x16 + q11*x17 + q12*x18 - q13*x19)
                     - x14*(q10*x18 + q11*x19 - q12*x16 + q13*x17)
                     - x15*(q10*x19 - q11*x18 + q12*x17 + q13*x16)
                     + x4*(x20 - x21 - x22 - x23);

    // We only need the w part of the quaternion for atan2
    // (see rot2to1() in PivotJointR.cpp, and sympy expression above.)
    double rot2to1w = -x20 + x21 + x22 + x23;

    // In case of joint constraints, it's okay to use dot product=0, but
    // in the case of the free axis we must measure the actual angle
    // using atan2 so that stops can be placed correctly.
    double wrappedAngle = piwrap(2*atan2(rot2to1w, Adot2to1) - _initialAngle);

    // Count the number of twists around the angle, and report the
    // unwrapped angle.  Needed to implement joint stops near pi.
    if (wrappedAngle < -M_PI*3/4 && _previousAngle > M_PI*3/4)
      _twistCount ++;
    else if (wrappedAngle > M_PI*3/4 && _previousAngle < -M_PI*3/4)
      _twistCount --;
    _previousAngle = wrappedAngle;
    double unwrappedAngle = wrappedAngle + 2*M_PI*_twistCount;

    /* axis 1 = Angular position */
    y.setValue(i, unwrappedAngle);
    i ++;
  }
}

/** Compute the jacobian of linear and angular DoF with respect to some q */
void CylindricalJointR::computeJachqDoF(double time, Interaction& inter,
                                        SP::BlockVector q0, SimpleMatrix& jachq,
                                        unsigned int axis)
{
  if (axis > 1)
    return;

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

  if (q0->numberOfBlocks()>1)
  {
    SP::SiconosVector q2 = (q0->getAllVect())[1];
    X2 = q2->getValue(0);
    Y2 = q2->getValue(1);
    Z2 = q2->getValue(2);
    q20 = q2->getValue(3);
    q21 = q2->getValue(4);
    q22 = q2->getValue(5);
    q23 = q2->getValue(6);
  }

  /*
   * sympy expression:
   *
   * HPos = np.dot(rot(A,q1), HP)
   * HAngDot = np.dot(rot(A,q1), q2to1)
   * HAng = 2*sympy.atan2(q2to1[0], HAngDot)
   *
   * # We call atan2 ourselves
   * HDoF = [HPos,HAngDot,q2to1[0]]
   *
   * # But calculate the full derivative symbolically
   * jachqDoF = [[h.diff(d) for d in dq] for h in [HPos,HAng]]
   */

  /* Pre-calculations */
  const double x0 = _axis0->getValue(0)*q11 + _axis0->getValue(1)*q12 + _axis0->getValue(2)*q13;
  const double x1 = q11*x0;
  const double x2 = _axis0->getValue(0)*q13 + _axis0->getValue(1)*q10 - _axis0->getValue(2)*q11;
  const double x3 = q13*x2;
  const double x4 = _axis0->getValue(0)*q10 - _axis0->getValue(1)*q13 + _axis0->getValue(2)*q12;
  const double x5 = q10*x4;
  const double x6 = -_axis0->getValue(0)*q12 + _axis0->getValue(1)*q11 + _axis0->getValue(2)*q10;
  const double x7 = q12*x6;
  const double x8 = q11*x6;
  const double x9 = q12*x0;
  const double x10 = q10*x2;
  const double x11 = q13*x4;
  const double x12 = q12*x4;
  const double x13 = q13*x0;
  const double x14 = q10*x6;
  const double x15 = q11*x2;
  const double x16 = _G1P0->getValue(0)*q10 - _G1P0->getValue(1)*q13 + _G1P0->getValue(2)*q12;
  const double x17 = x1 - x3 + x5 + x7;
  const double x18 = -_G1P0->getValue(0)*q12 + _G1P0->getValue(1)*q11 + _G1P0->getValue(2)*q10;
  const double x19 = -x12 + x13 + x14 + x15;
  const double x20 = _G1P0->getValue(0)*q13 + _G1P0->getValue(1)*q10 - _G1P0->getValue(2)*q11;
  const double x21 = x10 + x11 - x8 + x9;
  const double x22 = _G1P0->getValue(0)*q11 + _G1P0->getValue(1)*q12 + _G1P0->getValue(2)*q13;
  const double x23 = _G2P0->getValue(0)*q21 + _G2P0->getValue(1)*q22 + _G2P0->getValue(2)*q23;
  const double x24 = _G2P0->getValue(0)*q23 + _G2P0->getValue(1)*q20 - _G2P0->getValue(2)*q21;
  const double x25 = _G2P0->getValue(0)*q20 - _G2P0->getValue(1)*q23 + _G2P0->getValue(2)*q22;
  const double x26 = -_G2P0->getValue(0)*q22 + _G2P0->getValue(1)*q21 + _G2P0->getValue(2)*q20;
  const double x27 = X1 - X2 + q10*x16 + q11*x22 + q12*x18 - q13*x20 - q20*x25 - q21*x23 - q22*x26 + q23*x24;
  const double x28 = Z1 - Z2 + q10*x18 + q11*x20 - q12*x16 + q13*x22 - q20*x26 - q21*x24 + q22*x25 - q23*x23;
  const double x29 = Y1 - Y2 + q10*x20 - q11*x18 + q12*x22 + q13*x16 - q20*x24 + q21*x26 - q22*x23 - q23*x25;
  const double x30 = _cq2q103*q23;
  const double x31 = _cq2q101*q21;
  const double x32 = _cq2q102*q20;
  const double x33 = _cq2q104*q22;
  const double x34 = _cq2q104*q21;
  const double x35 = _cq2q101*q22;
  const double x36 = _cq2q102*q23;
  const double x37 = _cq2q103*q20;
  const double x38 = _cq2q102*q22;
  const double x39 = _cq2q101*q23;
  const double x40 = _cq2q103*q21;
  const double x41 = _cq2q104*q20;
  const double x42 = _cq2q101*q20;
  const double x43 = _cq2q102*q21;
  const double x44 = _cq2q103*q22;
  const double x45 = _cq2q104*q23;
  const double x46 = -q10*(x42 - x43 - x44 - x45) + q11*(x30 - x31 - x32 - x33) + q12*(x34 - x35 - x36 - x37) + q13*(x38 - x39 - x40 - x41);
  const double x47 = -x34 + x35 + x36 + x37;
  const double x48 = q13*x47;
  const double x49 = -x30 + x31 + x32 + x33;
  const double x50 = q10*x49;
  const double x51 = -x42 + x43 + x44 + x45;
  const double x52 = q11*x51;
  const double x53 = -x38 + x39 + x40 + x41;
  const double x54 = q12*x53;
  const double x55 = x48 - x50 - x52 - x54;
  const double x56 = q12*x49;
  const double x57 = q10*x53;
  const double x58 = q11*x47;
  const double x59 = q13*x51;
  const double x60 = x56 - x57 - x58 - x59;
  const double x61 = q11*x53;
  const double x62 = q10*x47;
  const double x63 = q12*x51;
  const double x64 = q13*x49;
  const double x65 = x61 - x62 - x63 - x64;
  const double x66 = -q10*x0 + q11*x4 + q12*x2 + q13*x6;
  const double x67 = x17*x55 + x19*x60 + x21*x65 + x46*x66;
  const double x68 = q10*x51;
  const double x69 = q11*x49;
  const double x70 = q12*x47;
  const double x71 = q13*x53;
  const double x72 = x68 - x69 - x70 - x71;
  const double x73 = 2*_axis0->getValue(0)*q10 - 2*_axis0->getValue(1)*q13 + 2*_axis0->getValue(2)*q12;
  const double x74 = -2*_axis0->getValue(0)*q12 + 2*_axis0->getValue(1)*q11 + 2*_axis0->getValue(2)*q10;
  const double x75 = 2*_axis0->getValue(0)*q13 + 2*_axis0->getValue(1)*q10 - 2*_axis0->getValue(2)*q11;
  const double x76 = -x17*(-x48 + x50 + x52 + x54) - x19*(-x56 + x57 + x58 + x59) - x21*(-x61 + x62 + x63 + x64) + x66*x72;
  const double x77 = 2/(pow(x72, 2) + pow(x76, 2));
  const double x78 = -x68 + x69 + x70 + x71;
  const double x79 = 2*_axis0->getValue(0)*q11 + 2*_axis0->getValue(1)*q12 + 2*_axis0->getValue(2)*q13;
  const double x80 = _cq2q101*q10 + _cq2q102*q11 + _cq2q103*q12 + _cq2q104*q13;
  const double x81 = _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12;
  const double x82 = _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11;
  const double x83 = _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10;

  /* Linear axis */
  unsigned int i = 0;
  if (axis+i == 0 && jachq.size(0) > i)
  {
    /* axis 0 = Linear position */
    jachq.setValue(i, 0, -x1 + x3 - x5 - x7);
    jachq.setValue(i, 1, -x10 - x11 + x8 - x9);
    jachq.setValue(i, 2, x12 - x13 - x14 - x15);
    jachq.setValue(i, 3, -2*x16*x17 - 2*x18*x19 - 2*x2*x29 - 2*x20*x21 - 2*x27*x4 - 2*x28*x6);
    jachq.setValue(i, 4, -2*x0*x27 - 2*x17*x22 + 2*x18*x21 - 2*x19*x20 - 2*x2*x28 + 2*x29*x6);
    jachq.setValue(i, 5, -2*x0*x29 + 2*x16*x19 - 2*x17*x18 - 2*x21*x22 - 2*x27*x6 + 2*x28*x4);
    jachq.setValue(i, 6, -2*x0*x28 - 2*x16*x21 + 2*x17*x20 - 2*x19*x22 + 2*x2*x27 - 2*x29*x4);

    if (q0->numberOfBlocks()>=2)
    {
      jachq.setValue(i, 7, x17);
      jachq.setValue(i, 8, x21);
      jachq.setValue(i, 9, x19);
      jachq.setValue(i, 10, 2*x17*x25 + 2*x19*x26 + 2*x21*x24);
      jachq.setValue(i, 11, 2*x17*x23 + 2*x19*x24 - 2*x21*x26);
      jachq.setValue(i, 12, 2*x17*x26 - 2*x19*x25 + 2*x21*x23);
      jachq.setValue(i, 13, -2*x17*x24 + 2*x19*x23 + 2*x21*x25);
    }

    i++;
  }

  /* Rotational axis */
  if (axis+i == 1 && jachq.size(0) > i)
  {

    jachq.setValue(i, 0, 0);
    jachq.setValue(i, 1, 0);
    jachq.setValue(i, 2, 0);
    jachq.setValue(i, 3, 2*(-x51*x67 + x72*(-x17*x49 - x19*x53 - x21*x47 + x51*x66 + x55*x73 + x60*x74 + x65*x75))/(pow(x46, 2) + pow(x67, 2)));
    jachq.setValue(i, 4, x77*(x49*x76 + x78*(x17*x51 + x19*x47 - x21*x53 + x49*x66 - x55*x79 - x60*x75 + x65*x74)));
    jachq.setValue(i, 5, x77*(x47*x76 + x78*(x17*x53 - x19*x49 + x21*x51 + x47*x66 - x55*x74 + x60*x73 - x65*x79)));
    jachq.setValue(i, 6, x77*(x53*x76 + x78*(-x17*x47 + x19*x51 + x21*x49 + x53*x66 + x55*x75 - x60*x79 - x65*x73)));

    if (q0->numberOfBlocks()>=2)
    {
      /*
       * sympy expression:
       *
       * for i in range(4): print('jachq.setValue(i, {}, {});'.format(i+10,e[i+4]))
       */

      jachq.setValue(i, 7, 0);
      jachq.setValue(i, 8, 0);
      jachq.setValue(i, 9, 0);
      jachq.setValue(i, 10, x77*(x76*x80 - x78*(x17*x81 + x19*x83 + x21*x82 - x66*x80)));
      jachq.setValue(i, 11, x77*(x76*x81 + x78*(x17*x80 - x19*x82 + x21*x83 + x66*x81)));
      jachq.setValue(i, 12, x77*(x76*x82 + x78*(-x17*x83 + x19*x81 + x21*x80 + x66*x82)));
      jachq.setValue(i, 13, x77*(x76*x83 + x78*(x17*x82 + x19*x80 - x21*x81 + x66*x83)));
    }

    i++;
  }
}
