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
/*! \file NewtonEulerR.hpp

*/

#include "PrismaticJointR.hpp"
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

/**axe is the axis of the prismatic joint, in the frame of the first DS, d1.*/
PrismaticJointR::PrismaticJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SiconosVector axis): NewtonEulerR()
{
  _axis0 = axis;
  computeFromInitialPosition(d1->q(),d2->q());
}
/*axis is the axis of the prismatic joint, in the absolute frame.*/
PrismaticJointR::PrismaticJointR(SP::NewtonEulerDS d1, SP::SiconosVector axis): NewtonEulerR()
{
  //    _d1=NULL;
  _axis0 = axis;
  computeFromInitialPosition(d1->q());
}
void PrismaticJointR::displayInitialPosition()
{
  std::cout << "Prismatic axis :\n";
  _axis0->display();
  std::cout << "V1 :" << _V1x << " " << _V1y << " " << _V1z << "\n";
  std::cout << "V2 :" << _V2x << " " << _V2y << " " << _V2z << "\n";
  std::cout << "G10G20d1 :" << _G10G20d1x << " " << _G10G20d1y << " " << _G10G20d1z << "\n";
  std::cout << "cq2c10 :" << _cq2q101 << " " << _cq2q102
            << " " << _cq2q103 << " " << _cq2q104 << "\n";
}

void PrismaticJointR::computeFromInitialPosition(SP::SiconosVector q1, SP::SiconosVector q2)
{
  computeV1V2FromAxis();
  SP::SiconosVector q2i(new SiconosVector(7));
  q2i->zero();
  q2i->setValue(3, 1);

  if(q2)
    *q2i = *q2;

  ::boost::math::quaternion<double>    quat1(q1->getValue(3), q1->getValue(4), q1->getValue(5), q1->getValue(6));
  ::boost::math::quaternion<double>    quat2(q2i->getValue(3), q2i->getValue(4), q2i->getValue(5), q2i->getValue(6));
  ::boost::math::quaternion<double>    quat1_inv(q1->getValue(3), -q1->getValue(4), -q1->getValue(5), -q1->getValue(6));
  ::boost::math::quaternion<double>    quatG10G20_abs(0, q2i->getValue(0) - q1->getValue(0), q2i->getValue(1) - q1->getValue(1), q2i->getValue(2) - q1->getValue(2));
  ::boost::math::quaternion<double>    quatBuff(0, 0, 0, 0);
  quatBuff = quat1_inv * (quatG10G20_abs * quat1);
  _G10G20d1x = quatBuff.R_component_2();
  _G10G20d1y = quatBuff.R_component_3();
  _G10G20d1z = quatBuff.R_component_4();
  quatBuff = 1.0/quat2 * quat1;
  _cq2q101 = quatBuff.R_component_1();
  _cq2q102 = quatBuff.R_component_2();
  _cq2q103 = quatBuff.R_component_3();
  _cq2q104 = quatBuff.R_component_4();
  //  displayInitialPosition();
}

void PrismaticJointR::computeV1V2FromAxis()
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

void PrismaticJointR::computeJachq(double time, Interaction& inter,  SP::BlockVector q0)
{

  DEBUG_PRINT("PrismaticJointR::computeJachq(double time, Interaction& inter, SP::BlockVector q0 ) \n");

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

  if(q0->getNumberOfBlocks()>1)
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
    Jd2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);

  DEBUG_END("PrismaticJointR::computeJachq(double time, Interaction& inter, SP::BlockVector q0 ) \n");
}

void PrismaticJointR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  DEBUG_PRINT("PrismaticJointR::computeh(double time, BlockVector& q0, SiconosVector& y) \n");
  SP::SiconosVector q1 = (q0.getAllVect())[0];
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

  if (q0.getNumberOfBlocks()>1)
  {
    SP::SiconosVector q2 = (q0.getAllVect())[1];
    X2 = q0.getValue(7);
    Y2 = q0.getValue(8);
    Z2 = q0.getValue(9);
    q20 = q0.getValue(10);
    q21 = q0.getValue(11);
    q22 = q0.getValue(12);
    q23 = q0.getValue(13);
  }
  y.setValue(0, H1(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(1, H2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(2, H3(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(3, H4(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(4, H5(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));

  DEBUG_EXPR(y.display());
  DEBUG_PRINTF(" y.normInf() = %12.8e \n", y.normInf());

  // double norm = 0;
  // for(int ii = 0; ii < 5; ii++)
  //   norm += y.getValue(ii) * y.getValue(ii);
  //std::cout<<"Prismatic norm computeH: "<<norm<<std::endl;



}

// The rest of the code is generated.
// we can disable some warning
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

/* sympy expreesion:
 *
 * G1G2d1 = unrot(G2-G1, q1)
 *
 * H  = lambda V: np.dot(V, G1G2d1)
 * H0 = lambda V: np.dot(V, G10G20d1)
 *
 * H1 = H(V1) - H0(V1)
 * H2 = H(V2) - H0(V2)
 */

/* The options were    : operatorarrow */
double PrismaticJointR::H1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
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

/* The options were    : operatorarrow */
double PrismaticJointR::H2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
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
 * H3 = q2to1[1]
 * H4 = q2to1[2]
 * H5 = q2to1[3]
 */

/* The options were    : operatorarrow */
double PrismaticJointR::H3(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  return q10*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
    + q11*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
    + q12*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
    - q13*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21);
}

/* The options were    : operatorarrow */
double PrismaticJointR::H4(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  return q10*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
    - q11*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
    + q12*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
    + q13*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22);
}

/* The options were    : operatorarrow */
double PrismaticJointR::H5(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  return q10*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
    + q11*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
    - q12*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
    + q13*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23);
}

void PrismaticJointR::Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  /*
   * sympy expression:
   *
   * H = [H1,H2,H3,H4,H5]
   * dq = list(G1[1:])+list(q1)+list(G2[1:])+list(q2)
   *
   * jachq = [[h.diff(d) for d in dq] for h in H]
   *
   */

  /* Prismatic constraints (H1, H2)
   */
  _jachq->setValue(0, 0, _V1x*(-pow(q10, 2) - pow(q11, 2) + pow(q12, 2) + pow(q13, 2))
                   + _V1y*(2*q10*q13 - 2*q11*q12)
                   + _V1z*(-2*q10*q12 - 2*q11*q13));
  _jachq->setValue(0, 1, _V1x*(-2*q10*q13 - 2*q11*q12)
                   + _V1y*(-pow(q10, 2) + pow(q11, 2) - pow(q12, 2) + pow(q13, 2))
                   + _V1z*(2*q10*q11 - 2*q12*q13));
  _jachq->setValue(0, 2, _V1x*(2*q10*q12 - 2*q11*q13)
                   + _V1y*(-2*q10*q11 - 2*q12*q13)
                   + _V1z*(-pow(q10, 2) + pow(q11, 2) + pow(q12, 2) - pow(q13, 2)));
  _jachq->setValue(0, 3, _V1x*(2*q10*(-X1 + X2) - 2*q12*(-Z1 + Z2) + 2*q13*(-Y1 + Y2))
                   + _V1y*(2*q10*(-Y1 + Y2) + 2*q11*(-Z1 + Z2) - 2*q13*(-X1 + X2))
                   + _V1z*(2*q10*(-Z1 + Z2) - 2*q11*(-Y1 + Y2) + 2*q12*(-X1 + X2)));
  _jachq->setValue(0, 4, _V1x*(q11*(-X1 + X2) - q11*(X1 - X2) + q12*(-Y1 + Y2)
                               - q12*(Y1 - Y2) + 2*q13*(-Z1 + Z2))
                   + _V1y*(2*q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q11*(Y1 - Y2)
                           + q12*(-X1 + X2) - q12*(X1 - X2))
                   + _V1z*(-q10*(-Y1 + Y2) + q10*(Y1 - Y2) - 2*q11*(-Z1 + Z2)
                           + q13*(-X1 + X2) - q13*(X1 - X2)));
  _jachq->setValue(0, 5, _V1x*(-q10*(-Z1 + Z2) + q10*(Z1 - Z2) + q11*(-Y1 + Y2)
                               - q11*(Y1 - Y2) - 2*q12*(-X1 + X2))
                   + _V1y*(2*q11*(-X1 + X2) + q12*(-Y1 + Y2) - q12*(Y1 - Y2)
                           + q13*(-Z1 + Z2) - q13*(Z1 - Z2))
                   + _V1z*(2*q10*(-X1 + X2) - q12*(-Z1 + Z2) + q12*(Z1 - Z2)
                           + q13*(-Y1 + Y2) - q13*(Y1 - Y2)));
  _jachq->setValue(0, 6, _V1x*(2*q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q11*(Z1 - Z2)
                               - q13*(-X1 + X2) + q13*(X1 - X2))
                   + _V1y*(-q10*(-X1 + X2) + q10*(X1 - X2) + q12*(-Z1 + Z2)
                           - q12*(Z1 - Z2) - 2*q13*(-Y1 + Y2))
                   + _V1z*(q11*(-X1 + X2) - q11*(X1 - X2) + 2*q12*(-Y1 + Y2)
                           + q13*(-Z1 + Z2) - q13*(Z1 - Z2)));
  _jachq->setValue(0, 7, _V1x*(pow(q10, 2) + pow(q11, 2) - pow(q12, 2) - pow(q13, 2))
                   + _V1y*(-2*q10*q13 + 2*q11*q12) + _V1z*(2*q10*q12 + 2*q11*q13));
  _jachq->setValue(0, 8, _V1x*(2*q10*q13 + 2*q11*q12)
                   + _V1y*(pow(q10, 2) - pow(q11, 2) + pow(q12, 2) - pow(q13, 2))
                   + _V1z*(-2*q10*q11 + 2*q12*q13));
  _jachq->setValue(0, 9, _V1x*(-2*q10*q12 + 2*q11*q13)
                   + _V1y*(2*q10*q11 + 2*q12*q13)
                   + _V1z*(pow(q10, 2) - pow(q11, 2) - pow(q12, 2) + pow(q13, 2)));
  _jachq->setValue(0, 10, 0);
  _jachq->setValue(0, 11, 0);
  _jachq->setValue(0, 12, 0);
  _jachq->setValue(0, 13, 0);
  _jachq->setValue(1, 0, _V2x*(-pow(q10, 2) - pow(q11, 2) + pow(q12, 2) + pow(q13, 2))
                   + _V2y*(2*q10*q13 - 2*q11*q12)
                   + _V2z*(-2*q10*q12 - 2*q11*q13));
  _jachq->setValue(1, 1, _V2x*(-2*q10*q13 - 2*q11*q12)
                   + _V2y*(-pow(q10, 2) + pow(q11, 2) - pow(q12, 2) + pow(q13, 2))
                   + _V2z*(2*q10*q11 - 2*q12*q13));
  _jachq->setValue(1, 2, _V2x*(2*q10*q12 - 2*q11*q13)
                   + _V2y*(-2*q10*q11 - 2*q12*q13)
                   + _V2z*(-pow(q10, 2) + pow(q11, 2) + pow(q12, 2) - pow(q13, 2)));
  _jachq->setValue(1, 3, _V2x*(2*q10*(-X1 + X2) - 2*q12*(-Z1 + Z2) + 2*q13*(-Y1 + Y2))
                   + _V2y*(2*q10*(-Y1 + Y2) + 2*q11*(-Z1 + Z2) - 2*q13*(-X1 + X2))
                   + _V2z*(2*q10*(-Z1 + Z2) - 2*q11*(-Y1 + Y2) + 2*q12*(-X1 + X2)));
  _jachq->setValue(1, 4, _V2x*(q11*(-X1 + X2) - q11*(X1 - X2) + q12*(-Y1 + Y2)
                               - q12*(Y1 - Y2) + 2*q13*(-Z1 + Z2))
                   + _V2y*(2*q10*(-Z1 + Z2) - q11*(-Y1 + Y2) + q11*(Y1 - Y2)
                           + q12*(-X1 + X2) - q12*(X1 - X2))
                   + _V2z*(-q10*(-Y1 + Y2) + q10*(Y1 - Y2) - 2*q11*(-Z1 + Z2)
                           + q13*(-X1 + X2) - q13*(X1 - X2)));
  _jachq->setValue(1, 5, _V2x*(-q10*(-Z1 + Z2) + q10*(Z1 - Z2) + q11*(-Y1 + Y2)
                               - q11*(Y1 - Y2) - 2*q12*(-X1 + X2))
                   + _V2y*(2*q11*(-X1 + X2) + q12*(-Y1 + Y2) - q12*(Y1 - Y2)
                           + q13*(-Z1 + Z2) - q13*(Z1 - Z2))
                   + _V2z*(2*q10*(-X1 + X2) - q12*(-Z1 + Z2) + q12*(Z1 - Z2)
                           + q13*(-Y1 + Y2) - q13*(Y1 - Y2)));
  _jachq->setValue(1, 6, _V2x*(2*q10*(-Y1 + Y2) + q11*(-Z1 + Z2) - q11*(Z1 - Z2)
                               - q13*(-X1 + X2) + q13*(X1 - X2))
                   + _V2y*(-q10*(-X1 + X2) + q10*(X1 - X2) + q12*(-Z1 + Z2)
                           - q12*(Z1 - Z2) - 2*q13*(-Y1 + Y2))
                   + _V2z*(q11*(-X1 + X2) - q11*(X1 - X2) + 2*q12*(-Y1 + Y2)
                           + q13*(-Z1 + Z2) - q13*(Z1 - Z2)));
  _jachq->setValue(1, 7, _V2x*(pow(q10, 2) + pow(q11, 2) - pow(q12, 2) - pow(q13, 2))
                   + _V2y*(-2*q10*q13 + 2*q11*q12) + _V2z*(2*q10*q12 + 2*q11*q13));
  _jachq->setValue(1, 8, _V2x*(2*q10*q13 + 2*q11*q12)
                   + _V2y*(pow(q10, 2) - pow(q11, 2) + pow(q12, 2) - pow(q13, 2))
                   + _V2z*(-2*q10*q11 + 2*q12*q13));
  _jachq->setValue(1, 9, _V2x*(-2*q10*q12 + 2*q11*q13)
                   + _V2y*(2*q10*q11 + 2*q12*q13)
                   + _V2z*(pow(q10, 2) - pow(q11, 2) - pow(q12, 2) + pow(q13, 2)));
  _jachq->setValue(1, 10, 0);
  _jachq->setValue(1, 11, 0);
  _jachq->setValue(1, 12, 0);
  _jachq->setValue(1, 13, 0);

  /* Orientation constraints (H3, H4, H5)
   */
  _jachq->setValue(2, 0, 0.0);
  _jachq->setValue(2, 1, 0.0);
  _jachq->setValue(2, 2, 0.0);
  _jachq->setValue(2, 3, -_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22);
  _jachq->setValue(2, 4, _cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23);
  _jachq->setValue(2, 5, -_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20);
  _jachq->setValue(2, 6, _cq2q101*q22 + _cq2q102*q23 + _cq2q103*q20 - _cq2q104*q21);
  _jachq->setValue(2, 7, 0.0);
  _jachq->setValue(2, 8, 0.0);
  _jachq->setValue(2, 9, 0.0);
  _jachq->setValue(2, 10, _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12);
  _jachq->setValue(2, 11, -_cq2q101*q10 - _cq2q102*q11 - _cq2q103*q12 - _cq2q104*q13);
  _jachq->setValue(2, 12, _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10);
  _jachq->setValue(2, 13, -_cq2q101*q12 + _cq2q102*q13 + _cq2q103*q10 - _cq2q104*q11);
  _jachq->setValue(3, 0, 0.0);
  _jachq->setValue(3, 1, 0.0);
  _jachq->setValue(3, 2, 0.0);
  _jachq->setValue(3, 3, -_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21);
  _jachq->setValue(3, 4, _cq2q101*q23 - _cq2q102*q22 + _cq2q103*q21 + _cq2q104*q20);
  _jachq->setValue(3, 5, _cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23);
  _jachq->setValue(3, 6, -_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22);
  _jachq->setValue(3, 7, 0.0);
  _jachq->setValue(3, 8, 0.0);
  _jachq->setValue(3, 9, 0.0);
  _jachq->setValue(3, 10, _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11);
  _jachq->setValue(3, 11, -_cq2q101*q13 - _cq2q102*q12 + _cq2q103*q11 + _cq2q104*q10);
  _jachq->setValue(3, 12, -_cq2q101*q10 - _cq2q102*q11 - _cq2q103*q12 - _cq2q104*q13);
  _jachq->setValue(3, 13, _cq2q101*q11 - _cq2q102*q10 + _cq2q103*q13 - _cq2q104*q12);
  _jachq->setValue(4, 0, 0.0);
  _jachq->setValue(4, 1, 0.0);
  _jachq->setValue(4, 2, 0.0);
  _jachq->setValue(4, 3, -_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20);
  _jachq->setValue(4, 4, -_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21);
  _jachq->setValue(4, 5, _cq2q101*q21 + _cq2q102*q20 - _cq2q103*q23 + _cq2q104*q22);
  _jachq->setValue(4, 6, _cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23);
  _jachq->setValue(4, 7, 0.0);
  _jachq->setValue(4, 8, 0.0);
  _jachq->setValue(4, 9, 0.0);
  _jachq->setValue(4, 10, _cq2q101*q13 + _cq2q102*q12 - _cq2q103*q11 - _cq2q104*q10);
  _jachq->setValue(4, 11, _cq2q101*q12 - _cq2q102*q13 - _cq2q103*q10 + _cq2q104*q11);
  _jachq->setValue(4, 12, -_cq2q101*q11 + _cq2q102*q10 - _cq2q103*q13 + _cq2q104*q12);
  _jachq->setValue(4, 13, -_cq2q101*q10 - _cq2q102*q11 - _cq2q103*q12 - _cq2q104*q13);
}

void PrismaticJointR::Jd2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  double df[16];
  double dfr0[16];
  double t1;
  double t102;
  double t104;
  double t109;
  double t11;
  double t12;
  double t122;
  double t123;
  double t124;
  double t125;
  double t127;
  double t129;
  double t135;
  double t137;
  double t14;
  double t142;
  double t15;
  double t17;
  double t18;
  double t2;
  double t22;
  double t24;
  double t27;
  double t28;
  double t29;
  double t3;
  double t33;
  double t36;
  double t38;
  double t4;
  double t40;
  double t41;
  double t48;
  double t49;
  double t5;
  double t55;
  double t56;
  double t6;
  double t62;
  double t63;
  double t64;
  double t65;
  double t66;
  double t72;
  double t73;
  double t79;
  double t8;
  double t80;
  double t86;
  double t87;
  double t88;
  double t89;
  double t9;
  double t90;
  double t91;
  double t92;
  double t94;
  double t96;
  {
    t1 = q11 * q11;
    t2 = q10 * q10;
    t3 = q13 * q13;
    t4 = q12 * q12;
    t6 = X2 - X1;
    t8 = q11 * q12;
    t9 = q10 * q13;
    t12 = Y2 - Y1;
    t14 = q11 * q13;
    t15 = q10 * q12;
    t18 = Z2 - Z1;
    t5 = t1 + t2 - t3 - t4;
    t11 = 2.0 * t8 + 2.0 * t9;
    t17 = 2.0 * t14 - 2.0 * t15;
    t27 = q12 * q13;
    t28 = q10 * q11;
    t22 = 2.0 * t8 - 2.0 * t9;
    t24 = t4 - t3 + t2 - t1;
    t29 = 2.0 * t27 + 2.0 * t28;
    t33 = 2.0 * t14 + 2.0 * t15;
    t36 = 2.0 * t27 - 2.0 * t28;
    t38 = t3 - t4 - t1 + t2;
    t40 = _V1z * t12;
    t41 = _V1y * t18;
    df[13] = -2.0 * t40 + 2.0 * t41;
    df[12] = 2.0 * t40 + 2.0 * t41;
    df[10] = _V1z * t38 + _V1y * t29 + _V1x * t17;
    t48 = _V1z * t6;
    t49 = _V1x * t18;
    df[9] = 2.0 * t48 - 2.0 * t49;
    df[8] = 2.0 * t48 + 2.0 * t49;
    df[7] = _V1z * t36 + _V1y * t24 + _V1x * t11;
    t55 = _V1y * t6;
    t56 = _V1x * t12;
    df[6] = -2.0 * t55 + 2.0 * t56;
    df[5] = 2.0 * t55 + 2.0 * t56;
    df[4] = _V1z * t33 + _V1y * t22 + _V1x * t5;
    t62 = _V1z * t18;
    t63 = _V1y * t12;
    t64 = _V1x * t6;
    df[3] = -t62 + t63 - t64;
    df[2] = t62 - t63 - t64;
    df[1] = t62 + t63 + t64;
    df[0] = -t62 - t63 + t64;
    t65 = _V2z * t12;
    t66 = _V2y * t18;
    dfr0[13] = -2.0 * t65 + 2.0 * t66;
    dfr0[12] = 2.0 * t65 + 2.0 * t66;
    dfr0[10] = _V2z * t38 + _V2y * t29 + _V2x * t17;
    t72 = _V2z * t6;
    t73 = _V2x * t18;
    dfr0[9] = 2.0 * t72 - 2.0 * t73;
    dfr0[8] = 2.0 * t72 + 2.0 * t73;
    dfr0[7] = _V2z * t36 + _V2y * t24 + _V2x * t11;
    t79 = _V2y * t6;
    t80 = _V2x * t12;
    dfr0[6] = -2.0 * t79 + 2.0 * t80;
    dfr0[5] = 2.0 * t79 + 2.0 * t80;
    dfr0[4] = _V2z * t33 + _V2y * t22 + _V2x * t5;
    t86 = _V2z * t18;
    t87 = _V2y * t12;
    t88 = _V2x * t6;
    dfr0[3] = -t86 + t87 - t88;
    dfr0[2] = t86 - t87 - t88;
    dfr0[1] = t86 + t87 + t88;
    dfr0[0] = -t86 - t87 + t88;
    t89 = df[4];
    //    _jachq->setValue(0,0, -t89);
    t90 = df[7];
    //    _jachq->setValue(0,1, -t90);
    t91 = df[10];
    //    _jachq->setValue(0,2, -t91);
    t92 = df[13];
    t94 = df[9];
    t96 = df[6];
    //    _jachq->setValue(0,3, t92*q11+t94*q12+t96*q13+2.0*df[1]*q10);
    t102 = df[8];
    t104 = df[5];
    //    _jachq->setValue(0,4, t92*q10+t102*q13+t104*q12+2.0*df[0]*q11);
    t109 = df[12];
    //    _jachq->setValue(0,5, t109*q13+t94*q10+t104*q11+2.0*df[3]*q12);
    //    _jachq->setValue(0,6, t109*q12+t102*q11+t96*q10+2.0*df[2]*q13);
    _jachq->setValue(0, 7 - 7, t89);
    _jachq->setValue(0, 8 - 7, t90);
    _jachq->setValue(0, 9 - 7, t91);
    _jachq->setValue(0, 10 - 7, 0.0);
    _jachq->setValue(0, 11 - 7, 0.0);
    _jachq->setValue(0, 12 - 7, 0.0);
    _jachq->setValue(0, 13 - 7, 0.0);
    t122 = dfr0[4];
    //    _jachq->setValue(1,0, -t122);
    t123 = dfr0[7];
    //    _jachq->setValue(1,1, -t123);
    t124 = dfr0[10];
    //    _jachq->setValue(1,2, -t124);
    t125 = dfr0[13];
    t127 = dfr0[9];
    t129 = dfr0[6];
    //    _jachq->setValue(1,3, t125*q11+t127*q12+t129*q13+2.0*dfr0[1]*q10);
    t135 = dfr0[8];
    t137 = dfr0[5];
    //    _jachq->setValue(1,4, t125*q10+t135*q13+t137*q12+2.0*dfr0[0]*q11);
    t142 = dfr0[12];
    //    _jachq->setValue(1,5, t142*q13+t127*q10+t137*q11+2.0*dfr0[3]*q12);
    //    _jachq->setValue(1,6, t142*q12+t135*q11+t129*q10+2.0*dfr0[2]*q13);
    _jachq->setValue(1, 7 - 7, t122);
    _jachq->setValue(1, 8 - 7, t123);
    _jachq->setValue(1, 9 - 7, t124);
    _jachq->setValue(1, 10 - 7, 0.0);
    _jachq->setValue(1, 11 - 7, 0.0);
    _jachq->setValue(1, 12 - 7, 0.0);
    _jachq->setValue(1, 13 - 7, 0.0);
    //    _jachq->setValue(2,0, 0.0);
    //    _jachq->setValue(2,1, 0.0);
    //    _jachq->setValue(2,2, 0.0);
    //    _jachq->setValue(2,3, -q21);
    //    _jachq->setValue(2,4, q20);
    //    _jachq->setValue(2,5, -q23);
    //    _jachq->setValue(2,6, q22);
    _jachq->setValue(2, 7 - 7, 0.0);
    _jachq->setValue(2, 8 - 7, 0.0);
    _jachq->setValue(2, 9 - 7, 0.0);
    _jachq->setValue(2, 10 - 7, q11);
    _jachq->setValue(2, 11 - 7, -q10);
    _jachq->setValue(2, 12 - 7, q13);
    _jachq->setValue(2, 13 - 7, -q12);
    //    _jachq->setValue(3,0, 0.0);
    //    _jachq->setValue(3,1, 0.0);
    //    _jachq->setValue(3,2, 0.0);
    //    _jachq->setValue(3,3, -q22);
    //    _jachq->setValue(3,4, q23);
    //    _jachq->setValue(3,5, q20);
    //    _jachq->setValue(3,6, _jachq->getValue(2,3));
    _jachq->setValue(3, 7 - 7, 0.0);
    _jachq->setValue(3, 8 - 7, 0.0);
    _jachq->setValue(3, 9 - 7, 0.0);
    _jachq->setValue(3, 10 - 7, q12);
    _jachq->setValue(3, 11 - 7, -q13);
    _jachq->setValue(3, 12 - 7, _jachq->getValue(2, 11 - 7));
    _jachq->setValue(3, 13 - 7, q11);
    //    _jachq->setValue(4,0, 0.0);
    //    _jachq->setValue(4,1, 0.0);
    //    _jachq->setValue(4,2, 0.0);
    //    _jachq->setValue(4,3, _jachq->getValue(2,5));
    //    _jachq->setValue(4,4, _jachq->getValue(3,3));
    //    _jachq->setValue(4,5, q21);
    //    _jachq->setValue(4,6, q20);
    _jachq->setValue(4, 7 - 7, 0.0);
    _jachq->setValue(4, 8 - 7, 0.0);
    _jachq->setValue(4, 9 - 7, 0.0);
    _jachq->setValue(4, 10 - 7, q13);
    _jachq->setValue(4, 11 - 7, q12);
    _jachq->setValue(4, 12 - 7, -q11);
    _jachq->setValue(4, 13 - 7, _jachq->getValue(3, 12 - 7));
    return;
  }
}




void PrismaticJointR::computeDotJachq(double time, BlockVector& workQ, BlockVector& workZ, BlockVector& workQdot)
{
  std::cout << "Warning:  PrismaticJointR::computeDotJachq(...) not yet implemented" << std::endl;
}

void PrismaticJointR::DotJd1d2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13,
                               double Xdot2, double Ydot2, double Zdot2, double qdot20, double qdot21, double qdot22, double qdot23)
{
}

void PrismaticJointR::DotJd2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13,
                             double X2, double Y2, double Z2, double qdot20, double qdot21, double qdot22, double qdot23)
{
}
