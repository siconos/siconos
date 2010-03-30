/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
/*! \file NewtonEulerR.h

*/

#include "PrismaticJoint.hpp"
#include <boost/math/quaternion.hpp>


int PrismaticJointR::_sNbEqualities = 5;


/*axe is the axis of the prismatic joint, in the frame of the first DS, d1.*/
PrismaticJointR::PrismaticJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SimpleVector axe): NewtonEulerR()
{
  _axe0 = axe;
  _d1 = d1;
  _d2 = d2;
  computeFromInitialPosition();
}
/*axe is the axis of the prismatic joint, in the absolute frame.*/
PrismaticJointR::PrismaticJointR(SP::NewtonEulerDS d2, SP::SimpleVector axe): NewtonEulerR()
{
  //    _d1=NULL;
  _axe0 = axe;
  _d2 = d2;
  computeFromInitialPosition();
}
void PrismaticJointR::displayInitialPosition()
{
  std::cout << "Prismatic axis :\n";
  _axe0->display();
  std::cout << "V1 :" << _V1x << " " << _V1y << " " << _V1z << "\n";
  std::cout << "V2 :" << _V2x << " " << _V2y << " " << _V2z << "\n";
  std::cout << "G10G20d1 :" << _G10G20d1x << " " << _G10G20d1y << " " << _G10G20d1z << "\n";
  std::cout << "q1cq2 :" << _q1cq202 << " " << _q1cq203 << " " << _q1cq204 << "\n";

}
void PrismaticJointR::computeFromInitialPosition()
{
  computeV1V2FromAxe();
  SP::SiconosVector q1;
  SP::SiconosVector q2;
  q2 = _d2->q0();

  if (_d1)
  {
    q1 = _d1->q0();
  }
  else
  {
    q1.reset(new SimpleVector(7));
    q1->zero();
    q1->setValue(3, 1);
  }
  ::boost::math::quaternion<float>    quat1(q1->getValue(3), q1->getValue(4), q1->getValue(5), q1->getValue(6));
  ::boost::math::quaternion<float>    quat2(q2->getValue(3), q2->getValue(4), q2->getValue(5), q2->getValue(6));
  ::boost::math::quaternion<float>    quat1_inv(q1->getValue(3), -q1->getValue(4), -q1->getValue(5), -q1->getValue(6));
  ::boost::math::quaternion<float>    quatG10G20_abs(0, q2->getValue(0) - q1->getValue(0), q2->getValue(1) - q1->getValue(1), q2->getValue(2) - q1->getValue(2));
  ::boost::math::quaternion<float>    quatBuff(0, 0, 0, 0);
  quatBuff = quat1_inv * quatG10G20_abs * quat1;
  _G10G20d1x = quatBuff.R_component_2();
  _G10G20d1y = quatBuff.R_component_3();
  _G10G20d1z = quatBuff.R_component_4();
  quatBuff = quat1 / quat2;
  _q1cq202 = quatBuff.R_component_2();
  _q1cq203 = quatBuff.R_component_3();
  _q1cq204 = quatBuff.R_component_4();
  //  displayInitialPosition();
}

void PrismaticJointR::computeV1V2FromAxe()
{
  _V1.reset(new SimpleVector(3));
  _V2.reset(new SimpleVector(3));
  _V1->zero();
  _V2->zero();
  //build _V1
  if (_axe0->getValue(0) > _axe0->getValue(1))
    if (_axe0->getValue(0) > _axe0->getValue(2))
    {
      _V1->setValue(1, -_axe0->getValue(0));
      _V1->setValue(0, _axe0->getValue(1));
    }
    else
    {
      _V1->setValue(1, -_axe0->getValue(2));
      _V1->setValue(2, _axe0->getValue(1));
    }
  else if (_axe0->getValue(2) > _axe0->getValue(1))
  {
    _V1->setValue(1, -_axe0->getValue(2));
    _V1->setValue(2, _axe0->getValue(1));
  }
  else
  {
    _V1->setValue(1, -_axe0->getValue(0));
    _V1->setValue(0, _axe0->getValue(1));
  }
  double aux = 1 / _V1->norm2();
  scal(aux, *_V1, *_V1);
  cross_product(*_axe0, *_V1, *_V2);
  _V1x = _V1->getValue(0);
  _V1y = _V1->getValue(1);
  _V1z = _V1->getValue(2);
  _V2x = _V2->getValue(0);
  _V2y = _V2->getValue(1);
  _V2z = _V2->getValue(2);

}



void PrismaticJointR::computeJachq(double t)
{
  _jachq->zero();
  SP::SiconosVector x2 = _d2->q();
  double X2 = x2->getValue(0);
  double Y2 = x2->getValue(1);
  double Z2 = x2->getValue(2);
  double q20 = x2->getValue(3);
  double q21 = x2->getValue(4);
  double q22 = x2->getValue(5);
  double q23 = x2->getValue(6);

  double X1 = 0;
  double Y1 = 0;
  double Z1 = 0;
  double q10 = 1;
  double q11 = 0;
  double q12 = 0;
  double q13 = 0;
  if (_d1)
  {
    SP::SiconosVector x1 = _d1->q();
    X1 = x1->getValue(0);
    Y1 = x1->getValue(1);
    Z1 = x1->getValue(2);
    q10 = x1->getValue(3);
    q11 = x1->getValue(4);
    q12 = x1->getValue(5);
    q13 = x1->getValue(6);
    Jd1d2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);
  }
  else
    Jd2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);



}



void PrismaticJointR::computeh(double t)
{

  SP::SiconosVector x2 = _d2->q();
  double X2 = x2->getValue(0);
  double Y2 = x2->getValue(1);
  double Z2 = x2->getValue(2);
  double q20 = x2->getValue(3);
  double q21 = x2->getValue(4);
  double q22 = x2->getValue(5);
  double q23 = x2->getValue(6);
  double X1 = 0;
  double Y1 = 0;
  double Z1 = 0;
  double q10 = 1;
  double q11 = 0;
  double q12 = 0;
  double q13 = 0;
  if (_d1)
  {
    SP::SiconosVector x1 = _d1->q();
    X1 = x1->getValue(0);
    Y1 = x1->getValue(1);
    Z1 = x1->getValue(2);
    q10 = x1->getValue(3);
    q11 = x1->getValue(4);
    q12 = x1->getValue(5);
    q13 = x1->getValue(6);
  }
  SP::SiconosVector y = interaction()->y(0);
  y->setValue(0, H1(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y->setValue(1, H2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y->setValue(2, H3(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y->setValue(3, H4(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y->setValue(4, H5(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  double norm = 0;
  for (int ii = 0; ii < 5; ii++)
    norm += y->getValue(ii) * y->getValue(ii);
  std::cout << "Prismatic norm computeH: " << norm << std::endl;



}



/* The options were    : operatorarrow */
double PrismaticJointR::H1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  double t1;
  double t12;
  double t14;
  double t15;
  double t18;
  double t2;
  double t27;
  double t28;
  double t3;
  double t4;
  double t6;
  double t8;
  double t9;
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
    t27 = q12 * q13;
    t28 = q10 * q11;
    return(((t1 + t2 - t3 - t4) * t6 + (2.0 * t8 + 2.0 * t9) * t12 + (2.0 * t14 - 2.0 * t15) * t18) * _V1x + ((
             2.0 * t8 - 2.0 * t9) * t6 + (t4 - t3 + t2 - t1) * t12 + (2.0 * t27 + 2.0 * t28) * t18) * _V1y + ((2.0 * t14 + 2.0 *
                 t15) * t6 + (2.0 * t27 - 2.0 * t28) * t12 + (t3 - t4 - t1 + t2) * t18) * _V1z - _G10G20d1x * _V1x - _G10G20d1y *
           _V1y - _G10G20d1z * _V1z);
  }
}

/* The options were    : operatorarrow */
double PrismaticJointR::H2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  double t1;
  double t12;
  double t14;
  double t15;
  double t18;
  double t2;
  double t27;
  double t28;
  double t3;
  double t4;
  double t6;
  double t8;
  double t9;
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
    t27 = q12 * q13;
    t28 = q10 * q11;
    return(((t1 + t2 - t3 - t4) * t6 + (2.0 * t8 + 2.0 * t9) * t12 + (2.0 * t14 - 2.0 * t15) * t18) * _V2x + ((
             2.0 * t8 - 2.0 * t9) * t6 + (t4 - t3 + t2 - t1) * t12 + (2.0 * t27 + 2.0 * t28) * t18) * _V2y + ((2.0 * t14 + 2.0 *
                 t15) * t6 + (2.0 * t27 - 2.0 * t28) * t12 + (t3 - t4 - t1 + t2) * t18) * _V2z - _G10G20d1x * _V2x - _G10G20d1y *
           _V2y - _G10G20d1z * _V2z);
  }
}

/* The options were    : operatorarrow */
double PrismaticJointR::H3(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  {
    return(q11 * q20 - q10 * q21 + q13 * q22 - q12 * q23 - _q1cq202);
  }
}

/* The options were    : operatorarrow */
double PrismaticJointR::H5(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  {
    return(q13 * q20 + q12 * q21 - q11 * q22 - q10 * q23 - _q1cq204);
  }
}

/* The options were    : operatorarrow */
double PrismaticJointR::H4(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  {
    return(q12 * q20 - q13 * q21 - q10 * q22 + q11 * q23 - _q1cq203);
  }
}

void PrismaticJointR::Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
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
    _jachq->setValue(0, 0, -t89);
    t90 = df[7];
    _jachq->setValue(0, 1, -t90);
    t91 = df[10];
    _jachq->setValue(0, 2, -t91);
    t92 = df[13];
    t94 = df[9];
    t96 = df[6];
    _jachq->setValue(0, 3, t92 * q11 + t94 * q12 + t96 * q13 + 2.0 * df[1]*q10);
    t102 = df[8];
    t104 = df[5];
    _jachq->setValue(0, 4, t92 * q10 + t102 * q13 + t104 * q12 + 2.0 * df[0]*q11);
    t109 = df[12];
    _jachq->setValue(0, 5, t109 * q13 + t94 * q10 + t104 * q11 + 2.0 * df[3]*q12);
    _jachq->setValue(0, 6, t109 * q12 + t102 * q11 + t96 * q10 + 2.0 * df[2]*q13);
    _jachq->setValue(0, 7, t89);
    _jachq->setValue(0, 8, t90);
    _jachq->setValue(0, 9, t91);
    _jachq->setValue(0, 10, 0.0);
    _jachq->setValue(0, 11, 0.0);
    _jachq->setValue(0, 12, 0.0);
    _jachq->setValue(0, 13, 0.0);
    t122 = dfr0[4];
    _jachq->setValue(1, 0, -t122);
    t123 = dfr0[7];
    _jachq->setValue(1, 1, -t123);
    t124 = dfr0[10];
    _jachq->setValue(1, 2, -t124);
    t125 = dfr0[13];
    t127 = dfr0[9];
    t129 = dfr0[6];
    _jachq->setValue(1, 3, t125 * q11 + t127 * q12 + t129 * q13 + 2.0 * dfr0[1]*q10);
    t135 = dfr0[8];
    t137 = dfr0[5];
    _jachq->setValue(1, 4, t125 * q10 + t135 * q13 + t137 * q12 + 2.0 * dfr0[0]*q11);
    t142 = dfr0[12];
    _jachq->setValue(1, 5, t142 * q13 + t127 * q10 + t137 * q11 + 2.0 * dfr0[3]*q12);
    _jachq->setValue(1, 6, t142 * q12 + t135 * q11 + t129 * q10 + 2.0 * dfr0[2]*q13);
    _jachq->setValue(1, 7, t122);
    _jachq->setValue(1, 8, t123);
    _jachq->setValue(1, 9, t124);
    _jachq->setValue(1, 10, 0.0);
    _jachq->setValue(1, 11, 0.0);
    _jachq->setValue(1, 12, 0.0);
    _jachq->setValue(1, 13, 0.0);
    _jachq->setValue(2, 0, 0.0);
    _jachq->setValue(2, 1, 0.0);
    _jachq->setValue(2, 2, 0.0);
    _jachq->setValue(2, 3, -q21);
    _jachq->setValue(2, 4, q20);
    _jachq->setValue(2, 5, -q23);
    _jachq->setValue(2, 6, q22);
    _jachq->setValue(2, 7, 0.0);
    _jachq->setValue(2, 8, 0.0);
    _jachq->setValue(2, 9, 0.0);
    _jachq->setValue(2, 10, q11);
    _jachq->setValue(2, 11, -q10);
    _jachq->setValue(2, 12, q13);
    _jachq->setValue(2, 13, -q12);
    _jachq->setValue(3, 0, 0.0);
    _jachq->setValue(3, 1, 0.0);
    _jachq->setValue(3, 2, 0.0);
    _jachq->setValue(3, 3, -q22);
    _jachq->setValue(3, 4, q23);
    _jachq->setValue(3, 5, q20);
    _jachq->setValue(3, 6, _jachq->getValue(2, 3));
    _jachq->setValue(3, 7, 0.0);
    _jachq->setValue(3, 8, 0.0);
    _jachq->setValue(3, 9, 0.0);
    _jachq->setValue(3, 10, q12);
    _jachq->setValue(3, 11, -q13);
    _jachq->setValue(3, 12, _jachq->getValue(2, 11));
    _jachq->setValue(3, 13, q11);
    _jachq->setValue(4, 0, 0.0);
    _jachq->setValue(4, 1, 0.0);
    _jachq->setValue(4, 2, 0.0);
    _jachq->setValue(4, 3, _jachq->getValue(2, 5));
    _jachq->setValue(4, 4, _jachq->getValue(3, 3));
    _jachq->setValue(4, 5, q21);
    _jachq->setValue(4, 6, q20);
    _jachq->setValue(4, 7, 0.0);
    _jachq->setValue(4, 8, 0.0);
    _jachq->setValue(4, 9, 0.0);
    _jachq->setValue(4, 10, q13);
    _jachq->setValue(4, 11, q12);
    _jachq->setValue(4, 12, -q11);
    _jachq->setValue(4, 13, _jachq->getValue(3, 12));
    return;
  }
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






