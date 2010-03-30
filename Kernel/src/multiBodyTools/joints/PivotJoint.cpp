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
/*! \file PivotJoint.cpp

*/

#include "PivotJoint.hpp"
#include <boost/math/quaternion.hpp>

int PivotJointR::_sNbEqualities = 5;

PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SimpleVector P, SP::SimpleVector A): NewtonEulerR()
{
  _P0.reset(new SimpleVector(3));
  *_P0 = *P;
  _d1 = d1;
  SP::SiconosVector q1 = d1->q0();
  _G1P0x = _P0->getValue(0);
  _G1P0y = _P0->getValue(1);
  _G1P0z = _P0->getValue(2);
  _d2 = d2;
  SP::SiconosVector q2 = d2->q0();
  SimpleVector G2_abs(3);
  G2_abs.setValue(0, q2->getValue(0));
  G2_abs.setValue(1, q2->getValue(1));
  G2_abs.setValue(2, q2->getValue(2));
  ::boost::math::quaternion<float>    quat1(q1->getValue(3), q1->getValue(4), q1->getValue(5), q1->getValue(6));
  ::boost::math::quaternion<float>    quat2_inv(q2->getValue(3), -q2->getValue(4), -q2->getValue(5), -q2->getValue(6));
  ::boost::math::quaternion<float>    quatG1P0(0, _G1P0x, _G1P0y, _G1P0z);
  ::boost::math::quaternion<float>    quatA(0, A->getValue(0), A->getValue(1), A->getValue(2));
  ::boost::math::quaternion<float>    quatBuff(0, 0, 0, 0);


  quatBuff = quat1 * quatG1P0 / quat1;
  SimpleVector P0_abs(3);
  P0_abs.setValue(0, quatBuff.R_component_2() + q1->getValue(0));
  P0_abs.setValue(1, quatBuff.R_component_3() + q1->getValue(1));
  P0_abs.setValue(2, quatBuff.R_component_4() + q1->getValue(2));
  std::cout << "PivotJoint: P0_abs\n";
  P0_abs.display();
  SimpleVector G2P0_abs(3);
  G2P0_abs = P0_abs - G2_abs;
  ::boost::math::quaternion<float>    quatG2P0_abs(0, G2P0_abs.getValue(0), G2P0_abs.getValue(1), G2P0_abs.getValue(2));
  quatBuff = quat2_inv * quatG2P0_abs / quat2_inv;
  _G2P0x = quatBuff.R_component_2();
  _G2P0y = quatBuff.R_component_3();
  _G2P0z = quatBuff.R_component_4();
  std::cout << "PivotJoint G1P0 :" << _G1P0x << " " << _G1P0y << " " << _G1P0z << std::endl;
  std::cout << "PivotJoint G2P0 :" << _G2P0x << " " << _G2P0y << " " << _G2P0z << std::endl;


  /*calcul of axis _A*/
  quatBuff = quat1 * quatA / quat1;
  _Ax = quatBuff.R_component_2();
  _Ay = quatBuff.R_component_3();
  _Az = quatBuff.R_component_4();

  buildA1A2();
}
/* constructor,
   \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
   \param a SP::SimpleVector P0, P0 contains the coordinates of the Pivot point, in the absolute frame.
*/
PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::SimpleVector P0, SP::SimpleVector A): NewtonEulerR()
{
  _P0.reset(new SimpleVector(3));
  *_P0 = *P0;
  _d1 = d1;
  SP::SiconosVector q1 = d1->q0();


  ::boost::math::quaternion<float>    quat1(q1->getValue(3), q1->getValue(4), q1->getValue(5), q1->getValue(6));
  ::boost::math::quaternion<float>    quat1_inv(q1->getValue(3), -q1->getValue(4), -q1->getValue(5), -q1->getValue(6));
  ::boost::math::quaternion<float>    quatG1P0_abs_init_position(0, _P0->getValue(0) - q1->getValue(0), _P0->getValue(1) - q1->getValue(1), _P0->getValue(2) - q1->getValue(2));
  ::boost::math::quaternion<float>    quatBuff(0, 0, 0, 0);
  /*quadBuff contains the vector G1P0 if the object has no orientation.*/
  quatBuff = quat1_inv * quatG1P0_abs_init_position * quat1;


  _G1P0x = quatBuff.R_component_2();
  _G1P0y = quatBuff.R_component_3();
  _G1P0z = quatBuff.R_component_4();
  _G2P0x = 0;
  _G2P0y = 0;
  _G2P0z = 0;
  std::cout << "PivotJoint G1P0 :" << _G1P0x << " " << _G1P0y << " " << _G1P0z << std::endl;
  _Ax = A->getValue(0);
  _Ay = A->getValue(1);
  _Az = A->getValue(2);
  buildA1A2();
}
void PivotJointR::buildA1A2()
{
  double normA = sqrt(_Ax * _Ax + _Ay * _Ay + _Az * _Az);
  assert(normA > 0.9 && "PivotJoint, normA to small\n");
  _Ax /= normA;
  _Ay /= normA;
  _Az /= normA;
  std::cout << "JointKnee: _Ax _Ay _Az :" << _Ax << " " << _Ay << " " << _Az << std::endl;
  /*build _A1*/
  if (_Ax > _Ay)
  {
    if (_Ax > _Az) /*_Ax is the bigest*/
    {
      _A1x = _Ay;
      _A1y = -_Ax;
      _A1z = 0;
      _A2x = _Az;
      _A2z = -_Ax;
      _A2y = 0;
    }
    else  /*_Az is the biggest*/
    {
      _A1z = _Ax;
      _A1y = -_Az;
      _A1x = 0;
      _A2z = _Ay;
      _A2x = -_Az;
      _A2y = 0;
    }
  }
  else if (_Ay > _Az)  /*_Ay is the bigest*/
  {
    _A1y = _Ax;
    _A1x = -_Ay;
    _A1z = 0;
    _A2y = _Az;
    _A2z = -_Ay;
    _A2x = 0;
  }
  else  /*_Az is the biggest*/
  {
    _A1z = _Ax;
    _A1y = -_Az;
    _A1x = 0;
    _A2z = _Ay;
    _A2x = -_Az;
    _A2y = 0;
  }
  assert(fabs(_A1x * _Ax + _A1y * _Ay + _A1z * _Az) < 1e-9 && "PivotJoint, _A1 wrong\n");
  assert(fabs(_A2x * _Ax + _A2y * _Ay + _A2z * _Az) < 1e-9 && "PivotJoint, _A2 wrong\n");
  assert(fabs(_A1x * _A2x + _A1y * _A2y + _A1z * _A2z) < 1e-9 && "PivotJoint, _A12 wrong\n");
  std::cout << "JointKnee: _A1x _A1y _A1z :" << _A1x << " " << _A1y << " " << _A1z << std::endl;
  std::cout << "JointKnee: _A2x _A2y _A2z :" << _A2x << " " << _A2y << " " << _A2z << std::endl;
}
void PivotJointR::Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  double df[20];
  double dfr0[20];
  double dfr1[20];
  double t104;
  double t109;
  double t11;
  double t119;
  double t125;
  double t129;
  double t41;
  double t5;
  double t6;
  double t60;
  double t70;
  double t79;
  double t89;
  double t9;
  double t99;
  {
    t5 = 2.0 * _G2P0z;
    df[15] = -t5;
    df[14] = df[15];
    df[13] = 2.0 * _G2P0y;
    df[12] = -df[13];
    df[9] = -_G2P0x;
    df[8] = df[9];
    df[7] = 2.0 * _G1P0z;
    df[6] = df[7];
    t6 = 2.0 * _G1P0y;
    df[5] = -t6;
    df[4] = t6;
    df[3] = -_G1P0x;
    df[2] = df[3];
    dfr0[19] = t5;
    dfr0[18] = df[14];
    dfr0[17] = -df[6];
    dfr0[16] = df[6];
    t9 = 2.0 * _G2P0x;
    dfr0[13] = -t9;
    dfr0[12] = dfr0[13];
    dfr0[11] = -_G2P0y;
    dfr0[9] = dfr0[11];
    dfr0[5] = 2.0 * _G1P0x;
    dfr0[4] = dfr0[5];
    dfr0[2] = -_G1P0y;
    dfr0[0] = dfr0[2];
    dfr1[19] = df[12];
    dfr1[18] = dfr1[19];
    dfr1[17] = df[4];
    dfr1[16] = dfr1[17];
    dfr1[15] = t9;
    dfr1[14] = dfr0[12];
    dfr1[10] = -_G2P0z;
    dfr1[9] = dfr1[10];
    dfr1[7] = -dfr0[4];
    dfr1[6] = dfr0[4];
    dfr1[3] = -_G1P0z;
    dfr1[0] = dfr1[3];
    _jachq->setValue(0, 0,  1.0);
    _jachq->setValue(0, 1, 0.0);
    _jachq->setValue(0, 2, 0.0);
    t11 = df[5];
    _jachq->setValue(0, 3, df[6]*q12 + t11 * q13 + 2.0 * q10 * _G1P0x);
    _jachq->setValue(0, 4, dfr0[16]*q13 + dfr1[17]*q12 + 2.0 * q11 * _G1P0x);
    _jachq->setValue(0, 5, df[6]*q10 + dfr1[17]*q11 + 2.0 * df[2]*q12);
    _jachq->setValue(0, 6, dfr0[16]*q11 + t11 * q10 + 2.0 * df[2]*q13);
    _jachq->setValue(0, 7, -1.0);
    _jachq->setValue(0, 8, 0.0);
    _jachq->setValue(0, 9, 0.0);
    t41 = df[13];
    _jachq->setValue(0, 10, df[14]*q22 + t41 * q23 + 2.0 * df[8]*q20);
    _jachq->setValue(0, 11, dfr0[18]*q23 + dfr1[19]*q22 + 2.0 * df[8]*q21);
    _jachq->setValue(0, 12, df[14]*q20 + dfr1[19]*q21 + 2.0 * q22 * _G2P0x);
    _jachq->setValue(0, 13, dfr0[18]*q21 + t41 * q20 + 2.0 * q23 * _G2P0x);
    _jachq->setValue(1, 0, 0.0);
    _jachq->setValue(1, 1, 1.0);
    _jachq->setValue(1, 2, 0.0);
    t60 = dfr0[17];
    _jachq->setValue(1, 3, t60 * q11 + dfr0[4]*q13 + 2.0 * q10 * _G1P0y);
    _jachq->setValue(1, 4, t60 * q10 + dfr1[6]*q12 + 2.0 * dfr0[0]*q11);
    t70 = dfr0[16];
    _jachq->setValue(1, 5, t70 * q13 + dfr1[6]*q11 + 2.0 * q12 * _G1P0y);
    _jachq->setValue(1, 6, t70 * q12 + dfr0[4]*q10 + 2.0 * dfr0[0]*q13);
    _jachq->setValue(1, 7, 0.0);
    _jachq->setValue(1, 8, -1.0);
    _jachq->setValue(1, 9, 0.0);
    t79 = dfr0[19];
    _jachq->setValue(1, 10, t79 * q21 + dfr0[12]*q23 + 2.0 * dfr0[9]*q20);
    _jachq->setValue(1, 11, t79 * q20 + dfr1[14]*q22 + 2.0 * q21 * _G2P0y);
    t89 = dfr0[18];
    _jachq->setValue(1, 12, t89 * q23 + dfr1[14]*q21 + 2.0 * dfr0[9]*q22);
    _jachq->setValue(1, 13, t89 * q22 + dfr0[12]*q20 + 2.0 * q23 * _G2P0y);
    _jachq->setValue(2, 0, 0.0);
    _jachq->setValue(2, 1, 0.0);
    _jachq->setValue(2, 2, 1.0);
    t99 = dfr1[7];
    _jachq->setValue(2, 3, dfr1[16]*q11 + t99 * q12 + 2.0 * q10 * _G1P0z);
    t104 = dfr1[6];
    _jachq->setValue(2, 4, dfr1[16]*q10 + t104 * q13 + 2.0 * dfr1[0]*q11);
    t109 = dfr1[16];
    _jachq->setValue(2, 5, t109 * q13 + t99 * q10 + 2.0 * dfr1[0]*q12);
    _jachq->setValue(2, 6, t109 * q12 + t104 * q11 + 2.0 * q13 * _G1P0z);
    _jachq->setValue(2, 7, 0.0);
    _jachq->setValue(2, 8, 0.0);
    _jachq->setValue(2, 9, -1.0);
    t119 = dfr1[15];
    _jachq->setValue(2, 10, dfr1[18]*q21 + t119 * q22 + 2.0 * dfr1[9]*q20);
    t125 = dfr1[14];
    _jachq->setValue(2, 11, dfr1[18]*q20 + t125 * q23 + 2.0 * q21 * _G2P0z);
    t129 = dfr1[18];
    _jachq->setValue(2, 12, t129 * q23 + t119 * q20 + 2.0 * q22 * _G2P0z);
    _jachq->setValue(2, 13, t129 * q22 + t125 * q21 + 2.0 * dfr1[9]*q23);


    _jachq->setValue(3, 0, 0);
    _jachq->setValue(3, 1, 0);
    _jachq->setValue(3, 2, 0);
    _jachq->setValue(3, 3, _A1x * (-q21) + _A1y * (-q22) + _A1z * (-q23));
    _jachq->setValue(3, 4, _A1x * (q20) + _A1y * (q23) + _A1z * (-q22));
    _jachq->setValue(3, 5, _A1x * (-q23) + _A1y * (q20) + _A1z * (q21));
    _jachq->setValue(3, 6, _A1x * (q22) + _A1y * (-q21) + _A1z * (q20));
    _jachq->setValue(3, 7, 0);
    _jachq->setValue(3, 8, 0);
    _jachq->setValue(3, 9, 0);
    _jachq->setValue(3, 10, _A1x * (q11) + _A1y * (q12) + _A1z * (q13));
    _jachq->setValue(3, 11, _A1x * (-q10) + _A1y * (-q13) + _A1z * (q12));
    _jachq->setValue(3, 12, _A1x * (q13) + _A1y * (-q10) + _A1z * (-q11));
    _jachq->setValue(3, 13, _A1x * (-q12) + _A1y * (q11) + _A1z * (-q10));

    _jachq->setValue(4, 0, 0);
    _jachq->setValue(4, 1, 0);
    _jachq->setValue(4, 2, 0);
    _jachq->setValue(4, 3, _A2x * (-q21) + _A2y * (-q22) + _A2z * (-q23));
    _jachq->setValue(4, 4, _A2x * (q20) + _A2y * (q23) + _A2z * (-q22));
    _jachq->setValue(4, 5, _A2x * (-q23) + _A2y * (q20) + _A2z * (q21));
    _jachq->setValue(4, 6, _A2x * (q22) + _A2y * (-q21) + _A2z * (q20));
    _jachq->setValue(4, 7, 0);
    _jachq->setValue(4, 8, 0);
    _jachq->setValue(4, 9, 0);
    _jachq->setValue(4, 10, _A2x * (q11) + _A2y * (q12) + _A2z * (q13));
    _jachq->setValue(4, 11, _A2x * (-q10) + _A2y * (-q13) + _A2z * (q12));
    _jachq->setValue(4, 12, _A2x * (q13) + _A2y * (-q10) + _A2z * (-q11));
    _jachq->setValue(4, 13, _A2x * (-q12) + _A2y * (q11) + _A2z * (-q10));


  }
  //_jachq->display();
}

void PivotJointR::Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13)
{
  double df[20];
  double dfr0[20];
  double dfr1[20];
  double t104;
  double t109;
  double t11;
  double t119;
  double t125;
  double t129;
  double t41;
  double t5;
  double t6;
  double t60;
  double t70;
  double t79;
  double t89;
  double t9;
  double t99;

  double X2;
  double Y2;
  double Z2;
  double q20;
  double q21;
  double q22;
  double q23;
  X2 = 0;
  Y2 = 0;
  Z2 = 0;
  q20 = 0;
  q21 = 0;
  q22 = 0;
  q23 = 0;
  {
    t5 = 2.0 * _G2P0z;
    df[15] = -t5;
    df[14] = df[15];
    df[13] = 2.0 * _G2P0y;
    df[12] = -df[13];
    df[9] = -_G2P0x;
    df[8] = df[9];
    df[7] = 2.0 * _G1P0z;
    df[6] = df[7];
    t6 = 2.0 * _G1P0y;
    df[5] = -t6;
    df[4] = t6;
    df[3] = -_G1P0x;
    df[2] = df[3];
    dfr0[19] = t5;
    dfr0[18] = df[14];
    dfr0[17] = -df[6];
    dfr0[16] = df[6];
    t9 = 2.0 * _G2P0x;
    dfr0[13] = -t9;
    dfr0[12] = dfr0[13];
    dfr0[11] = -_G2P0y;
    dfr0[9] = dfr0[11];
    dfr0[5] = 2.0 * _G1P0x;
    dfr0[4] = dfr0[5];
    dfr0[2] = -_G1P0y;
    dfr0[0] = dfr0[2];
    dfr1[19] = df[12];
    dfr1[18] = dfr1[19];
    dfr1[17] = df[4];
    dfr1[16] = dfr1[17];
    dfr1[15] = t9;
    dfr1[14] = dfr0[12];
    dfr1[10] = -_G2P0z;
    dfr1[9] = dfr1[10];
    dfr1[7] = -dfr0[4];
    dfr1[6] = dfr0[4];
    dfr1[3] = -_G1P0z;
    dfr1[0] = dfr1[3];
    _jachq->setValue(0, 0,  1.0);
    _jachq->setValue(0, 1, 0.0);
    _jachq->setValue(0, 2, 0.0);
    t11 = df[5];
    _jachq->setValue(0, 3, df[6]*q12 + t11 * q13 + 2.0 * q10 * _G1P0x);
    _jachq->setValue(0, 4, dfr0[16]*q13 + dfr1[17]*q12 + 2.0 * q11 * _G1P0x);
    _jachq->setValue(0, 5, df[6]*q10 + dfr1[17]*q11 + 2.0 * df[2]*q12);
    _jachq->setValue(0, 6, dfr0[16]*q11 + t11 * q10 + 2.0 * df[2]*q13);
    //      _jachq->setValue(0,7, -1.0);
    //      _jachq->setValue(0,8, 0.0);
    //      _jachq->setValue(0,9, 0.0);
    t41 = df[13];
    //      _jachq->setValue(0,10, df[14]*q22+t41*q23+2.0*df[8]*q20);
    //      _jachq->setValue(0,11, dfr0[18]*q23+dfr1[19]*q22+2.0*df[8]*q21);
    //      _jachq->setValue(0,12, df[14]*q20+dfr1[19]*q21+2.0*q22*_G2P0x);
    //      _jachq->setValue(0,13, dfr0[18]*q21+t41*q20+2.0*q23*_G2P0x);
    _jachq->setValue(1, 0, 0.0);
    _jachq->setValue(1, 1, 1.0);
    _jachq->setValue(1, 2, 0.0);
    t60 = dfr0[17];
    _jachq->setValue(1, 3, t60 * q11 + dfr0[4]*q13 + 2.0 * q10 * _G1P0y);
    _jachq->setValue(1, 4, t60 * q10 + dfr1[6]*q12 + 2.0 * dfr0[0]*q11);
    t70 = dfr0[16];
    _jachq->setValue(1, 5, t70 * q13 + dfr1[6]*q11 + 2.0 * q12 * _G1P0y);
    _jachq->setValue(1, 6, t70 * q12 + dfr0[4]*q10 + 2.0 * dfr0[0]*q13);
    //      _jachq->setValue(1,7, 0.0);
    //      _jachq->setValue(1,8, -1.0);
    //      _jachq->setValue(1,9, 0.0);
    t79 = dfr0[19];
    //      _jachq->setValue(1,10, t79*q21+dfr0[12]*q23+2.0*dfr0[9]*q20);
    //      _jachq->setValue(1,11, t79*q20+dfr1[14]*q22+2.0*q21*_G2P0y);
    t89 = dfr0[18];
    //      _jachq->setValue(1,12, t89*q23+dfr1[14]*q21+2.0*dfr0[9]*q22);
    //      _jachq->setValue(1,13, t89*q22+dfr0[12]*q20+2.0*q23*_G2P0y);
    _jachq->setValue(2, 0, 0.0);
    _jachq->setValue(2, 1, 0.0);
    _jachq->setValue(2, 2, 1.0);
    t99 = dfr1[7];
    _jachq->setValue(2, 3, dfr1[16]*q11 + t99 * q12 + 2.0 * q10 * _G1P0z);
    t104 = dfr1[6];
    _jachq->setValue(2, 4, dfr1[16]*q10 + t104 * q13 + 2.0 * dfr1[0]*q11);
    t109 = dfr1[16];
    _jachq->setValue(2, 5, t109 * q13 + t99 * q10 + 2.0 * dfr1[0]*q12);
    _jachq->setValue(2, 6, t109 * q12 + t104 * q11 + 2.0 * q13 * _G1P0z);
    //      _jachq->setValue(2,7, 0.0);
    //      _jachq->setValue(2,8, 0.0);
    //      _jachq->setValue(2,9, -1.0);
    t119 = dfr1[15];
    //      _jachq->setValue(2,10, dfr1[18]*q21+t119*q22+2.0*dfr1[9]*q20);
    t125 = dfr1[14];
    //      _jachq->setValue(2,11, dfr1[18]*q20+t125*q23+2.0*q21*_G2P0z);
    t129 = dfr1[18];
    //      _jachq->setValue(2,12, t129*q23+t119*q20+2.0*q22*_G2P0z);
    //      _jachq->setValue(2,13, t129*q22+t125*q21+2.0*dfr1[9]*q23);
    _jachq->setValue(3, 0, 0);
    _jachq->setValue(3, 1, 0);
    _jachq->setValue(3, 2, 0);
    _jachq->setValue(3, 3, 0);
    _jachq->setValue(3, 4, _A1x);
    _jachq->setValue(3, 5, _A1y);
    _jachq->setValue(3, 6, _A1z);

    _jachq->setValue(4, 0, 0);
    _jachq->setValue(4, 1, 0);
    _jachq->setValue(4, 2, 0);
    _jachq->setValue(4, 3, 0);
    _jachq->setValue(4, 4, _A2x);
    _jachq->setValue(4, 5, _A2y);
    _jachq->setValue(4, 6, _A2z);


  }
}



void PivotJointR::computeJachq(double t)
{
  _jachq->zero();
  SP::SiconosVector x1 = _d1->q();
  double X1 = x1->getValue(0);
  double Y1 = x1->getValue(1);
  double Z1 = x1->getValue(2);
  double q10 = x1->getValue(3);
  double q11 = x1->getValue(4);
  double q12 = x1->getValue(5);
  double q13 = x1->getValue(6);

  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;
  if (_d2)
  {
    SP::SiconosVector x2 = _d2->q();
    X2 = x2->getValue(0);
    Y2 = x2->getValue(1);
    Z2 = x2->getValue(2);
    q20 = x2->getValue(3);
    q21 = x2->getValue(4);
    q22 = x2->getValue(5);
    q23 = x2->getValue(6);
    Jd1d2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);
  }
  else
    Jd1(X1, Y1, Z1, q10, q11, q12, q13);



}


/** to compute p
 *  \param double : current time
 *  \param unsigned int: "derivative" order of lambda used to compute input
 */
double PivotJointR::Hx(double X1, double Y1, double Z1, double  q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)

{
  double t1;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t3;
  double t4;
  {
    t1 = q11 * q11;
    t2 = q10 * q10;
    t3 = q13 * q13;
    t4 = q12 * q12;
    t17 = q21 * q21;
    t18 = q20 * q20;
    t19 = q23 * q23;
    t20 = q22 * q22;
    return(X1 - X2 + (t1 + t2 - t3 - t4) * _G1P0x + (2.0 * q11 * q12 - 2.0 * q10 * q13) * _G1P0y + (2.0 * q11 *
           q13 + 2.0 * q10 * q12) * _G1P0z - (t17 + t18 - t19 - t20) * _G2P0x - (2.0 * q21 * q22 - 2.0 * q20 * q23) * _G2P0y -
           (2.0 * q21 * q23 + 2.0 * q20 * q22) * _G2P0z);
  }
}

/* The options were    : operatorarrow */
double PivotJointR::Hy(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)

{
  double t22;
  double t23;
  double t24;
  double t25;
  double t6;
  double t7;
  double t8;
  double t9;
  {
    t6 = q12 * q12;
    t7 = q13 * q13;
    t8 = q10 * q10;
    t9 = q11 * q11;
    t22 = q22 * q22;
    t23 = q23 * q23;
    t24 = q20 * q20;
    t25 = q21 * q21;
    return(Y1 - Y2 + (2.0 * q11 * q12 + 2.0 * q10 * q13) * _G1P0x + (t6 - t7 + t8 - t9) * _G1P0y + (2.0 * q12 *
           q13 - 2.0 * q10 * q11) * _G1P0z - (2.0 * q21 * q22 + 2.0 * q20 * q23) * _G2P0x - (t22 - t23 + t24 - t25) * _G2P0y -
           (2.0 * q22 * q23 - 2.0 * q20 * q21) * _G2P0z);
  }
}

/* The options were    : operatorarrow */
double PivotJointR::Hz(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)

{
  double t11;
  double t12;
  double t13;
  double t14;
  double t27;
  double t28;
  double t29;
  double t30;
  {
    t11 = q13 * q13;
    t12 = q12 * q12;
    t13 = q11 * q11;
    t14 = q10 * q10;
    t27 = q23 * q23;
    t28 = q22 * q22;
    t29 = q21 * q21;
    t30 = q20 * q20;
    return(Z1 - Z2 + (2.0 * q11 * q13 - 2.0 * q10 * q12) * _G1P0x + (2.0 * q12 * q13 + 2.0 * q10 * q11) *
           _G1P0y + (t11 - t12 - t13 + t14) * _G1P0z - (2.0 * q21 * q23 - 2.0 * q20 * q22) * _G2P0x - (2.0 * q22 * q23 + 2.0 *
               q20 * q21) * _G2P0y - (t27 - t28 - t29 + t30) * _G2P0z);
  }
}
double PivotJointR::AscalA1(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23)
{
  ::boost::math::quaternion<float>    quat1(q10, q11, q12, q13);
  ::boost::math::quaternion<float>    quat2_inv(q20, -q21, -q22, -q23);
  ::boost::math::quaternion<float>    quatBuff = quat1 * quat2_inv;
  double aX = quatBuff.R_component_2();
  double aY = quatBuff.R_component_3();
  double aZ = quatBuff.R_component_4();
  return _A1x * aX + _A1y * aY + _A1z * aZ;
}
double PivotJointR::AscalA2(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23)
{
  ::boost::math::quaternion<float>    quat1(q10, q11, q12, q13);
  ::boost::math::quaternion<float>    quat2_inv(q20, -q21, -q22, -q23);
  ::boost::math::quaternion<float>    quatBuff = quat1 * quat2_inv;
  double aX = quatBuff.R_component_2();
  double aY = quatBuff.R_component_3();
  double aZ = quatBuff.R_component_4();
  return _A2x * aX + _A2y * aY + _A2z * aZ;
}
void PivotJointR::computeh(double t)
{

  SP::SiconosVector x1 = _d1->q();
  double X1 = x1->getValue(0);
  double Y1 = x1->getValue(1);
  double Z1 = x1->getValue(2);
  double q10 = x1->getValue(3);
  double q11 = x1->getValue(4);
  double q12 = x1->getValue(5);
  double q13 = x1->getValue(6);
  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;
  if (_d2)
  {
    SP::SiconosVector x2 = _d2->q();
    X2 = x2->getValue(0);
    Y2 = x2->getValue(1);
    Z2 = x2->getValue(2);
    q20 = x2->getValue(3);
    q21 = x2->getValue(4);
    q22 = x2->getValue(5);
    q23 = x2->getValue(6);
  }

  SP::SiconosVector y = interaction()->y(0);
  y->setValue(0, Hx(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y->setValue(1, Hy(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y->setValue(2, Hz(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y->setValue(3, AscalA1(q10, q11, q12, q13, q20, q21, q22, q23));
  y->setValue(4, AscalA2(q10, q11, q12, q13, q20, q21, q22, q23));

  std::cout << "PivotJoint computeH:\n";
  y->display();
}


