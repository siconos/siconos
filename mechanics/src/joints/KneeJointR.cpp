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
/*! \file KneeJointR.cpp

*/

#include "KneeJointR.hpp"
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

void KneeJointR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  NewtonEulerR::initComponents(inter, DSlink, workV, workM);
  if (!_dotjachq)
  {
    unsigned int sizeY = inter.getSizeOfY();
    unsigned int xSize = inter.getSizeOfDS();
    unsigned int qSize = 7 * (xSize / 6);

    _dotjachq.reset(new SimpleMatrix(sizeY, qSize));
  }
}


void KneeJointR::checkInitPos( SP::SiconosVector x1 ,  SP::SiconosVector x2 )
{
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
  if(x2)
  {
    X2 = x2->getValue(0);
    Y2 = x2->getValue(1);
    Z2 = x2->getValue(2);
    q20 = x2->getValue(3);
    q21 = x2->getValue(4);
    q22 = x2->getValue(5);
    q23 = x2->getValue(6);
  }

  assert(Hx(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23)
         < DBL_EPSILON);
  assert(Hy(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23)
         < DBL_EPSILON);
  assert(Hz(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23)
         < DBL_EPSILON);
}

KneeJointR::KneeJointR()
  : NewtonEulerJointR()
  , _P0(std11::make_shared<SiconosVector>(3))
{
  _points.resize(1);
}

KneeJointR::KneeJointR(SP::SiconosVector P, bool absoluteRef,
                       SP::NewtonEulerDS d1, SP::NewtonEulerDS d2)
  : NewtonEulerJointR()
  , _P0(std11::make_shared<SiconosVector>(3))
{
  _points.resize(1);
  setAbsolute(absoluteRef);
  setPoint(0, P);
  if (d1) {
    setInitialConditions(d1->q(), d2 ? d2->q() : SP::SiconosVector());
    checkInitPos(d1->q(), d2 ? d2->q() : SP::SiconosVector());
  }
}

static ::boost::math::quaternion<double> rotquat(const SP::SiconosVector& v)
{
  if (v)
    return ::boost::math::quaternion<double>((*v)(3),(*v)(4),(*v)(5),(*v)(6));
  else
    return ::boost::math::quaternion<double>(1, 0, 0, 0);
}

static ::boost::math::quaternion<double> rotquat(const SiconosVector& v)
{
  return ::boost::math::quaternion<double>(v(3),v(4),v(5),v(6));
}

static ::boost::math::quaternion<double> posquat(const SP::SiconosVector& v)
{
  return ::boost::math::quaternion<double>(0, (*v)(0),(*v)(1),(*v)(2));
}

static ::boost::math::quaternion<double> posquat(const SiconosVector& v)
{
  return ::boost::math::quaternion<double>(0, v(0),v(1),v(2));
}

void KneeJointR::setInitialConditions(SP::SiconosVector q1,
                                      SP::SiconosVector q2)
{
  *_P0 = *_points[0];
  boost::math::quaternion<double> rot1(rotquat(q1));
  boost::math::quaternion<double> quatBuff;

  if (q2)
  {
    boost::math::quaternion<double> rot2(rotquat(q2));

    if (_absoluteRef)
    {
      (*_P0)(0) -= (*q1)(0);
      (*_P0)(1) -= (*q1)(1);
      (*_P0)(2) -= (*q1)(2);
    }

    /** Computation of _G1P0 and _G2P0 */
    _G1P0x = _P0->getValue(0);
    _G1P0y = _P0->getValue(1);
    _G1P0z = _P0->getValue(2);

    boost::math::quaternion<double> quatG1P0(0, _G1P0x, _G1P0y, _G1P0z);
    quatBuff = rot1 * quatG1P0 / rot1;

    SiconosVector P0_abs(3);
    P0_abs.setValue(0, quatBuff.R_component_2() + q1->getValue(0));
    P0_abs.setValue(1, quatBuff.R_component_3() + q1->getValue(1));
    P0_abs.setValue(2, quatBuff.R_component_4() + q1->getValue(2));

    boost::math::quaternion<double> quatG2P0_abs(posquat(P0_abs) - posquat(q2));
    quatBuff = (1.0/rot2) * quatG2P0_abs * rot2;
    _G2P0x = quatBuff.R_component_2();
    _G2P0y = quatBuff.R_component_3();
    _G2P0z = quatBuff.R_component_4();
  }
  else
  {
    if (_absoluteRef)
    {
      /*quadG1P0_abs contains the vector _G1P0 if the object has no orientation.*/
      boost::math::quaternion<double> quatG1P0_abs(posquat(_P0) - posquat(q1));
      quatBuff = (1.0/rot1) * quatG1P0_abs * rot1;

      _G1P0x = quatBuff.R_component_2();
      _G1P0y = quatBuff.R_component_3();
      _G1P0z = quatBuff.R_component_4();

      _G2P0x = _P0->getValue(0);
      _G2P0y = _P0->getValue(1);
      _G2P0z = _P0->getValue(2);
    }
    else
    {
      _G1P0x = _P0->getValue(0);
      _G1P0y = _P0->getValue(1);
      _G1P0z = _P0->getValue(2);

      /*d2 is look as the ref frame. G2 is the origine, and d2 has no orientation.
        Where is P0 in the frame of d2 ?:*/
      /*Use initial value of q1 to place P0 in the absolute frame.*/
      /*quatG1P0_abs_ without any orientation*/
      ::boost::math::quaternion<double> quatG1P0_abs(posquat(_P0));

      /*quatBuff contains the vector G1P at the initial position*/
      quatBuff = rot1 * quatG1P0_abs / rot1;

      _G2P0x = q1->getValue(0) + quatBuff.R_component_2();
      _G2P0y = q1->getValue(1) + quatBuff.R_component_3();
      _G2P0z = q1->getValue(2) + quatBuff.R_component_4();
    }
  }

  // std::cout << "KneeJoint G1P0 :" << _G1P0x << " " << _G1P0y << " " << _G1P0z << std::endl;
  // std::cout << "KneeJoint G2P0 :" << _G2P0x << " " << _G2P0y << " " << _G2P0z << std::endl;
}

void KneeJointR::Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
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
  }
  //_jachq->display();
}

void KneeJointR::Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13)
{
  double df[20];
  double dfr0[20];
  double dfr1[20];
  double t104;
  double t109;
  double t11;
  double t5;
  double t6;
  double t60;
  double t70;
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
    _jachq->setValue(1, 0, 0.0);
    _jachq->setValue(1, 1, 1.0);
    _jachq->setValue(1, 2, 0.0);
    t60 = dfr0[17];
    _jachq->setValue(1, 3, t60 * q11 + dfr0[4]*q13 + 2.0 * q10 * _G1P0y);
    _jachq->setValue(1, 4, t60 * q10 + dfr1[6]*q12 + 2.0 * dfr0[0]*q11);
    t70 = dfr0[16];
    _jachq->setValue(1, 5, t70 * q13 + dfr1[6]*q11 + 2.0 * q12 * _G1P0y);
    _jachq->setValue(1, 6, t70 * q12 + dfr0[4]*q10 + 2.0 * dfr0[0]*q13);
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
  }
}

void KneeJointR::computeJachq(double time, Interaction& inter, SP::BlockVector q0)
{
  DEBUG_BEGIN("KneeJointR::computeJachq(double time, Interaction& inter,  SP::BlockVector q0) \n");
  
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
  DEBUG_END("KneeJointR::computeJachq(double time, Interaction& inter,  SP::BlockVector q0 ) \n");

}


void KneeJointR::DotJd1(double Xdot1, double Ydot1, double Zdot1,
                        double qdot10, double qdot11, double qdot12, double qdot13)
{

  double t1 = _G1P0y*qdot10;
  double t2 = _G1P0x*qdot10;
  double t4 = -t1+t2+_G1P0z*qdot12;
  double t7 = _G1P0z*qdot10;
  double t8 = _G1P0x*qdot11+_G1P0y*qdot12+t7;
  double t11 = t7-_G1P0x*qdot12+_G1P0y*qdot11;
  double t13 = -t2+_G1P0z*qdot11-t1;
  //double t17 = -_G2P0x*qdot20+_G2P0y*qdot23-_G2P0z*qdot22;
  //double t21 = -_G2P0x*qdot21-_G2P0y*qdot22-_G2P0z*qdot23;
  //double t25 = -_G2P0z*qdot20+_G2P0x*qdot22-_G2P0y*qdot21;
  //double t29 = _G2P0y*qdot20-_G2P0z*qdot21+_G2P0x*qdot23;
  _dotjachq->setValue(0,0, 0.0);
  _dotjachq->setValue(0,1, 0.0);
  _dotjachq->setValue(0,2, 0.0);
  _dotjachq->setValue(0,3, 2.0*t4);
  _dotjachq->setValue(0,4, 2.0*t8);
  _dotjachq->setValue(0,5, 2.0*t11);
  _dotjachq->setValue(0,6, 2.0*t13);

  _dotjachq->setValue(1,0, 0.0);
  _dotjachq->setValue(1,1, 0.0);
  _dotjachq->setValue(1,2, 0.0);
  _dotjachq->setValue(1,3, -2.0*t13);
  _dotjachq->setValue(1,4, -2.0*t11);
  _dotjachq->setValue(1,5, 2.0*t8);
  _dotjachq->setValue(1,6, 2.0*t4);

  _dotjachq->setValue(2,0, 0.0);
  _dotjachq->setValue(2,1, 0.0);
  _dotjachq->setValue(2,2, 0.0);
  _dotjachq->setValue(2,3, 2.0*t11);
  _dotjachq->setValue(2,4, -2.0*t13);
  _dotjachq->setValue(2,5, -2.0*t4);
  _dotjachq->setValue(2,6, 2.0*t8);


}


void KneeJointR::DotJd1d2(double Xdot1, double Ydot1, double Zdot1,
                          double qdot10, double qdot11, double qdot12, double qdot13,
                          double Xdot2, double Ydot2, double Zdot2,
                          double qdot20, double qdot21, double qdot22, double qdot23)
{
  double t1 = _G1P0y*qdot10;
  double t2 = _G1P0x*qdot10;
  double t4 = -t1+t2+_G1P0z*qdot12;
  double t7 = _G1P0z*qdot10;
  double t8 = _G1P0x*qdot11+_G1P0y*qdot12+t7;
  double t11 = t7-_G1P0x*qdot12+_G1P0y*qdot11;
  double t13 = -t2+_G1P0z*qdot11-t1;
  double t17 = -_G2P0x*qdot20+_G2P0y*qdot23-_G2P0z*qdot22;
  double t21 = -_G2P0x*qdot21-_G2P0y*qdot22-_G2P0z*qdot23;
  double t25 = -_G2P0z*qdot20+_G2P0x*qdot22-_G2P0y*qdot21;
  double t29 = _G2P0y*qdot20-_G2P0z*qdot21+_G2P0x*qdot23;
  _dotjachq->setValue(0,0, 0.0);
  _dotjachq->setValue(0,1, 0.0);
  _dotjachq->setValue(0,2, 0.0);
  _dotjachq->setValue(0,3, 2.0*t4);
  _dotjachq->setValue(0,4, 2.0*t8);
  _dotjachq->setValue(0,5, 2.0*t11);
  _dotjachq->setValue(0,6, 2.0*t13);
  _dotjachq->setValue(0,7, 0.0);
  _dotjachq->setValue(0,8, 0.0);
  _dotjachq->setValue(0,9, 0.0);
  _dotjachq->setValue(0,10, 2.0*t17);
  _dotjachq->setValue(0,11, 2.0*t21);
  _dotjachq->setValue(0,12, 2.0*t25);
  _dotjachq->setValue(0,13, 2.0*t29);
  _dotjachq->setValue(1,0, 0.0);
  _dotjachq->setValue(1,1, 0.0);
  _dotjachq->setValue(1,2, 0.0);
  _dotjachq->setValue(1,3, -2.0*t13);
  _dotjachq->setValue(1,4, -2.0*t11);
  _dotjachq->setValue(1,5, 2.0*t8);
  _dotjachq->setValue(1,6, 2.0*t4);
  _dotjachq->setValue(1,7, 0.0);
  _dotjachq->setValue(1,8, 0.0);
  _dotjachq->setValue(1,9, 0.0);
  _dotjachq->setValue(1,10, -2.0*t29);
  _dotjachq->setValue(1,11, -2.0*t25);
  _dotjachq->setValue(1,12, 2.0*t21);
  _dotjachq->setValue(1,13, 2.0*t17);
  _dotjachq->setValue(2,0, 0.0);
  _dotjachq->setValue(2,1, 0.0);
  _dotjachq->setValue(2,2, 0.0);
  _dotjachq->setValue(2,3, 2.0*t11);
  _dotjachq->setValue(2,4, -2.0*t13);
  _dotjachq->setValue(2,5, -2.0*t4);
  _dotjachq->setValue(2,6, 2.0*t8);
  _dotjachq->setValue(2,7, 0.0);
  _dotjachq->setValue(2,8, 0.0);
  _dotjachq->setValue(2,9, 0.0);
  _dotjachq->setValue(2,10, 2.0*t25);
  _dotjachq->setValue(2,11, -2.0*t29);
  _dotjachq->setValue(2,12, -2.0*t17);
  _dotjachq->setValue(2,13, 2.0*t21);

}
void KneeJointR::computeDotJachq(double time, BlockVector& workQ, BlockVector& workZ, BlockVector& workQdot)
{
  DEBUG_PRINT("KneeJointR::computeDotJachq(double time, Interaction& inter) starts \n");
  if (workQdot.numberOfBlocks()>1)
  {
    computeDotJachq(time, (workQdot.getAllVect())[0], (workQdot.getAllVect())[1]);
  }
  else
  {
    computeDotJachq(time, (workQdot.getAllVect())[0]);
  }
  DEBUG_PRINT("KneeJointR::computeDotJachq(double time, Interaction& inter) ends \n");
}

void KneeJointR::computeDotJachq(double time, SP::SiconosVector qdot1, SP::SiconosVector qdot2 )
{
  DEBUG_BEGIN("KneeJointR::computeDotJachq(double time, SP::SiconosVector qdot1, SP::SiconosVector qdot2) \n");
  _dotjachq->zero();
  
  double Xdot1 = qdot1->getValue(0);
  double Ydot1 = qdot1->getValue(1);
  double Zdot1 = qdot1->getValue(2);
  double qdot10 = qdot1->getValue(3);
  double qdot11 = qdot1->getValue(4);
  double qdot12 = qdot1->getValue(5);
  double qdot13 = qdot1->getValue(6);

  double Xdot2 = 0;
  double Ydot2 = 0;
  double Zdot2 = 0;
  double qdot20 = 1;
  double qdot21 = 0;
  double qdot22 = 0;
  double qdot23 = 0;

  if(qdot2)
  {
    Xdot2 = qdot2->getValue(0);
    Ydot2 = qdot2->getValue(1);
    Zdot2 = qdot2->getValue(2);
    qdot20 = qdot2->getValue(3);
    qdot21 = qdot2->getValue(4);
    qdot22 = qdot2->getValue(5);
    qdot23 = qdot2->getValue(6);
    DotJd1d2(Xdot1, Ydot1, Zdot1, qdot10, qdot11, qdot12, qdot13, Xdot2, Ydot2, Zdot2, qdot20, qdot21, qdot22, qdot23);
  }
  else
    DotJd1(Xdot1, Ydot1, Zdot1, qdot10, qdot11, qdot12, qdot13);

  DEBUG_END("KneeJointR::computeDotJachq(double time, SP::SiconosVector qdot1, SP::SiconosVector qdot2 ) \n");
}


/** to compute p
 *  \param double : current time
 *  \param unsigned int: "derivative" order of lambda used to compute input
 */
double KneeJointR::Hx(double X1, double Y1, double Z1, double  q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)

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
double KneeJointR::Hy(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)

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
double KneeJointR::Hz(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)

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

void KneeJointR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  DEBUG_BEGIN("KneeJointR::computeh(double time, BlockVector& q0, SiconosVector& y)\n");
  DEBUG_EXPR(q0.display());

  double X1 = q0.getValue(0);
  double Y1 = q0.getValue(1);
  double Z1 = q0.getValue(2);
  double q10 = q0.getValue(3);
  double q11 = q0.getValue(4);
  double q12 = q0.getValue(5);
  double q13 = q0.getValue(6);
  DEBUG_PRINTF("X1 = %12.8e,\t Y1 = %12.8e,\t Z1 = %12.8e,\n",X1,Y1,Z1);
  DEBUG_PRINTF("q10 = %12.8e,\t q11 = %12.8e,\t q12 = %12.8e,\t q13 = %12.8e,\n",q10,q11,q12,q13);
  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;
  if(q0.numberOfBlocks()>1)
  {
    // SP::SiconosVector x2 = _d2->q();
    // DEBUG_EXPR( _d2->q()->display(););
    X2 = q0.getValue(7);
    Y2 = q0.getValue(8);
    Z2 = q0.getValue(9);
    q20 = q0.getValue(10);
    q21 = q0.getValue(11);
    q22 = q0.getValue(12);
    q23 = q0.getValue(13);
  }
  y.setValue(0, Hx(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(1, Hy(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  y.setValue(2, Hz(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23));
  DEBUG_EXPR(y.display());
  DEBUG_END("KneeJointR::computeh(double time, BlockVector& q0, SiconosVector& y)\n");
    
}
