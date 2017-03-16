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

PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SiconosVector P, SP::SiconosVector A): KneeJointR(d1, d2, P)
{
  // SP::SiconosVector q1 = d1->q0();
  // ::boost::math::quaternion<double>    quat1(q1->getValue(3),-q1->getValue(4),-q1->getValue(5),-q1->getValue(6));
  // ::boost::math::quaternion<double>    quatA(0,A->getValue(0),A->getValue(1),A->getValue(2));
  // ::boost::math::quaternion<double>    quatBuff(0,0,0,0);
  // /*calcul of axis _A*/
  // quatBuff=quat1*quatA/quat1;
  // _Ax=quatBuff.R_component_2();
  // _Ay=quatBuff.R_component_3();
  // _Az=quatBuff.R_component_4();
  _A.reset( new SiconosVector(*A) );
  buildA1A2();
}
/* constructor,
   \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
   \param a SP::SiconosVector P0, see KneeJointR documentation.
   \param a SP::SiconosVector A, axis in the frame of the object.
   \param a bool, used only by the KneeJointR constructor see KneeJointR documentation.
*/
PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::SiconosVector P0, SP::SiconosVector A, bool absolutRef): KneeJointR(d1, P0, absolutRef)
{
  _A.reset( new SiconosVector(*A) );
  buildA1A2();
}

void PivotJointR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  KneeJointR::initComponents(inter,DSlink,workV,workM);
  //if (_d2){
  //proj_with_q  _jachqProj.reset(new SimpleMatrix(7,14));
  //proj_with_q    _yProj.reset(new SiconosVector(7));
  //_yProj.reset(new SiconosVector(5));
  // }else{
  //proj_with_q  _jachqProj.reset(new SimpleMatrix(6,7));
  //proj_with_q    _yProj.reset(new SiconosVector(6));
  //_yProj.reset(new SiconosVector(5));
  // }
  //proj_with_q  _jachqProj->zero();
  //_yProj->zero();
}
void PivotJointR::buildA1A2()
{
  double _Ax = (*_A)(0);
  double _Ay = (*_A)(1);
  double _Az = (*_A)(2);
  if (orthoBaseFromVector(&_Ax, &_Ay, &_Az, &_A1x, &_A1y, &_A1z, &_A2x, &_A2y, &_A2z))
    RuntimeException::selfThrow("PivotJointR::initComponents. Problem in calling orthoBaseFromVector");
  assert(fabs(_A1x * _Ax + _A1y * _Ay + _A1z * _Az) < 1e-9 && "PivotJoint, _A1 wrong\n");
  assert(fabs(_A2x * _Ax + _A2y * _Ay + _A2z * _Az) < 1e-9 && "PivotJoint, _A2 wrong\n");
  assert(fabs(_A1x * _A2x + _A1y * _A2y + _A1z * _A2z) < 1e-9 && "PivotJoint, _A12 wrong\n");
  std::cout << "JointPivot: _A1x _A1y _A1z :" << _A1x << " " << _A1y << " " << _A1z << std::endl;
  std::cout << "JointPivot: _A2x _A2y _A2z :" << _A2x << " " << _A2y << " " << _A2z << std::endl;
}
void PivotJointR::Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23)
{
  KneeJointR::Jd1d2(X1, Y1, Z1, q10, q11, q12, q13, X2, Y2, Z2, q20, q21, q22, q23);

  _jachq->setValue(3, 0, 0);
  _jachq->setValue(3, 1, 0);
  _jachq->setValue(3, 2, 0);
  _jachq->setValue(3, 3, _A1x * (-q21) + _A1y * (-q22) + _A1z * (-q23));
  _jachq->setValue(3, 4, _A1x * (q20) + _A1y * (-q23) + _A1z * (q22));
  _jachq->setValue(3, 5, _A1x * (q23) + _A1y * (q20) + _A1z * (-q21));
  _jachq->setValue(3, 6, _A1x * (-q22) + _A1y * (q21) + _A1z * (q20));
  _jachq->setValue(3, 7, 0);
  _jachq->setValue(3, 8, 0);
  _jachq->setValue(3, 9, 0);
  _jachq->setValue(3, 10, _A1x * (q11) + _A1y * (q12) + _A1z * (q13));
  _jachq->setValue(3, 11, _A1x * (-q10) + _A1y * (q13) + _A1z * (-q12));
  _jachq->setValue(3, 12, _A1x * (-q13) + _A1y * (-q10) + _A1z * (q11));
  _jachq->setValue(3, 13, _A1x * (q12) + _A1y * (-q11) + _A1z * (-q10));

  _jachq->setValue(4, 0, 0);
  _jachq->setValue(4, 1, 0);
  _jachq->setValue(4, 2, 0);
  _jachq->setValue(4, 3, _A2x * (-q21) + _A2y * (-q22) + _A2z * (-q23));
  _jachq->setValue(4, 4, _A2x * (q20) + _A2y * (-q23) + _A2z * (q22));
  _jachq->setValue(4, 5, _A2x * (q23) + _A2y * (q20) + _A2z * (-q21));
  _jachq->setValue(4, 6, _A2x * (-q22) + _A2y * (q21) + _A2z * (q20));
  _jachq->setValue(4, 7, 0);
  _jachq->setValue(4, 8, 0);
  _jachq->setValue(4, 9, 0);
  _jachq->setValue(4, 10, _A2x * (q11) + _A2y * (q12) + _A2z * (q13));
  _jachq->setValue(4, 11, _A2x * (-q10) + _A2y * (q13) + _A2z * (-q12));
  _jachq->setValue(4, 12, _A2x * (-q13) + _A2y * (-q10) + _A2z * (q11));
  _jachq->setValue(4, 13, _A2x * (q12) + _A2y * (-q11) + _A2z * (-q10));

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






double PivotJointR::AscalA1(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23)
{
  ::boost::math::quaternion<double>    quat1(q10, q11, q12, q13);
  ::boost::math::quaternion<double>    quat2_inv(q20, -q21, -q22, -q23);
  ::boost::math::quaternion<double>    quatBuff = quat2_inv * quat1;
  double aX = quatBuff.R_component_2();
  double aY = quatBuff.R_component_3();
  double aZ = quatBuff.R_component_4();
  return _A1x * aX + _A1y * aY + _A1z * aZ;
}
double PivotJointR::AscalA2(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23)
{
  ::boost::math::quaternion<double>    quat1(q10, q11, q12, q13);
  ::boost::math::quaternion<double>    quat2_inv(q20, -q21, -q22, -q23);
  ::boost::math::quaternion<double>    quatBuff = quat2_inv * quat1;
  double aX = quatBuff.R_component_2();
  double aY = quatBuff.R_component_3();
  double aZ = quatBuff.R_component_4();
  return _A2x * aX + _A2y * aY + _A2z * aZ;
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
  if(q0.getNumberOfBlocks()>1)
  {
    q20 = q0.getValue(10);
    q21 = q0.getValue(11);
    q22 = q0.getValue(12);
    q23 = q0.getValue(13);
  }

  y.setValue(3, AscalA1(q10, q11, q12, q13, q20, q21, q22, q23));
  y.setValue(4, AscalA2(q10, q11, q12, q13, q20, q21, q22, q23));

}
