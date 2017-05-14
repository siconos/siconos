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

PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SiconosVector P, SP::SiconosVector A): KneeJointR(d1, d2, P)
{
  _A.reset( new SiconosVector(*A) );

  ::boost::math::quaternion<double> q1((*d1->q())(3), (*d1->q())(4),
                                       (*d1->q())(5), (*d1->q())(6));
  ::boost::math::quaternion<double> q2((*d2->q())(3), (*d2->q())(4),
                                       (*d2->q())(5), (*d2->q())(6));
  ::boost::math::quaternion<double> cq2q10(1.0 / q2 * q1);

  _cq2q101 = cq2q10.R_component_1();
  _cq2q102 = cq2q10.R_component_2();
  _cq2q103 = cq2q10.R_component_3();
  _cq2q104 = cq2q10.R_component_4();

  buildA1A2();

  _initial_AscalA1 = AscalA1((*d1->q())(3), (*d1->q())(4),
                             (*d1->q())(5), (*d1->q())(6),
                             (*d2->q())(3), (*d2->q())(4),
                             (*d2->q())(5), (*d2->q())(6));

  _initial_AscalA2 = AscalA2((*d1->q())(3), (*d1->q())(4),
                             (*d1->q())(5), (*d1->q())(6),
                             (*d2->q())(3), (*d2->q())(4),
                             (*d2->q())(5), (*d2->q())(6));
}

/* constructor,
   \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
   \param a SP::SiconosVector P0, see KneeJointR documentation.
   \param a SP::SiconosVector A, axis in the frame of the object.
   \param a bool, used only by the KneeJointR constructor see KneeJointR documentation.
*/
PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::SiconosVector P0, SP::SiconosVector A, bool absolutRef): KneeJointR(d1, P0, absolutRef)
{
  ::boost::math::quaternion<double> q1((*d1->q())(3), (*d1->q())(4),
                                       (*d1->q())(5), (*d1->q())(6));
  ::boost::math::quaternion<double> cq2q10(q1);

  _cq2q101 = cq2q10.R_component_1();
  _cq2q102 = cq2q10.R_component_2();
  _cq2q103 = cq2q10.R_component_3();
  _cq2q104 = cq2q10.R_component_4();

  _A.reset( new SiconosVector(*A) );
  buildA1A2();

  _initial_AscalA1 = AscalA1((*d1->q())(3), (*d1->q())(4),
                             (*d1->q())(5), (*d1->q())(6),
                             1, 0, 0, 0);

  _initial_AscalA2 = AscalA2((*d1->q())(3), (*d1->q())(4),
                             (*d1->q())(5), (*d1->q())(6),
                             1, 0, 0, 0);
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






double PivotJointR::AscalA1(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23)
{
  /*
   * The angle between A1 and rotation q2-to-q1 must be zero,
   * (taking into account original difference in orientation q2to1).
   *
   * sympy expression:
   * AscalA1 = np.dot(A1,qmul(qinv(qmul(q2,cq2q10)),q1)) - initial_AscalA1
   */

  return _A1x*(q10*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
               + q11*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
               - q12*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
               + q13*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21))
    + _A1y*(q10*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
            + q11*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
            + q12*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
            - q13*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22))
    + _A1z*(q10*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
            - q11*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
            + q12*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
            + q13*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23));
}

double PivotJointR::AscalA2(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23)
{
  /*
   * The angle between A2 and rotation q2-to-q1 must be zero,
   * (taking into account original difference in orientation q2to1).
   *
   * sympy expression:
   * AscalA2 = np.dot(A2,qmul(qinv(qmul(q2,cq2q10)),q1)) - initial_AscalA2
   */

  return _A2x*(q10*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
               + q11*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
               - q12*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
               + q13*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21))
    + _A2y*(q10*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
            + q11*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
            + q12*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23)
            - q13*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22))
    + _A2z*(q10*(-_cq2q101*q23 + _cq2q102*q22 - _cq2q103*q21 - _cq2q104*q20)
            - q11*(-_cq2q101*q22 - _cq2q102*q23 - _cq2q103*q20 + _cq2q104*q21)
            + q12*(-_cq2q101*q21 - _cq2q102*q20 + _cq2q103*q23 - _cq2q104*q22)
            + q13*(_cq2q101*q20 - _cq2q102*q21 - _cq2q103*q22 - _cq2q104*q23));
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

  y.setValue(3, AscalA1(q10, q11, q12, q13, q20, q21, q22, q23) - _initial_AscalA1);
  y.setValue(4, AscalA2(q10, q11, q12, q13, q20, q21, q22, q23) - _initial_AscalA2);

}
