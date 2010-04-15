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

PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SimpleVector P, SP::SimpleVector A): KneeJointR(d1, d2, P)
{
  SP::SiconosVector q1 = d1->q0();
  ::boost::math::quaternion<float>    quat1(q1->getValue(3), q1->getValue(4), q1->getValue(5), q1->getValue(6));
  ::boost::math::quaternion<float>    quatA(0, A->getValue(0), A->getValue(1), A->getValue(2));
  ::boost::math::quaternion<float>    quatBuff(0, 0, 0, 0);
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
PivotJointR::PivotJointR(SP::NewtonEulerDS d1, SP::SimpleVector P0, SP::SimpleVector A): KneeJointR(d1, P0)
{

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
  if (fabs(_Ax) > fabs(_Ay))
  {
    if (fabs(_Ax) > fabs(_Az)) /*_Ax is the bigest*/
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
  else if (fabs(_Ay) > fabs(_Az))  /*_Ay is the bigest*/
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
  KneeJointR::computeh(t);
  SP::SiconosVector x1 = _d1->q();
  //std::cout<<"PivotJoint computeH d1->q:\n";
  //x1->display();
  double q10 = x1->getValue(3);
  double q11 = x1->getValue(4);
  double q12 = x1->getValue(5);
  double q13 = x1->getValue(6);
  double q20 = 1;
  double q21 = 0;
  double q22 = 0;
  double q23 = 0;
  if (_d2)
  {
    SP::SiconosVector x2 = _d2->q();
    q20 = x2->getValue(3);
    q21 = x2->getValue(4);
    q22 = x2->getValue(5);
    q23 = x2->getValue(6);
  }

  SP::SiconosVector y = interaction()->y(0);
  y->setValue(3, AscalA1(q10, q11, q12, q13, q20, q21, q22, q23));
  y->setValue(4, AscalA2(q10, q11, q12, q13, q20, q21, q22, q23));

  std::cout << "PivotJoint computeH:\n";
  y->display();
}


