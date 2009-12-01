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
#include "myDS.h"





MyDS::MyDS(const SiconosVector& x0): FirstOrderNonLinearDS(x0)
{
  _jacobianXF.reset(new SimpleMatrix(4, 4));
  mf.reset(new SimpleVector(4));
  _M.reset(new SimpleMatrix(4, 4));
  _M->eye();
}

void MyDS::computeF(double t)
{

  SP::SimpleMatrix Q(new SimpleMatrix(2, 2));
  Q->eye();

  SP::SimpleVector QX(new SimpleVector(2));
  SP::SimpleVector X(new SimpleVector(2));

  X->setValue(0, (_x[0]->getValue(0) - 2.0));
  X->setValue(1, (_x[0]->getValue(1) + 1.0));

  prod(*Q, *X, *QX, true);



  SP::SimpleMatrix K1(new SimpleMatrix(2, 2));
  K1->setValue(0, 0, 0.0);
  K1->setValue(0, 1, -1.0 / 2.0);
  K1->setValue(1, 0, +1.0 / 2.0);
  K1->setValue(1, 1, +1.0);

  SP::SimpleVector K1P(new SimpleVector(2));
  SP::SimpleVector P(new SimpleVector(2));
  P->setValue(0, _x[0]->getValue(2));
  P->setValue(1, _x[0]->getValue(3));
  prod(*K1, *P, *K1P, true);


  mf->setValue(0, -1.0 / 2.0 * _x[0]->getValue(1) - 1.0 / 2.0);
  mf->setValue(1, -1.0 / 2.0 * _x[0]->getValue(0) - _x[0]->getValue(1));
  mf->setValue(2, -QX->getValue(0) + K1P->getValue(0));
  mf->setValue(2, -QX->getValue(1) + K1P->getValue(1));



}
void  MyDS::computeF(double t, SP::SiconosVector _xvalue)
{

  SP::SimpleMatrix Q(new SimpleMatrix(2, 2));
  Q->eye();

  SP::SimpleVector QX(new SimpleVector(2));
  SP::SimpleVector X(new SimpleVector(2));

  X->setValue(0, (_xvalue->getValue(0) - 2.0));
  X->setValue(1, (_xvalue->getValue(1) + 1.0));

  prod(*Q, *X, *QX, true);



  SP::SimpleMatrix K1(new SimpleMatrix(2, 2));
  K1->setValue(0, 0, 0.0);
  K1->setValue(0, 1, -1.0 / 2.0);
  K1->setValue(1, 0, +1.0 / 2.0);
  K1->setValue(1, 1, +1.0);

  SP::SimpleVector K1P(new SimpleVector(2));
  SP::SimpleVector P(new SimpleVector(2));
  P->setValue(0, _xvalue->getValue(2));
  P->setValue(1, _xvalue->getValue(3));
  prod(*K1, *P, *K1P, true);


  mf->setValue(0, -1.0 / 2.0 * _xvalue->getValue(1) - 1.0 / 2.0);
  mf->setValue(1, -1.0 / 2.0 * _xvalue->getValue(0) - _xvalue->getValue(1));
  mf->setValue(2, -QX->getValue(0) + K1P->getValue(0));
  mf->setValue(2, -QX->getValue(1) + K1P->getValue(1));

}

void MyDS::computeJacobianXF(double t, bool  b)
{




  SP::SimpleMatrix Q(new SimpleMatrix(2, 2));
  Q->eye();
  SP::SimpleMatrix K1(new SimpleMatrix(2, 2));
  K1->setValue(0, 0, 0.0);
  K1->setValue(0, 1, -1.0 / 2.0);
  K1->setValue(1, 0, +1.0 / 2.0);
  K1->setValue(1, 1, +1.0);



  _jacobianXF->setValue(0, 0, 0);
  _jacobianXF->setValue(0, 1, -1.0 / 2.0);
  _jacobianXF->setValue(0, 2, 0.0);
  _jacobianXF->setValue(0, 3, 0.0);
  _jacobianXF->setValue(1, 0, 1.0 / 2.0);
  _jacobianXF->setValue(1, 1, -1.0);
  _jacobianXF->setValue(1, 2, 0.0);
  _jacobianXF->setValue(1, 3, 0.0);
  _jacobianXF->setValue(2, 0, -Q->getValue(0, 0));
  _jacobianXF->setValue(2, 1, -Q->getValue(0, 1));
  _jacobianXF->setValue(2, 2, -Q->getValue(1, 0));
  _jacobianXF->setValue(2, 3, -Q->getValue(1, 1));
  _jacobianXF->setValue(3, 0, K1->getValue(0, 0));
  _jacobianXF->setValue(3, 1, K1->getValue(0, 1));
  _jacobianXF->setValue(3, 2, K1->getValue(1, 0));
  _jacobianXF->setValue(3, 3, K1->getValue(1, 1));

}

void MyDS::computeJacobianXF(double t, SP::SiconosVector v)
{

  SP::SimpleMatrix Q(new SimpleMatrix(2, 2));
  Q->eye();
  SP::SimpleMatrix K1(new SimpleMatrix(2, 2));
  K1->setValue(0, 0, 0.0);
  K1->setValue(0, 1, -1.0 / 2.0);
  K1->setValue(1, 0, +1.0 / 2.0);
  K1->setValue(1, 1, +1.0);



  _jacobianXF->setValue(0, 0, 0);
  _jacobianXF->setValue(0, 1, -1.0 / 2.0);
  _jacobianXF->setValue(0, 2, 0.0);
  _jacobianXF->setValue(0, 3, 0.0);
  _jacobianXF->setValue(1, 0, 1.0 / 2.0);
  _jacobianXF->setValue(1, 1, -1.0);
  _jacobianXF->setValue(1, 2, 0.0);
  _jacobianXF->setValue(1, 3, 0.0);
  _jacobianXF->setValue(2, 0, -Q->getValue(0, 0));
  _jacobianXF->setValue(2, 1, -Q->getValue(0, 1));
  _jacobianXF->setValue(2, 2, -Q->getValue(1, 0));
  _jacobianXF->setValue(2, 3, -Q->getValue(1, 1));
  _jacobianXF->setValue(3, 0, K1->getValue(0, 0));
  _jacobianXF->setValue(3, 1, K1->getValue(0, 1));
  _jacobianXF->setValue(3, 2, K1->getValue(1, 0));
  _jacobianXF->setValue(3, 3, K1->getValue(1, 1));

}

void MyDS::computeRhs(double t, bool  b)
{
  ;
}
void MyDS::resetNonSmoothPart()
{
}
