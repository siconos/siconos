/* Siconos-Kernel , Copyright INRIA 2005-2011.
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



//#define SICONOS_DEBUG

MyDS::MyDS(SP::SiconosVector x0): FirstOrderNonLinearDS(x0)
{
  _jacobianfx.reset(new SimpleMatrix(4, 4));
  _f.reset(new SiconosVector(4));
  _M.reset(new SimpleMatrix(4, 4));
  _M->eye();

  Q = new SimpleMatrix(2, 2);
  Q->eye();
  K1 = new SimpleMatrix(2, 2);
  K1->setValue(0, 0, 0.0);
  K1->setValue(0, 1, 1.0 / 2.0);
  K1->setValue(1, 0, -1.0 / 2.0);
  K1->setValue(1, 1, +1.0);
}

MyDS::~MyDS()
{
  delete Q;
  delete K1;
}


void MyDS::computef(double t)
{


  SP::SiconosVector QX(new SiconosVector(2));
  SP::SiconosVector X(new SiconosVector(2));

  X->setValue(0, (_x[0]->getValue(0) - 2.0));
  X->setValue(1, (_x[0]->getValue(1) + 1.0));

  prod(*Q, *X, *QX, true);

  SP::SiconosVector K1P(new SiconosVector(2));
  SP::SiconosVector P(new SiconosVector(2));
  P->setValue(0, _x[0]->getValue(2));
  P->setValue(1, _x[0]->getValue(3));
  prod(*K1, *P, *K1P, true);


  SP::SiconosVector alphatmp(new SiconosVector(2));

  alpha(t, _x[0], alphatmp);


  _f->setValue(0, alphatmp->getValue(0));
  _f->setValue(1, alphatmp->getValue(1));
  _f->setValue(2, -QX->getValue(0) + K1P->getValue(0));
  _f->setValue(3, -QX->getValue(1) + K1P->getValue(1));
}
void  MyDS::computef(double t, SiconosVector& xvalue) {}

void MyDS::computeJacobianfx(double t, bool  b)
{




  SP::SiconosMatrix jacXalpha(new SimpleMatrix(2, 2));

  JacobianXalpha(t, _x[0], jacXalpha);

  _jacobianfx->setValue(0, 0, jacXalpha->getValue(0, 0));
  _jacobianfx->setValue(0, 1, jacXalpha->getValue(0, 1));
  _jacobianfx->setValue(0, 2, 0.0);
  _jacobianfx->setValue(0, 3, 0.0);
  _jacobianfx->setValue(1, 0, jacXalpha->getValue(1, 0));
  _jacobianfx->setValue(1, 1, jacXalpha->getValue(1, 1));
  _jacobianfx->setValue(1, 2, 0.0);
  _jacobianfx->setValue(1, 3, 0.0);
  _jacobianfx->setValue(2, 0, -Q->getValue(0, 0));
  _jacobianfx->setValue(2, 1, -Q->getValue(0, 1));
  _jacobianfx->setValue(2, 2, K1->getValue(0, 0));
  _jacobianfx->setValue(2, 3, K1->getValue(0, 1));
  _jacobianfx->setValue(3, 0, -Q->getValue(1, 0));
  _jacobianfx->setValue(3, 1, -Q->getValue(1, 1));
  _jacobianfx->setValue(3, 2, K1->getValue(1, 0));
  _jacobianfx->setValue(3, 3, K1->getValue(1, 1));



}

void MyDS::computeJacobianfx(double t, const SiconosVector& v) {}

void MyDS::alpha(double t, SP::SiconosVector _xvalue, SP::SiconosVector _alpha)
{
  _alpha->setValue(0, 1.0 / 2.0 * _xvalue->getValue(1) + 1.0 / 2.0) ;
  _alpha->setValue(1,  -1.0 / 2.0 * _xvalue->getValue(0) - _xvalue->getValue(1)) ;

}



void MyDS::JacobianXalpha(double t, SP::SiconosVector _xvalue, SP::SiconosMatrix JacXalpha)
{

  JacXalpha->setValue(0, 0, 0.0) ;
  JacXalpha->setValue(0, 1, 1.0 / 2.0) ;
  JacXalpha->setValue(1, 0, -1.0 / 2.0) ;
  JacXalpha->setValue(1, 1, -1.0) ;

#ifdef SICONOS_DEBUG
  std::cout << "JacXalpha\n" << std::endl;;
  JacXalpha->display();
#endif

}
