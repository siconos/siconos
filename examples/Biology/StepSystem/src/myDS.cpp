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

#define SICONOS_DEBUG

MyDS::MyDS(SP::SiconosVector x0): FirstOrderNonLinearDS(x0)
{
  _jacobianfx.reset(new SimpleMatrix(2, 2));
  _f.reset(new SiconosVector(2));

  _M.reset(new SimpleMatrix(2, 2));
  _M->zero();
  _M->setValue(0, 0, 1);
  _M->setValue(1, 1, 1);
}

void  MyDS::computef(double t)
{
  //SP::SiconosVector x=x();
  _f->setValue(0, -4.5 * (x()->getValue(0)));
  _f->setValue(1, -1.5 * (x()->getValue(1)));

  /*
  #ifdef SICONOS_DEBUG
    std::cout<<"MyDS::computeF with x="<<std::endl;
    x()->display();
    std::cout<<std::endl;
    std::cout<<"F(x)="<<std::endl;
    _f->display();
    std::cout<<std::endl;
  #endif
  */

}

void MyDS::computeJacobianfx(double t, bool)
{
  _jacobianfx->setValue(0, 0, -4.5);
  _jacobianfx->setValue(1, 0, 0);
  _jacobianfx->setValue(0, 1, 0);
  _jacobianfx->setValue(1, 1, -1.5);

  /*
  #ifdef SICONOS_DEBUG
    std::cout<<"MyDS::computeJacobianfx."<<std::endl;
  std::cout<<"Nabla f="<<std::endl;
    _jacobianfx->display();
    std::cout<<std::endl;
  #endif
  */

}

// void MyDS::computeRhs(double t, bool  b)
// {
//   ;
// }


void MyDS::resetNonSmoothPart()
{
}
