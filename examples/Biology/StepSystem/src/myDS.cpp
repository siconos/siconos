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
#include "myDS.h"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

MyDS::MyDS(SP::SiconosVector x0): FirstOrderNonLinearDS(x0)
{
  _jacobianfx.reset(new SimpleMatrix(2, 2));
  _f.reset(new SiconosVector(2));

  _M.reset(new SimpleMatrix(2, 2));
  _M->zero();
  _M->setValue(0, 0, 1);
  _M->setValue(1, 1, 1);
}

void  MyDS::computef(double t, SP::SiconosVector x)
{
  //SP::SiconosVector x=x();
  _f->setValue(0, -4.5 * x->getValue(0));
  _f->setValue(1, -1.5 * x->getValue(1));
  DEBUG_PRINT("MyDS::computeF");
  DEBUG_EXPR(x->display(););
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

void MyDS::computeJacobianfx(double t, SP::SiconosVector x)
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
