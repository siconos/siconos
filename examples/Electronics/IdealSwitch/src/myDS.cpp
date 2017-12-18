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



MyDS::MyDS(SP::SiconosVector x0): FirstOrderNonLinearDS(x0)
{
  _jacobianfx.reset(new SimpleMatrix(1, 1));
  _f.reset(new SiconosVector(1));
  _M.reset(new SimpleMatrix(1, 1));
  _M->eye();
}

void MyDS::computeF(double t)
{
  _f->setValue(0, 0);
}
void  MyDS::computeF(double, SP::SiconosVector)
{
  _f->setValue(0, 0);
}

void MyDS::computeJacobianfx(double t, bool  b)
{
  _jacobianfx->setValue(0, 0, 0);
}

void MyDS::computeJacobianfx(double t, SP::SiconosVector v)
{
  _jacobianfx->setValue(0, 0, 0);
}

void MyDS::computeRhs(double t, bool  b)
{
  ;
}
void MyDS::resetNonSmoothPart(unsigned int level)
{
  _r->zero();
}
