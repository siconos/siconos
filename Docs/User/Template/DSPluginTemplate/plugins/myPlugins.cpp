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
#include <stdio.h>

extern "C" void myF(double time, unsigned int sizeX, const double* x, double* f, unsigned int sizeZ, double* z)
{
  // Do your stuff ...

  printf("call myF ... \n");

}

extern "C" void myJacobianF(double time, unsigned int sizeX, const double* x, double* jacobF, unsigned int sizeZ, double* z)
{
  // Do your stuff ...

  printf("call myJacobianF ... \n ");

}

extern "C" void myH(unsigned int sizeX, const double* x, unsigned int sizeY, double* Y, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myH ... \n ");
}

extern "C" void myG(unsigned int sizeY, const double* lambda, unsigned int sizeX, double* g, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myG ... \n ");
}

extern "C" void myJacobianH(unsigned int sizeX, const double* x, unsigned int sizeY, double* jacob, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myJacobianH ... \n ");
}

extern "C" void myJacobianG(unsigned int sizeY, const double* lambda, unsigned int sizeX, double* jacob, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myJacobianG ... \n ");
}

