/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "odepack.h"

extern "C" void CNAME(f1)(integer *sizeOfX, doublereal *time, doublereal *x, doublereal *xdot);
extern "C" void CNAME(f1)(integer *sizeOfX, doublereal *time, doublereal *x, doublereal *xdot)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output parameter : xdot (pointer to Xdot vector)
   */

  //        printf("Call of the function 'f1' of the file funC.\n");

  if(sizeOfX[0] == 2)
  {
    xdot[0] = x[1];
    xdot[1] = 3.0 * (1.0 - x[0] * x[0]) * x[1] - x[0];
  }
  else
  {
    printf("OscillatorPlugin:vectorField --- Bad size of x. %ld -- %p\n", (long int) *sizeOfX, sizeOfX);
  }

}

extern "C" void CNAME(jac1)(integer *sizeOfX, doublereal *time, doublereal *x, integer* ml, integer *mu,  doublereal *jacob, integer *nrowpd);
extern "C" void CNAME(jac1)(integer *sizeOfX, doublereal *time, doublereal *x, integer* ml, integer *mu,  doublereal *jacob, integer *nrowpd)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */

  //      printf("Call of the function 'jac1' of the  the file funC.c .\n");

  if(*sizeOfX == 2)
  {
    jacob[0] = 0.0;
    jacob[1] = -6.0 * x[0] * x[1] - 1.0;
    jacob[2] = 1.0;
    jacob[3] =  3.0 * (1.0 - x[0] * x[0]);
  }
  else
  {
    printf("OscillatorPlugin:computeJacobianX --- Bad size of x. \n");
  }
}
