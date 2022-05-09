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
#include "odepack.h"
#include <stdio.h>

extern "C" void CNAME(f1)(integer *sizeOfX, doublereal *time, doublereal *x, doublereal *xdot);
extern "C" void CNAME(jac1)(integer *sizeOfX, doublereal *time, doublereal *x, integer* ml, integer *mu,  doublereal *jacob, integer *nrowpd);

int main(void)
{

  const doublereal tout1 = 1.39283880203;
  const doublereal dtout = 2.214773875;
  integer  itol = 1;
  doublereal  rtol = 0.0;
  doublereal  atol = 1.0e-6;
  integer  lrw = 697;
  integer  liw = 45;
  integer  iopt = 0;
  //
  //First problem
  //
  integer neq = 2;
  integer nout = 4;


  doublereal t = 0.0;
  doublereal y[2] = {2.0, 0.0};
  doublereal  rwork[697];
  integer  iwork[45];

  integer   itask = 1;
  integer  istate = 1;
  doublereal  tout = tout1;

  integer mf ;
  doublereal hu;
  integer nqu;

  for(int meth = 1; meth <= 2; meth++)
  {
    for(int miter = 0; miter <= 3; miter++)
    {
      mf = 10 * meth + miter;
      printf(" Solution with mf = %3d\n", (int)mf);
      printf("     t               x               xdot       nq      h\n");

      t = 0.0;
      y[0] = 2.0;
      y[1] = 0.0;


      itask = 1;
      istate = 1;
      tout = tout1;

      for(int i = 1; i <= nout; i++)
      {
        CNAME(dlsode)(CNAME(f1), &neq, y, &t, &tout, &itol, &rtol, &atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, CNAME(jac1), &mf);
        hu = rwork[10];
        nqu = iwork[13];
        printf("%15.4e\t%16.4e\t%14.2e\t%5d\t%14.2e\n", t, y[0], y[1], (int)nqu, hu);
        tout = tout + dtout;
      }

    }
  }


}
