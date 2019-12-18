/* Siconos-Numerics, Copyright INRIA 2005-2018.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
