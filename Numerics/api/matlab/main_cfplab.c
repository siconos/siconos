/* Siconos-Numerics version 2.1.0, Copyright INRIA 2005-2006.
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
/*!\file main_cfplab.c


 This subroutine allows the primal resolution of contact problems with friction.

  Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

 This system of equations and inequalities is solved thanks to @ref pfc solvers.
 The routine's call is due to the function solve_cfp.c.

  This subroutine is a matlab interface that allows you to use an @ref pfc solver.
  You only have to enter the M matrix, the q vector, the maximum of iterations required, the tolerance value, the name of the solver 'Gsnl', 'Lemke', 'Gcp', or 'Latin', the friction coefficient $\mu$ and a search direction for the Latin  if you wanna use it .
\author Nineb Sheherazade.

*/


#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SiconosNumerics_lab.h"
#include "solve_cfp.c"
#include "cfp_gsnl.c"
#include "cfp_gcp.c"
#include "cfp_latin.c"
#include "projf.c"
#include "projc.c"



double ddot_(int *, double [], int *, double [], int*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double *vec, *qq, *z, *w, *it_end, *info;
  int mrows, ncols, *nn, n, nl, nc, toto;
  int i, j, *itt;
  double *tol, *itermax, *k_lat, *mu;
  double *res;
  char *Str, *mot1, *mot2, *mot3, *mot4;
  methode meth_cfp;



  mot1 = "Gsnl";
  mot2 = "Gcp";
  mot3 = "Latin";


  if (nlhs != 3)
  {
    mexErrMsgTxt("3 output required.");
  }


  if (nrhs <= 5)
  {
    mexErrMsgTxt("on less 6 input required.");
  }

  vec = mxGetPr(prhs[0]);
  qq = mxGetPr(prhs[1]);
  itermax = mxGetPr(prhs[2]);
  tol = mxGetPr(prhs[3]);
  mu = mxGetPr(prhs[4]);



  if (nrhs == 7)
  {
    k_lat = mxGetPr(prhs[6]);
  }
  else *k_lat = 1.;



  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);



  nl = mxGetM(prhs[5]);
  nc = mxGetN(prhs[5]);
  Str = mxCalloc(nl * nc + 1, 8);
  mxGetString(prhs[5], Str, nl * nc + 1);




  meth_cfp.cfp.nom_method = Str;
  meth_cfp.cfp.itermax = *itermax;
  meth_cfp.cfp.k_latin = *k_lat;
  meth_cfp.cfp.tol = *tol;
  meth_cfp.cfp.mu = *mu;



  plhs[0] = mxCreateDoubleMatrix(mrows, ncols, mxREAL); /* z */
  plhs[1] = mxCreateDoubleMatrix(mrows, ncols, mxREAL); /* w */
  plhs[2] = mxCreateDoubleMatrix(ncols, ncols, mxREAL); /* info */



  z = mxGetPr(plhs[0]);
  w = mxGetPr(plhs[1]);
  info = mxGetPr(plhs[2]);



  printf("we enter in solve_cfp using the %s method \n\n", Str);


  if (strcmp(Str, mot1) == 0)
    toto = solve_cfp(vec, qq, &mrows, &meth_cfp, z, w);
  else if (strcmp(Str, mot2) == 0)
    toto = solve_cfp(vec, qq, &mrows, &meth_cfp, z, w);
  else if (strcmp(Str, mot3) == 0)
    toto = solve_cfp(vec, qq, &mrows, &meth_cfp, z, w);
  else printf("Warning : Unknown solving method : %s\n", Str);

  *info = toto;


}












