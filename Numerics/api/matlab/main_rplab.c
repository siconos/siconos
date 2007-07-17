/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
/*!\file main_rplab.c


 This subroutine allows the primal resolution of relay problems.

  Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
-w \in \partial\psi_{[b, a]}(z)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 This system of equations and inequalities is solved thanks to @ref pr solvers.
 The routine's call is due to the function solve_rp.c.

  This subroutine is a matlab interface that allows you to use an @ref pr solver.
  You only have to enter the M matrix, the q vector, the boundary a and b, the maximum of iterations required, the tolerance value, the name of the solver 'Gsnl', or 'Latin', and a search direction for the Latin  if you wanna use it .
\author Nineb Sheherazade.

*/


#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SiconosNumerics_lab.h"
#include "solve_rp.c"
#include "rp_gsnl.c"
#include "rp_latin.c"


double ddot_(int *, double [], int *, double [], int*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double *vec, *qq, *q, *z, *zt, *w, *it_end, *info, *a, *b, *c, *cc;
  int mrows, ncols, n, nl, nc;
  int i, j, *itt, incx = 1, incy = 1, toto;
  double *tol, *itermax, *k_lat, *mu, alpha, beta;
  double *res;
  char *Str, *mot1, *mot2, *mot3, trans = 'T';
  methode meth_rp;



  mot1 = "Gsnl";
  mot2 = "Latin";
  mot3 = "Gcp";



  if (nlhs != 3)
  {
    mexErrMsgTxt("3 output required.");
  }

  if (nrhs <= 6)
  {
    mexErrMsgTxt("on less 7 input required.");
  }

  vec = mxGetPr(prhs[0]);
  q = mxGetPr(prhs[1]);
  a = mxGetPr(prhs[2]);
  b = mxGetPr(prhs[3]);
  itermax = mxGetPr(prhs[4]);
  tol = mxGetPr(prhs[5]);



  if (nrhs == 8)
  {
    k_lat = mxGetPr(prhs[7]);
  }
  else *k_lat = 1.;




  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);


  nl = mxGetM(prhs[6]);
  nc = mxGetN(prhs[6]);
  Str = mxCalloc(nl * nc + 1, 8);
  mxGetString(prhs[6], Str, nl * nc + 1);




  meth_rp.rp.nom_method = Str;
  meth_rp.rp.itermax = *itermax;
  meth_rp.rp.k_latin = *k_lat;
  meth_rp.rp.tol = *tol;
  meth_rp.rp.a = (double*)malloc(mrows * sizeof(double));
  meth_rp.rp.b = (double*)malloc(mrows * sizeof(double));

  for (i = 0; i <= mrows - 1; i++)
  {
    meth_rp.rp.a[i] = a[i] ;
    meth_rp.rp.b[i] = -b[i] ;
  }




  plhs[0] = mxCreateDoubleMatrix(mrows, ncols, mxREAL); /* z */
  plhs[1] = mxCreateDoubleMatrix(mrows, ncols, mxREAL); /* w */
  plhs[2] = mxCreateDoubleMatrix(ncols, ncols, mxREAL); /* info */



  z = mxGetPr(plhs[0]);
  w = mxGetPr(plhs[1]);
  info = mxGetPr(plhs[2]);



  printf("we enter in solve_rp using the %s method \n\n", Str);


  if (strcmp(Str, mot1) == 0)
    toto = solve_rp(vec, q, &mrows, &meth_rp, z, w);
  else if (strcmp(Str, mot2) == 0)
    toto = solve_rp(vec, q, &mrows, &meth_rp, z, w);
  else if (strcmp(Str, mot3) == 0)
    toto = solve_rp(vec, q, &mrows, &meth_rp, z, w);
  else printf("Warning : Unknown solving method : %s\n", Str);

  *info = toto;


  free(meth_rp.rp.a);
  free(meth_rp.rp.b);

}












