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
/*!\file main_lcplab.c


   This main subroutine allows the resolution of LCP (Linear Complementary Problem).
   Try \f$(z,w)\f$ such that:

\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z \perp w \ge 0\\
\end{array}
\right.
\f$

  here M is an n by n  matrix, q an n-dimensional vector, w an n-dimensional  vector and z an n-dimensional vector.
  This system of equalities and inequalities is solved thanks to @ref lcp solvers.

  This subroutine is a matlab interface that allows you to use an @ref lcp solver.
  You only have to enter the M matrix, the q vector, the dimension row of q, the maximum of iterations required, the tolerance value, the name of the solver 'Gsnl', 'Lemke', 'Gcp', or 'Latin' and a search direction for the Latin lcp if you want to use it .
\author Nineb Sheherazade.
*/


#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SiconosNumerics_lab.h"
#include "solve_lcp.c"
#include "gsnl_lcp.c"
#include "latin_lcp.c"
#include "lemke_lcp.c"
#include "gcp_lcp.c"



double ddot_(int *, double [], int *, double [], int*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double *vec, *qq, *z, *w, *it_end, *info;
  int mrows, ncols, *nn, n, nl, nc, toto;
  int i, j, *itt;
  double *tol, *itermax, *k_lat;
  double *res;
  char *Str, *mot1, *mot2, *mot3, *mot4;
  methode meth_lcp;


  mot1 = "Lemke";
  mot2 = "Gsnl";
  mot3 = "Latin";
  mot4 = "Gcp";




  if (nlhs != 3)
  {
    mexErrMsgTxt("3 output required.");
  }


  if (nrhs <= 4)
  {
    mexErrMsgTxt("on less 6 input required.");
  }


  vec = mxGetPr(prhs[0]);
  qq = mxGetPr(prhs[1]);
  itermax = mxGetPr(prhs[2]);
  tol = mxGetPr(prhs[3]);





  if (nrhs == 6)
  {
    k_lat = mxGetPr(prhs[5]);
  }
  else *k_lat = 1.;





  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);



  nl = mxGetM(prhs[4]);
  nc = mxGetN(prhs[4]);
  Str = mxCalloc(nl * nc + 1, 8);
  mxGetString(prhs[4], Str, nl * nc + 1);




  meth_lcp.lcp.nom_method = Str;
  meth_lcp.lcp.itermax = *itermax;
  meth_lcp.lcp.k_latin = *k_lat;
  meth_lcp.lcp.tol = *tol;




  plhs[0] = mxCreateDoubleMatrix(mrows, ncols, mxREAL); /* z */
  plhs[1] = mxCreateDoubleMatrix(mrows, ncols, mxREAL); /* w */
  plhs[2] = mxCreateDoubleMatrix(ncols, ncols, mxREAL); /* info */



  z = mxGetPr(plhs[0]);
  w = mxGetPr(plhs[1]);
  info = mxGetPr(plhs[2]);



  printf("we enter in solve_lcp using the %s method \n\n", Str);


  if (strcmp(Str, mot1) == 0)
  {
    toto = solve_lcp(vec, qq, &mrows, &meth_lcp, z, w);
  }
  else if (strcmp(Str, mot2) == 0)
    toto = solve_lcp(vec, qq, &mrows, &meth_lcp, z, w);
  else if (strcmp(Str, mot3) == 0)
    toto = solve_lcp(vec, qq, &mrows, &meth_lcp, z, w);
  else if (strcmp(Str, mot4) == 0)
    toto = solve_lcp(vec, qq, &mrows, &meth_lcp, z, w);
  else printf("Warning : Unknown solving method : %s\n", Str);

  *info = toto;

}












