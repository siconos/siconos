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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif

int lcp_driver(double *vec, double *q , int *n , method *pt , double *z , double *w)
{

  const char lcpkey1[10] = "Lemke", lcpkey2[10] = "PGS", lcpkey3[10] = "CPG";
  const char lcpkey4[10] = "Latin", lcpkey5[10] = "QP", lcpkey6[10] = "NSQP";
  const char lcpkey8[15] = "NewtonMin";
  const char lcpkey9[15] = "Latin_w", lcpkey10[15] = "NewtonFB", lcpkey11[15] = "PSOR";
  const char lcpkey13[10] = "RPGS";
  const char lcpkey14[10] = "Path";

  /* Output result
     0: ok
     >0: problem (depends on solver)
  */
  int info = 1;

  /******************************************
   *  1 - Check for trivial solution
   ******************************************/

  /*  limqpos = -1e-16 / sqrt((double) *n); */
  int i = 0;
  while ((i < (*n - 1)) && (q[i] >= 0.)) i++;
  if ((i == (*n - 1)) && (q[*n - 1] >= 0.))
  {
    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (int j = 0 ; j < *n; j++)
    {
      z[j] = 0.0;
      w[j] = q[j];
    }
    info = 0;
    pt->lcp.iter = 0;
    pt->lcp.err  = 0.;
    if (pt->lcp.chat > 0) printf("LCP_driver: trivial solution for LCP (positive vector q)  \n");
    return info;
  }

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Lists of parameters for solvers
     iparamLCP[0] = maximum iterations number
     iparamLCP[1] = verbose mode (0:off, >0:on)
     iparamLCP[2] = number of iterations done (output)

     dparamLCP[0] = tolerance
     dparamLCP[1] = depends on solver
     dparamLCP[2] = depends on solver
     dparamLCP[3] = depends on solver
  */
  int     iparamLCP[3];
  double  dparamLCP[4];
  iparamLCP[1] = pt->lcp.chat;

  /****** Lemke algorithm ******/
  /* IN: itermax
     OUT: iter */
  if (strcmp(pt->lcp.name , lcpkey1) == 0)
  {
    iparamLCP[0] = pt->lcp.itermax;
    lcp_lexicolemke(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
  }

  /****** PGS Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else if (strcmp(pt->lcp.name , lcpkey2) == 0)
  {
    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;

    lcp_pgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];
  }

  /****** CPG Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else if (strcmp(pt->lcp.name , lcpkey3) == 0)
  {
    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;

    lcp_cpg(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];
  }

  /****** Latin Solver ******/
  /* IN: itermax, tolerance, k_latin
     OUT: iter, error */
  else if (strcmp(pt->lcp.name , lcpkey4) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;

    lcp_latin(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }

  /****** Latin_w Solver ******/
  /* IN: itermax, tolerance, k_latin, relax
     OUT: iter, error */
  else if (strcmp(pt->lcp.name , lcpkey9) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;
    dparamLCP[3] = pt->lcp.relax;


    lcp_latin_w(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }

  /****** QP Solver ******/
  /* IN: tolerance
     OUT:
  */
  else if (strcmp(pt->lcp.name , lcpkey5) == 0)
  {
    /* We assume that the LCP matrix M is symmetric*/

    dparamLCP[0] = pt->lcp.tol;

    lcp_qp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }

  /****** NSQP Solver ******/
  /* IN: tolerance
     OUT:
  */
  else if (strcmp(pt->lcp.name , lcpkey6) == 0)
  {
    /* We assume that the LCP matrix M is not symmetric*/

    dparamLCP[0] = pt->lcp.tol;

    lcp_nsqp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);
  }

  /****** Newton min ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else if (strcmp(pt->lcp.name , lcpkey8) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;

    lcp_newton_min(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];
  }

  /****** Newton Fischer-Burmeister ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else if (strcmp(pt->lcp.name , lcpkey10) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;

    lcp_newton_FB(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];
  }

  /****** PSOR Solver ******/
  /* IN: itermax, tolerance, relax
     OUT: iter, error
  */
  else if (strcmp(pt->lcp.name , lcpkey11) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.relax;

    lcp_psor(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];
  }

  /****** RPGS (Regularized Projected Gauss-Seidel) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else if (strcmp(pt->lcp.name , lcpkey13) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.rho;
    /* dparamLCP[2] = pt->lcp.relax;*/

    lcp_rpgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[3];
  }


  /****** PATH (Ferris) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else if (strcmp(pt->lcp.name , lcpkey14) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.rho;

    lcp_path(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[3];
  }

  else
    printf("LCP_driver error: unknown solver named: %s\n", pt->lcp.name);

  /*************************************************
   *  3 - Check solution validity
   *************************************************/

  /* Warning: it depends on the chosen solver */

  /* Not done for:  PGS, RPGS */
  if ((strcmp(pt->lcp.name , lcpkey2) != 0) && (strcmp(pt->lcp.name , lcpkey13) != 0))
    info = filter_result_LCP(*n, vec, q, z, pt->lcp.tol, pt->lcp.chat, w);


  return info;

}
