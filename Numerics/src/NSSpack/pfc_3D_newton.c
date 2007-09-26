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
/*!\file pfc_3D_newton.c
 *
 * \fn  pfc_3D_newton( int n , double *C , double *b , double *zz , \n
 *                     double *ww , double mu , int itermax , double tol )
 *
 * Generic pfc_3D parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param COEF    Unchanged parameter which represents the friction coefficient
 *
 *
 * \author houari khenous 14/09/2007
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>
#include <NSSpack.h>
#include <pfc_3D_Alart_Curnier.h>
#include <pfc_3D_Fischer_Burmeister.h>



void pfc_3D_newton(int n , double *C , double *b , double *zz , double *ww , double mu , Compute_G_function(*Compute_G), Compute_JacG_function(*Compute_JacG), int *iparam_local , double *dparam_local)
{

  int i, j, niter, mm;
  double nerr, nerr1, an, at;
  double a1, qs, alpha, beta, det;
  int incx, incy;
  double *www, *G, *JacG, *AA , *A, *B, *AC, *wwww, *zzzz;

  Linesearch_function Linesearch;

  int local_formulation = 0;

  if (local_formulation == 0)
    Linesearch = &Linesearch_AC;
  else
    Linesearch = &Linesearch_FB;

  mm   = n * n;
  incx =  1;
  incy =  1;
  a1 = -1.;
  qs = 1.;

  JacG = (double*)malloc(mm * sizeof(double));
  AA   = (double*)malloc(mm * sizeof(double));
  AC   = (double*)malloc(mm * sizeof(double));
  A    = (double*)malloc(mm * sizeof(double));
  B    = (double*)malloc(mm * sizeof(double));

  G    = (double*)malloc(n * sizeof(double));
  www  = (double*)malloc(n * sizeof(double));

  wwww = (double*)malloc(n * sizeof(double));
  zzzz = (double*)malloc(n * sizeof(double));

  /* Intialization of G, JacG, ww and www */
  for (i = 0 ; i < n ; ++i)
  {
    G[i] = www[i] = 0.;
    for (j = 0 ; j < n ; ++j)
      JacG[j * n + i] = AA[j * n + i] =  A[j * n + i] = B[j * n + i] = 0.;
  }

  an = 1. / C[0 * n + 0];
  alpha = C[1 * n + 1] + C[2 * n + 2];
  det = C[1 * n + 1] * C[2 * n + 2] - C[2 * n + 1] + C[1 * n + 2];
  beta = alpha * alpha - 4 * det;
  if (beta > 0.)
    beta = sqrt(beta);
  else
    beta = 0.;

  at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

  niter = 0;
  nerr  = 1.;

  /* Newton loop */

  while ((niter < 2000) && (nerr > 1.e-3))
  {
    ++niter;

    (*Compute_G)(n , G , C , zz , ww , b , A , B , AC , an , at , mu);

    (*Compute_JacG)(n , JacG , C , zz , ww , b , A , B , an , at , mu);

    nerr1 = DNRM2(n, G , incx);

    /***** Criterium convergence *****/
    /* compute the direction www s.t X^{k+1} = X^{k} + www, where X^{k} = za */
    for (i = 0 ; i < n ; ++i)
      www[i] = 0.;
    //matrix_inv( n,JacG, AA);
    matrix_inv3(JacG, AA);

    DGEMV(LA_NOTRANS , n , n , a1 , AA , n , G , incx , qs , www , incy);

    (*Linesearch)(n , G , zz , ww , www , b , C , zzzz, wwww , an , at , mu , nerr1);
    nerr = nerr1;
  }


  free(AA);
  free(A);
  free(B);
  free(AC);
  free(G);
  free(JacG);
  free(www);
  free(wwww);
  free(zzzz);
}
