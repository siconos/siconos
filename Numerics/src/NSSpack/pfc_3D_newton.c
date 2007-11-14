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
 * \fn  pfc_3D_newton( int n , double *C , double *b ,  double *zz ,double *ww , double mu ,
 *        Compute_G_function (*Compute_G), Compute_JacG_function (*Compute_JacG),
 *        int *iparam_local , double *dparam_local )
 *
 *
 * \author houari khenous last last modification 07/11/2007.
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



void pfc_3D_newton(int n , double *C , double *b ,  double *zz , double *ww , double mu ,
                   pfc3D_fPtr(*Compute_G), pfc3D_fPtr(*Compute_JacG),
                   double *param1, double *param2, double *param3, int *iparam_local , double *dparam_local)
{

  int i, j, niter, mm, local_itermax;
  double nerr, nerr1, an, at, local_tol;
  double a1, qs, alpha, beta, det;
  int incx, incy;
  double *www, *G, *JacG, *AA , *wwww, *zzzz;
  int nrhs = 1, infoDGESV;
  int *ipiv;
  ipiv = (int *)malloc(n * sizeof(int));

  Linesearch_function Linesearch;

  if (iparam_local[3] == 0)
    Linesearch = &Linesearch_AC;
  else
    Linesearch = &Linesearch_FB;


  local_itermax = iparam_local[0];
  local_tol     = dparam_local[0];



  mm   = n * n;
  incx =  1;
  incy =  1;
  a1 = -1.;
  qs = 1.;

  JacG = (double*)malloc(mm * sizeof(double));
  AA   = (double*)malloc(mm * sizeof(double));

  G    = (double*)malloc(n * sizeof(double));
  www  = (double*)malloc(n * sizeof(double));

  wwww = (double*)malloc(n * sizeof(double));
  zzzz = (double*)malloc(n * sizeof(double));

  /* Intialization of G, JacG, ww and www */
  for (i = 0 ; i < n ; ++i)
  {
    G[i] = www[i] = 0.;
    for (j = 0 ; j < n ; ++j)
      JacG[j * n + i] = AA[j * n + i] = 0.;
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

  while ((niter < local_itermax) && (nerr > local_tol))
  {
    ++niter;
    /*     printf("-----------------------------------Iteration Newton %i \n",niter); */
    (*Compute_G)(n , G , zz , C , ww , b , param1, param2, param3 , an , at , mu);

    (*Compute_JacG)(n , JacG , zz , C , ww , b , param1, param2, param3 , an , at , mu);

    nerr1 = DNRM2(n, G , incx);

    /***** Criterium convergence *****/
    /* compute the direction www s.t X^{k+1} = X^{k} + www, where X^{k} = za */

    incx =  1;
    incy =  1;
    a1   = -1.0;
    qs   = 1.0;

    /* compute the direction www s.t X^{k+1} = X^{k} + www, where X^{k} = ww = (wa,za) */

    DCOPY(n , G , incx , www , incy);
    DSCAL(n , a1 , www, incx);
    DCOPY(mm , JacG , incx , AA , incy);


    DGESV(n, nrhs, AA, n, ipiv, www, n, infoDGESV);

    /*    printf("INFO =  %i \n",infoDGESV); */

    if (infoDGESV)
    {
      /*  printf("Problem in DGESV\n");   */

      free(AA);
      free(G);
      free(JacG);
      free(www);
      free(wwww);
      free(zzzz);
      free(ipiv);

      return;

    }

    /*   /\* direct inverse of 3x3 matrix *\/ */
    /*     matrix_inv3( JacG, AA); */
    /*     DGEMV( LA_NOTRANS , n , n , a1 , AA , n , G , incx , qs , www , incy ); */

    (*Linesearch)(n , zz , ww , www , b , C , param1, param2, param3 , an , at , mu , nerr1);
    nerr = nerr1;

    for (i = 0 ; i < n ; ++i)
      www[i] = 0.;
    /*    printf("-----------------------------------Iteration Newton %i --------- Newton Error = %14.7e\n",niter,nerr); */
  }


  free(AA);
  free(G);
  free(JacG);
  free(www);
  free(wwww);
  free(zzzz);
  free(ipiv);
}
