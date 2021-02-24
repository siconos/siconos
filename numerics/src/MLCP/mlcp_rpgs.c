/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include <float.h>                              // for DBL_EPSILON
#ifndef __cplusplus
#include <stdbool.h>                       // for false
#endif
#include <stdio.h>                              // for printf
#include <stdlib.h>                             // for free, malloc
#include "MLCP_Solvers.h"                       // for mlcp_compute_error
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "NumericsFwd.h"                        // for SolverOptions, MixedL...
#include "SiconosBlas.h"                        // for cblas_ddot
#include "SolverOptions.h"                      // for SolverOptions, SICONO...
#include "mlcp_cst.h"                           // for SICONOS_DPARAM_MLCP_RHO
#include "numerics_verbose.h"                     // for numerics_printf

#define EPSDIAG DBL_EPSILON

/*
 * double *z : size n+m
 * double *w : size n+m
 */


void mlcp_rpgs(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  double* A = problem->A;
  double* B = problem->B;
  double* C = problem->C;
  double* D = problem->D;
  double* a = problem->a;
  double* b = problem->b;
  int n = problem->n;
  int m = problem->m;
  double *u = &z[0];
  double *v = &z[n];

  int incy, incAx, incBx;
  int i, iter;
  int itermax, verbose;
  double err, vi, viprev, uiprev;
  double tol, rho;
  double *diagA, *diagB;
  verbose = 0;
  incy = 1;
  /* Recup input */

  itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  tol   = options->dparam[SICONOS_DPARAM_TOL];
  rho   = options->dparam[SICONOS_DPARAM_MLCP_RHO];

  /* Initialize output */

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  /* Allocation */

  diagA = (double*)malloc(n * sizeof(double));
  diagB = (double*)malloc(m * sizeof(double));

  incy = 1;

  /* Preparation of the diagonal of the inverse matrix */

  for(i = 0 ; i < n ; ++i)
  {
    if(A[i * n + i] < -EPSDIAG)
    {

      numerics_printf_verbose(1," Vanishing diagonal term A[%i,%i]= %14.8e", i, i,  A[i * n + i] );
      numerics_printf_verbose(1," The local problem cannot be solved");

      *info = 2;
      free(diagA);
      free(diagB);

      return;
    }
    else
    {
      diagA[i] = 1.0 / (A[i * n + i] + rho);

    }
  }
  for(i = 0 ; i < m ; ++i)
  {
    if(B[i * m + i] < -EPSDIAG)
    {

      numerics_printf_verbose(1," Vanishing diagonal term \n");
      numerics_printf_verbose(1," The local problem cannot be solved \n");

      *info = 2;
      free(diagA);
      free(diagB);

      return;
    }
    else
    {
      diagB[i] = 1.0 / (B[i * m + i] + rho);

    }
  }
  /*start iterations*/

  iter = 0;
  err  = 1.;

  incy = 1;
  incAx = n;
  incBx = m;

  mlcp_compute_error(problem, z, w, tol, &err);
  //printf("Error = %12.8e\n", err);

  while((iter < itermax) && (err > tol))
  {

    ++iter;
    incy = 1;


    for(i = 0 ; i < n ; ++i)
    {
      uiprev = u[i];
      u[i] = 0.0;
      //zi = -( q[i] + cblas_ddot( n , &vec[i] , 1 , z , incy ))*diag[i];
      //u[i] = -( a[i]  - (rho*uiprev) +cblas_ddot( n , &A[i] , n , u , 1 )   + cblas_ddot( m , &C[i] , n , v , 1 )         )*diagA[i];
      u[i] = -(a[i]   - (rho * uiprev) + cblas_ddot(n, &A[i], incAx, u, incy)   + cblas_ddot(m, &C[i], incAx, v, incy)) * diagA[i];
    }

    for(i = 0 ; i < m ; ++i)
    {
      viprev = v[i];
      v[i] = 0.0;
      //zi = -( q[i] + cblas_ddot( n , &vec[i] , 1, z , incy ))*diag[i];
      //v[i] = -( b[i] -(rho*viprev) + cblas_ddot( n , &D[i] , m , u , 1 )   + cblas_ddot( m , &B[i] , m , v , 1 )         )*diagB[i];
      vi = -(b[i] - (rho * viprev) + cblas_ddot(n, &D[i], incBx, u, incy)   + cblas_ddot(m, &B[i], incBx, v, incy)) * diagB[i];
      if(vi > 0)
        v[i] = vi;
    }



    /* **** Criterium convergence compliant with filter_result_MLCP **** */
    mlcp_compute_error(problem, z, w, tol, &err);
    numerics_printf_verbose(1,"---- MLCP - RPGS  - Iteration %i rho = %8.4e, residual = %14.7e, tol = %14.7e", rho, iter, err, tol);

    if(verbose == 2)
    {
      for(i = 0 ; i < n ; ++i) printf(" %g", u[i]);
      for(i = 0 ; i < m ; ++i) printf(" %g", v[i]);
      for(i = 0 ; i < m ; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* **** ********************* **** */

  }
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if(err > tol)
  {
    numerics_printf_verbose(1,"---- MLCP - RPGS  - No convergence after %d iterations with error = %14.7e ", iter, err);
    *info = 1;
  }
  else
  {
    numerics_printf_verbose(1,"---- MLCP - RPGS  - Convergence after %d iterations with error = %14.7e ", iter, err);
    *info = 0;
  }
  free(diagA);
  free(diagB);
  return;
}

void mlcp_rpgs_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_IPARAM_MAX_ITER]  = 50000;
  options->dparam[SICONOS_DPARAM_MLCP_RHO] = 1.0; /*rho*/
  options->filterOn = false;
}

