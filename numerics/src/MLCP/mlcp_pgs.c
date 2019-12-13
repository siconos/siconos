/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <math.h>                               // for fabs
#ifndef __cplusplus
#include <stdbool.h>                       // for false
#endif
#include <stdio.h>                              // for printf
#include <stdlib.h>                             // for free, malloc, exit
#include "MLCP_Solvers.h"                       // for mlcp_compute_error
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "NumericsFwd.h"                        // for SolverOptions, MixedL...
#include "SiconosBlas.h"                        // for cblas_ddot, cblas_dcopy
#include "SolverOptions.h"                      // for SolverOptions, SICONO...
#include "mlcp_cst.h"                           // for SICONOS_IPARAM_MLCP_P...


/*
 *
 * double *z : size n+m
 * double *w : size n+m
 */

void mlcp_pgs(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{

  if (!problem->isStorageType2)
  {
    printf("Siconos/Numerics: mlcp_pgs: Wrong Storage (!isStorageType2) for PGS solver\n");
    exit(EXIT_FAILURE);
  }


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
  double *Buf;

  int incx, incy, incAx, incAy, incBx, incBy;
  int i, iter;
  int itermax, verbose;
  int pgsExplicit;
  double err, vi;
  double tol;
  double prev;
  double *diagA, *diagB;
  verbose = 0;

  incx = 1;
  incy = 1;
  /* Recup input */

  itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  pgsExplicit = options->iparam[SICONOS_IPARAM_MLCP_PGS_EXPLICIT];
  tol   = options->dparam[SICONOS_DPARAM_TOL];

  /* Initialize output */

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  /* Allocation */

  diagA = (double*)malloc(n * sizeof(double));
  diagB = (double*)malloc(m * sizeof(double));



  incx = 1;
  incy = 1;

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if ((fabs(A[i * n + i]) < DBL_EPSILON))
    {

      if (verbose > 0)
      {
        printf(" Vanishing diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);
      *info = 1;
      return;
    }
    else
    {
      diagA[i] = 1.0 / A[i * n + i];

    }
  }
  for (i = 0 ; i < m ; ++i)
  {
    if ((fabs(B[i * m + i]) < DBL_EPSILON))
    {

      if (verbose > 0)
      {
        printf(" Vanishing diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);

      return;
    }
    else
    {
      diagB[i] = 1.0 / B[i * m + i];

    }
  }
  /*start iterations*/

  iter = 0;
  err  = 1.;

  incx = 1;
  incy = 1;
  incAx = n;
  incAy = 1;
  incBx = m;
  incBy = 1;


  mlcp_compute_error(problem, z, w, tol, &err);

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    if (pgsExplicit)
    {
      /*Use w like a buffer*/
      cblas_dcopy(n , w , incx , u , incy);  //w <- q
      Buf = w;

      for (i = 0 ; i < n ; ++i)
      {
        prev = Buf[i];
        Buf[i] = 0;
        //zi = -( q[i] + cblas_ddot( n , &vec[i] , incx , z , incy ))*diag[i];
        u[i] =  - (a[i] + cblas_ddot(n , &A[i] , incAx , Buf , incAy)   + cblas_ddot(m , &C[i] , incAx , v , incBy)) * diagA[i];
        Buf[i] = prev;
      }
      for (i = 0 ; i < m ; ++i)
      {
        v[i] = 0.0;
        //zi = -( q[i] + cblas_ddot( n , &vec[i] , incx , z , incy ))*diag[i];
        vi = -(b[i] + cblas_ddot(n , &D[i] , incBx , u , incAy)   + cblas_ddot(m , &B[i] , incBx , v , incBy)) * diagB[i];

        if (vi < 0) v[i] = 0.0;
        else v[i] = vi;
      }
    }
    else
    {

      for (i = 0 ; i < n ; ++i)
      {
        u[i] = 0.0;

        //zi = -( q[i] + cblas_ddot( n , &vec[i] , incx , z , incy ))*diag[i];
        u[i] =  - (a[i] + cblas_ddot(n , &A[i] , incAx , u , incAy)   + cblas_ddot(m , &C[i] , incAx , v , incBy)) * diagA[i];
      }

      for (i = 0 ; i < m ; ++i)
      {
        v[i] = 0.0;
        //zi = -( q[i] + cblas_ddot( n , &vec[i] , incx , z , incy ))*diag[i];
        vi = -(b[i] + cblas_ddot(n , &D[i] , incBx , u , incAy)   + cblas_ddot(m , &B[i] , incBx , v , incBy)) * diagB[i];

        if (vi < 0) v[i] = 0.0;
        else v[i] = vi;
      }
    }

    /* **** Criterium convergence compliant with filter_result_MLCP **** */

    mlcp_compute_error(problem, z, w, tol, &err);

    if (verbose == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", u[i]);
      for (i = 0 ; i < m ; ++i) printf(" %g", v[i]);
      for (i = 0 ; i < m ; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* **** ********************* **** */

  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if (err > tol)
  {
    printf("Siconos/Numerics: mlcp_pgs: No convergence of PGS after %d iterations\n" , iter);
    printf("Siconos/Numerics: mlcp_pgs: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: mlcp_pgs: Convergence of PGS after %d iterations\n" , iter);
      printf("Siconos/Numerics: mlcp_pgs: The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(diagA);
  free(diagB);
  return;
}

void mlcp_pgs_set_default(SolverOptions* options)
{
  options->filterOn = false;
  options->iparam[SICONOS_IPARAM_MLCP_PGS_EXPLICIT] = 0; //implicit
}

