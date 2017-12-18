
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include "NumericsMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SiconosBlas.h"
#include "numerics_verbose.h"

/*
 * (input) double *z : size n+m
 * (output)double *w : size n+m
 *
 *
 */
int mlcp_compute_error(MixedLinearComplementarityProblem* problem, double *z, double *w, double tolerance, double * error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numerics_error("mlcp_compute_error", "null input for problem and/or z and/or w");

  int param = 1;
  int NbLines = problem->M->size0; /* Equalities */
  int n = problem->n; /* Equalities */
  int m = problem->m; /* Inequalities */
  int incx = 1, incy = 1;

  /* Computation of w: depends on the way the problem is written */

  /* Problem in the form (M,q) */
  if (problem->isStorageType1)
  {
    if (problem->M == NULL)
      numerics_error("mlcp_compute_error", "null input for M");

    /* Computes w = Mz + q */
    cblas_dcopy(NbLines , problem->q , incx , w , incy);
    NM_gemv(1.0, problem->M, z, 1.0, w);

  }
  /* Problem in the form ABCD */
  else //if (problem->isStorageType2)
  {



    /* Checks inputs */
    if (problem->A == NULL || problem->B == NULL || problem->C == NULL  || problem->D == NULL)
    {
      numerics_error("mlcp_compute_error: ", "null input for A, B, C or D");
    }

    /* Links to problem data */
    double *a = &problem->q[0];
    double *b = &problem->q[NbLines - m];
    double *A = problem->A;
    double *B = problem->B;
    double *C = problem->C;
    double *D = problem->D;

    /* Compute "equalities" part, we = Au + Cv + a - Must be equal to 0 */
    cblas_dcopy(NbLines - m , a , incx , w , incy); //  we = w[0..n-1] <-- a
    cblas_dgemv(CblasColMajor,CblasNoTrans , NbLines - m, n , 1.0 , A , NbLines - m , &z[0] , incx , 1.0 , w , incy); // we <-- A*u + we
    cblas_dgemv(CblasColMajor,CblasNoTrans , NbLines - m, m , 1.0 , C , NbLines - m , &z[n] , incx , 1.0 , w , incy); // we <-- C*v + we

    /* Computes part which corresponds to complementarity */
    double * pwi = w + NbLines - m; // No copy!!
    cblas_dcopy(m , b , incx , pwi , incy); //  wi = w[n..m] <-- b
    // following int param, we recompute the product wi = Du+BV +b and we = Au+CV +a
    // The test is then more severe if we compute w because it checks that the linear equation is satisfied
    if (param == 1)
    {
      cblas_dgemv(CblasColMajor,CblasNoTrans , m, n , 1.0 , D , m , &z[0] , incx , 1.0 , pwi , incy);   // wi <-- D*u+ wi
      cblas_dgemv(CblasColMajor,CblasNoTrans , m , m , 1.0 , B , m , &z[n] , incx , 1.0 , pwi , incy);  // wi <-- B*v + wi
    }
  }

  /* Error on equalities part */
  double error_e = 0;
  /* Checks complementarity (only for rows number n to size) */
  double error_i = 0.;
  double zi, wi;
  double *q = problem->q;
  double norm_e = 1;
  double norm_i = 1;
  if (problem->blocksRows)
  {
    int numBlock = 0;
    while (problem->blocksRows[numBlock] < n + m)
    {
      if (!problem->blocksIsComp[numBlock])
      {
        error_e += cblas_dnrm2(problem->blocksRows[numBlock + 1] - problem->blocksRows[numBlock], w + problem->blocksRows[numBlock] , incx);
        norm_e += cblas_dnrm2(problem->blocksRows[numBlock + 1] - problem->blocksRows[numBlock], q + problem->blocksRows[numBlock] , incx);
      }
      else
      {
        for (int numLine = problem->blocksRows[numBlock]; numLine < problem->blocksRows[numBlock + 1] ; numLine++)
        {
          zi = z[numLine];
          wi = w[numLine];
          if (zi < 0.0)
          {
            error_i += -zi;
            if (wi < 0.0) error_i += zi * wi;
          }
          if (wi < 0.0) error_i += -wi;
          if ((zi > 0.0) && (wi > 0.0)) error_i += zi * wi;
        }
        norm_i += cblas_dnrm2(problem->blocksRows[numBlock + 1] - problem->blocksRows[numBlock], w + problem->blocksRows[numBlock] , incx);
      }
      numBlock++;
    }
  }
  else
  {
    printf("WARNING, DEPRECATED MLCP API\n");
    /* Error on equalities part */
    error_e = cblas_dnrm2(NbLines - m , w , incx);;

    /* Checks complementarity (only for rows number n to size) */
    error_i = 0.;

    for (int i = 0 ; i < m ; i++)
    {
      zi = z[n + i];
      wi = w[(NbLines - m) + i];
      if (zi < 0.0)
      {
        error_i += -zi;
        if (wi < 0.0) error_i += zi * wi;
      }
      if (wi < 0.0) error_i += -wi;
      if ((zi > 0.0) && (wi > 0.0)) error_i += zi * wi;
    }


    /* Computes error */
    norm_i += cblas_dnrm2(m , q + NbLines - m , incx);
    norm_e += cblas_dnrm2(NbLines - m , q , incx);
  }

  if (error_i / norm_i >= error_e / norm_e)
  {
    *error = error_i / (1.0 + norm_i);
  }
  else
  {
    *error = error_e / (1.0 + norm_e);
  }

  if (*error > tolerance)
  {
    /*if (isVerbose > 0) printf(" Numerics - mlcp_compute_error failed: error = %g > tolerance = %g.\n",*error, tolerance);*/
    if (verbose)
      printf(" Numerics - mlcp_compute_error failed: error = %g > tolerance = %g.\n", *error, tolerance);
    /* displayMLCP(problem);*/
    return 1;
  }
  else
  {
    if (verbose > 0) printf("Siconos/Numerics: mlcp_compute_error: Error evaluation = %g \n", *error);
    return 0;
  }
}
