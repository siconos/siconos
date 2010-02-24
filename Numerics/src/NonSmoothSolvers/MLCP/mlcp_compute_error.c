/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "LA.h"
#include "MLCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
    numericsError("mlcp_compute_error", "null input for problem and/or z and/or w");

  int param = 1;
  int NbLines = problem->M->size0; /* Equalities */
  int n = problem->n; /* Equalities */
  int m = problem->m; /* Inequalities */
  int incx = 1, incy = 1;

  /* Computation of w: depends on the way the problem is written */

  /* Problem in the form (M,q) */
  if (problem->problemType == 0)
  {
    if (problem->M == NULL)
      numericsError("mlcp_compute_error", "null input for M");

    /* Computes w = Mz + q */
    DCOPY(NbLines , problem->q , incx , w , incy);
    prodNumericsMatrix(problem->M->size1, problem->M->size0, 1.0, problem->M, z, 1.0, w);

  }
  /* Problem in the form ABCD */
  else // if(problem->problemType == 1)
  {

    /* Checks inputs */
    if (problem->A == NULL || problem->B == NULL || problem->C == NULL  || problem->D == NULL);
    numericsError("mlcp_compute_error", "null input for A, B,C or D");

    /* Links to problem data */
    double *a = &problem->q[0];
    double *b = &problem->q[NbLines - m];
    double *A = problem->A;
    double *B = problem->B;
    double *C = problem->C;
    double *D = problem->D;

    /* Compute "equalities" part, we = Au + Cv + a - Must be equal to 0 */
    DCOPY(NbLines - m , a , incx , w , incy); //  we = w[0..n-1] <-- a
    DGEMV(LA_NOTRANS , NbLines - m, n , 1.0 , A , NbLines - m , &z[0] , incx , 1.0 , w , incy); // we <-- A*u + we
    DGEMV(LA_NOTRANS , NbLines - m, m , 1.0 , C , NbLines - m , &z[n] , incx , 1.0 , w , incy); // we <-- C*v + we

    /* Computes part which corresponds to complementarity */
    double * pwi = w + NbLines - m; // No copy!!
    DCOPY(m , b , incx , pwi , incy); //  wi = w[n..m] <-- b
    // following int param, we recompute the product wi = Du+BV +b and we = Au+CV +a
    // The test is then more severe if we compute w because it checks that the linear equation is satisfied
    if (param == 1)
    {
      DGEMV(LA_NOTRANS , m, n , 1.0 , D , m , &z[0] , incx , 1.0 , pwi , incy);   // wi <-- D*u+ wi
      DGEMV(LA_NOTRANS , m , m , 1.0 , B , m , &z[n] , incx , 1.0 , pwi , incy);  // wi <-- B*v + wi
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
  if (problem->blocksLine)
  {
    int numBlock = 0;
    while (problem->blocksLine[numBlock] < n + m)
    {
      if (!problem->blocksIsComp[numBlock])
      {
        error_e += DNRM2(problem->blocksLine[numBlock + 1] - problem->blocksLine[numBlock], w + problem->blocksLine[numBlock] , incx);
        norm_e += DNRM2(problem->blocksLine[numBlock + 1] - problem->blocksLine[numBlock], q + problem->blocksLine[numBlock] , incx);
      }
      else
      {
        for (int numLine = problem->blocksLine[numBlock]; numLine < problem->blocksLine[numBlock + 1] ; numLine++)
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
        norm_i += DNRM2(problem->blocksLine[numBlock + 1] - problem->blocksLine[numBlock], w + problem->blocksLine[numBlock] , incx);
      }
      numBlock++;
    }
  }
  else
  {
    printf("WARNING, DEPRECETED MLCP API\n");
    /* Error on equalities part */
    error_e = DNRM2(NbLines - m , w , incx);;

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
    norm_i += DNRM2(m , q + NbLines - m , incx);
    norm_e += DNRM2(NbLines - m , q , incx);
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
