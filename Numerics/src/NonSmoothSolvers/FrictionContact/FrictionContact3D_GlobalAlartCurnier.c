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
#include "FrictionContact3D_Solvers.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_compute_error.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define GAMMA frictionContact3D_gamma
#define FAC frictionContact3D_GlobalAlartCurnierFunction

void frictionContact3D_gamma(double* x, int i1, int i2, double *r,
                             int i11, int i12, int i21, int i22,
                             double *p)
{

  assert(x);
  assert(r);
  // p optional

  double invp_, p_, p3;

  if (p)
  {
    p_ = p[0];
  }
  else
  {
    p_ = hypot(x[i1], x[i2]);
  }

  p3 = p_ * p_ * p_;
  invp_ = 1. / p_;
  r[i11] = invp_ - x[i1] * x[i1] / p3;
  r[i12] = - x[i1] * x[i2] / p3;
  r[i21] = r[i12];
  r[i22] = invp_ - x[i2] * x[i2] / p3;

}



void frictionContact3D_GlobalAlartCurnierFunction(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *rho,
  double *mu,
  double *result,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);
  assert(result);
  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int ip3, ip9, i, i0, i1, i2,
           i00, i01, i02, i10, i11, i12, i20, i21, i22;

  for (i = 0, ip3 = 0, ip9 = 0; ip3 < problemSize; ++i)
  {

    i0 = ip3++;
    i1 = ip3++;
    i2 = ip3++;

    if (A && B)
    {
      // column major
      i00 = ip9++;
      i10 = ip9++;
      i20 = ip9++;
      i01 = ip9++;
      i11 = ip9++;
      i21 = ip9++;
      i02 = ip9++;
      i12 = ip9++;
      i22 = ip9++;

      A[i01] = 0.;
      A[i02] = 0.;
      A[i10] = 0.;
      A[i20] = 0.;

      B[i01] = 0.;
      B[i02] = 0.;

    }


    result[i0] = reaction[i0] - rho[i0] * velocity[i0];
    result[i1] = reaction[i1] - rho[i1] * velocity[i1];
    result[i2] = reaction[i2] - rho[i2] * velocity[i2];

    if (result[i0] > 0)
    {
      result[i0] -= reaction[i0]; // note : this is -PHI p425

      if (A && B)
      {
        // DUnPHI2
        A[i00] = rho[i0];
        B[i00] = 0.;
      }


    }

    else
    {
      result[i0] = -reaction[i0];

      if (A && B)
      {
        A[i00] = 0.;
        B[i00] = 1.;
      }
    }


    double rmu = reaction[i0] * mu[i];
    double p = hypot(result[i1], result[i2]);

    if (p > rmu)
    {

      // outside disk
      if (A && B)
      {
        double mureact, rho1mureact, rho2mureact;
        mureact = mu[i] * reaction[i0];
        rho1mureact = rho[i1] * mureact;
        rho2mureact = rho[i2] * mureact;


        GAMMA(result, i1, i2, A, i11, i12, i21, i22, &p);
        A[i11] *= rho1mureact;
        A[i12] *= rho1mureact;
        A[i21] *= rho2mureact;
        A[i22] *= rho2mureact;

        B[i10] = - mureact * result[i1] / p;
        B[i20] = - mureact * result[i2] / p;

        B[i11] = 1. - A[i11];
        B[i12] = - A[i12];
        B[i21] = - A[i21];
        B[i22] = 1. - A[i22];

      }

      if (rmu <= 0.)
      {
        result[i1] = 0.;
        result[i2] = 0.;
      }
      else
      {
        assert(p > 0.);

        result[i1] *= rmu / p;
        result[i2] *= rmu / p;
      }

    }

    else
    {
      // inside disk
      if (A && B)
      {
        A[i11] = rho[i1];
        A[i12] = 0.;
        A[i21] = 0.;
        A[i22] = rho[i2];

        B[i10] = 0.;
        B[i20] = 0.;

        B[i11] = 0.;
        B[i12] = 0.;
        B[i21] = 0.;
        B[i22] = 0.;

      }

    }

    result[i1] -= reaction[i1];
    result[i2] -= reaction[i2];

    /*
    assert ( ! isnan(result[i0]) );
    assert ( ! isnan(result[i1]) );
    assert ( ! isnan(result[i2]) );

    assert ( ! isnan(A[i00]) );
    assert ( ! isnan(A[i01]) );
    assert ( ! isnan(A[i02]) );
    assert ( ! isnan(A[i10]) );
    assert ( ! isnan(A[i11]) );
    assert ( ! isnan(A[i12]) );
    assert ( ! isnan(A[i20]) );
    assert ( ! isnan(A[i21]) );
    assert ( ! isnan(A[i22]) );


    assert ( ! isnan(B[i00]) );
    assert ( ! isnan(B[i01]) );
    assert ( ! isnan(B[i02]) );
    assert ( ! isnan(B[i10]) );
    assert ( ! isnan(B[i11]) );
    assert ( ! isnan(B[i12]) );
    assert ( ! isnan(B[i20]) );
    assert ( ! isnan(B[i21]) );
    assert ( ! isnan(B[i22]) );

    */
  }



}

#define _00 0
#define _10 1
#define _20 2
#define _01 3
#define _11 4
#define _21 5
#define _02 6
#define _12 7
#define _22 8

/* c = a*b */
/* simple 3x3 multiplication */
/* note: Strassen alg. and others better for size >> 3 */
void mult3x3(double *a, double *b, double *c)
{
  c[_00] = a[_00] * b[_00] + a[_01] * b[_10] + a[_02] * b[_20];
  c[_01] = a[_00] * b[_01] + a[_01] * b[_11] + a[_02] * b[_21];
  c[_02] = a[_00] * b[_02] + a[_01] * b[_12] + a[_02] * b[_22];

  c[_10] = a[_10] * b[_00] + a[_11] * b[_10] + a[_12] * b[_20];
  c[_11] = a[_10] * b[_01] + a[_11] * b[_11] + a[_12] * b[_21];
  c[_12] = a[_10] * b[_02] + a[_11] * b[_12] + a[_12] * b[_22];

  c[_20] = a[_20] * b[_00] + a[_21] * b[_10] + a[_22] * b[_20];
  c[_21] = a[_20] * b[_01] + a[_21] * b[_11] + a[_22] * b[_21];
  c[_22] = a[_20] * b[_02] + a[_21] * b[_12] + a[_22] * b[_22];
}

/* b += a */
void add3x3(double *a, double *b)
{
  b[_00] += a[_00];
  b[_01] += a[_01];
  b[_02] += a[_02];

  b[_10] += a[_10];
  b[_11] += a[_11];
  b[_12] += a[_12];

  b[_20] += a[_20];
  b[_21] += a[_21];
  b[_22] += a[_22];
}


/* copy a sub 3x3 matrix of *a into *b */
/* n row numbers of matrix a
 * i0, j0 indices of first element of the 3x3 submatrix
 * at most 3 cache miss
 */
void subextract3x3(int n, int i0, int j0, double *a, double *b)
{
  int k0 = j0 + n * i0;
  int nm2 = n - 2;


  b[_00] = a[k0];
  b[_10] = a[++k0];
  b[_20] = a[++k0];
  k0 += nm2;
  b[_01] = a[k0];
  b[_11] = a[++k0];
  b[_21] = a[++k0];
  k0 += nm2;
  b[_02] = a[k0];
  b[_12] = a[++k0];
  b[_22] = a[++k0];

  assert((n > 3) ? k0 < n * n : 1);

}


void subinsert3x3(double *a, int n, int i0, int j0, double *b)
{
  int k0 = j0 + n * i0;
  int nm2 = n - 2; // n - 3 + 1

  b[k0]   = a[_00];
  b[++k0] = a[_10];
  b[++k0] = a[_20];
  k0 += nm2;
  b[k0]   = a[_01];
  b[++k0] = a[_11];
  b[++k0] = a[_21];
  k0 += nm2;
  b[k0]   = a[_02];
  b[++k0] = a[_12];
  b[++k0] = a[_22];

  assert(k0 < n * n);

}



void frictionContact3D_GlobalAlartCurnier(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
{
  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);
  assert(problem->M->matrix0);

  unsigned int problemSize = 3 * problem->numberOfContacts;

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  assert(itermax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  /* Check for trivial case */
  *info = checkTrivialCase(problemSize, problem->q,
                           velocity, reaction,
                           options->iparam, options->dparam);

  if (!*info) return;

  /* => all mallocs in SolverOptions dwork + iwork */
  double *facValue = malloc(problemSize * sizeof(double));
  double *A = malloc(3 * problemSize * sizeof(double));
  double *B = malloc(3 * problemSize * sizeof(double));
  double *rho = malloc(problemSize * sizeof(double));
  double *R = malloc(problemSize * problemSize * sizeof(double));
  int *ipiv = malloc(problemSize * sizeof(int));


  for (unsigned int i = 0; i < problemSize; ++i) rho[i] = 1.;

  info[0] = 1;

  while (iter++ < itermax)
  {
    // velocity <- M*reaction + q
    DCOPY(problemSize, problem->q, 1, velocity, 1);
    DGEMV(LA_NOTRANS, problemSize, problemSize, 1.,
          problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);

    FAC(problemSize, reaction, velocity, rho, problem->mu, facValue, A, B);

    // AW + B
    unsigned int jp3 = 0;
    unsigned int jp9 = 0;

    double Wij[9], Aj[9], Bj[9], tmp[9];

    for (unsigned int j = 0; j < problem->numberOfContacts; ++j)
    {
      assert(jp9 < 3 * problemSize - 8);

      subextract3x3(3, 0, jp9, A, Aj);
      subextract3x3(3, 0, jp9, B, Bj);

      unsigned int ip3 = 0;

      for (unsigned int i = 0; i < problem->numberOfContacts; ++i)
      {
        assert(ip3 < problemSize - 2);
        assert(jp3 < problemSize - 2);


        subextract3x3(problemSize, ip3, jp3, problem->M->matrix0, Wij);
        mult3x3(Aj, Wij, tmp);
        add3x3(Bj, tmp);
        subinsert3x3(tmp, problemSize, ip3, jp3, R);

        ip3 += 3;

      }

      jp3 += 3;
      jp9 += 9;

    }

    int fail;
    DGESV(problemSize, 1, R, problemSize, ipiv, facValue, problemSize, fail);
    //    assert ( !fail );

    DAXPY(problemSize, 1, facValue, 1, reaction, 1);

    DCOPY(problemSize, problem->q, 1, velocity, 1);
    DGEMV(LA_NOTRANS, problemSize, problemSize, 1., problem->M->matrix0,
          problemSize, reaction, 1, 1., velocity, 1);

    options->dparam[1] = INFINITY;

    // erriter
    FrictionContact3D_compute_error(problem, reaction, velocity, tolerance, options, &(options->dparam[1]));

    if (verbose > 2)
      printf("error:%g\n", options->dparam[1]);

    if (options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }

  if (verbose > 0)
  {
    if (!info[0])
      printf("convergence after %d iterations, error : %g\n", iter, options->dparam[1]);
    else
      printf("no convergence after %d iterations, error : %g\n", iter, options->dparam[1]);
  }

  // => in SolverOptions
  free(facValue);
  free(A);
  free(B);
  free(rho);
  free(R);
  free(ipiv);

}

int frictionContact3D_GlobalAlartCurnier_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the default solver options for the ACGLOBAL Solver\n");
  }

  strcpy(options->solverName, "ACGLOBAL");

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *) malloc(options->iSize * sizeof(int));
  options->dparam = (double *) malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 20000;
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1e-3;

  options->internalSolvers = NULL;

  return 0;
}
