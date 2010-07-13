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
#include "Friction_cst.h"

#define FAC frictionContact3D_GlobalAlartCurnierFunction

#define _00 0
#define _10 1
#define _20 2
#define _01 3
#define _11 4
#define _21 5
#define _02 6
#define _12 7
#define _22 8


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

void printm(char *mes, unsigned int n, unsigned int m, double *mat)
{
  printf("%s\n", mes);

  for (unsigned int j = 0; j < m; ++j)
  {
    for (unsigned int i = 0; i < n; ++i)
    {
      printf(" %10.4g ", mat[i + n * j]);
    }
    printf("\n");
  }
}


void frictionContact3D_LocalAlartCurnierFunction(
  double un,
  double ut1,
  double ut2,
  double rn,
  double rt1,
  double rt2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  assert(result);

  double x0 = rhon * un;
  double x1 = rn - x0;
  double x4 = rhot1 * ut1;
  double x5 = rt1 - x4;
  double x6 = mu * x1;
  double x7 = rhot2 * ut2;
  double x8 = rt2 - x7;
  double x9 = pow(x8, 2);
  double x10 = pow(x5, 2);
  double x11 = x10 + x9;
  double x12 = pow(x11, (1.0 / 2.0));
  double x16 = 1.0 / x12;
  double x19 = pow(x16, 3);
  double x20 = x4 - rt1;
  double x22 = x7 - rt2;
  double x23 = mu * x1 * x16;
  result[6] = 0;
  result[9] = 0;
  result[15] = 0;
  result[18] = 0;
  if (x12 <= x6)
  {
    result[1] = x5 - rt1;
    result[4] = 0;
    result[7] = rhot1;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[19] = 0;
    result[2] = x8 - rt2;
    result[5] = 0;
    result[8] = 0;
    result[11] = rhot2;
    result[14] = 0;
    result[17] = 0;
    result[20] = 0;
  };
  if (x1 < 0)
  {
    result[0] = -rn;
    result[3] = 0;
    result[12] = 1;
  };
  if (x6 < x12)
  {
    result[1] = -rt1 + x23 * x5;
    result[4] = mu * rhon * x16 * x5;
    result[7] = rhot1 * x23 - mu * rhot1 * x1 * x10 * x19;
    result[10] = -mu * rhot2 * x1 * x19 * x5 * x8;
    result[13] = -mu * x16 * x5;
    result[16] = 1 - x23 - mu * x1 * x19 * x20 * x5;
    result[19] = -mu * x1 * x19 * x22 * x5;
    result[2] = -rt2 + x23 * x8;
    result[5] = mu * rhon * x16 * x8;
    result[8] = -mu * rhot1 * x1 * x19 * x5 * x8;
    result[11] = rhot2 * x23 - mu * rhot2 * x1 * x19 * x9;
    result[14] = -mu * x16 * x8;
    result[17] = -mu * x1 * x19 * x20 * x8;
    result[20] = 1 - x23 - mu * x1 * x19 * x22 * x8;
  };
  if (x6 <= 0)
  {
    result[1] = -rt1;
    result[4] = 0;
    result[7] = 0;
    result[10] = 0;
    result[13] = 0;
    result[16] = 1;
    result[19] = 0;
    result[2] = -rt2;
    result[5] = 0;
    result[8] = 0;
    result[11] = 0;
    result[14] = 0;
    result[17] = 0;
    result[20] = 1;
  };
  if (0 <= x1)
  {
    result[0] = -x0;
    result[3] = rhon;
    result[12] = 0;
  };
};

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
  assert(A);
  assert(B);

  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);


  unsigned int i;
  double localf[21];
  double* rn;
  double* rt1;
  double* rt2;
  double* un;
  double* ut1;
  double* ut2;
  double* rhon;
  double* rhot1;
  double* rhot2;
  double* plocalf;


  for (i = 0; i < problemSize; i += 3)
  {
    rn = reaction++;
    rt1 = reaction++;
    rt2 = reaction++;
    un = velocity++;
    ut1 = velocity++;
    ut2 = velocity++;
    rhon = rho++;
    rhot1 = rho++;
    rhot2 = rho++;

    frictionContact3D_LocalAlartCurnierFunction(*un, *ut1, *ut2,
        *rn, *rt1, *rt2,
        *mu++,
        *rhon, *rhot1, *rhot2,
        localf);
    plocalf = localf;


    *result++ = *plocalf++;
    *result++ = *plocalf++;
    *result++ = *plocalf++;

    *A++ = *plocalf++;
    *A++ = *plocalf++;
    *A++ = *plocalf++;

    *A++ = *plocalf++;
    *A++ = *plocalf++;
    *A++ = *plocalf++;

    *A++ = *plocalf++;
    *A++ = *plocalf++;
    *A++ = *plocalf++;

    *B++ = *plocalf++;
    *B++ = *plocalf++;
    *B++ = *plocalf++;

    *B++ = *plocalf++;
    *B++ = *plocalf++;
    *B++ = *plocalf++;

    *B++ = *plocalf++;
    *B++ = *plocalf++;
    *B++ = *plocalf++;

  }

}

void frictionContact3D_GlobalAlartCurnierFunctionO(
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

    // note : rmu was reaction[i0]*mu[i]
    double rmu = mu[i] * fmax(0, result[i0]);

    double p = hypot(result[i1], result[i2]);

    if (result[i0] >= 0)
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
        // DUnPHI2
        A[i00] = 0.;
        B[i00] = 1.;
      }
    }


    if (p > rmu && p > 1e-20)
    {

      // outside disk
      if (A && B)
      {

        double rho1rmu, rho2rmu;
        rho1rmu = rho[i1] * rmu;
        rho2rmu = rho[i2] * rmu;


        GAMMA(result, i1, i2, A, i11, i12, i21, i22, &p);

        int positivepnv = result[i0] > 0;

        // DRnPHI3
        B[i10] = - mu[i] * positivepnv * result[i1] / p; // not rmu
        B[i20] = - mu[i] * positivepnv * result[i2] / p; // not rmu

        // DRtPHI3
        B[i11] = 1. - rmu * A[i11];
        B[i12] = - rmu * A[i12];
        B[i21] = B[i12];
        B[i22] = 1. - rmu * A[i22];

        // DUtPHI3
        A[i11] *= rho1rmu;
        A[i12] *= rho1rmu;
        A[i21] *= rho2rmu;
        A[i22] *= rho2rmu;

      }

      assert(p > 0.);

      result[i1] *= rmu / p;
      result[i2] *= rmu / p;

    }

    else
    {
      // inside disk
      if (A && B)
      {

        // DUtPHI3
        A[i11] = rho[i1];
        A[i12] = 0.;
        A[i21] = 0.;
        A[i22] = rho[i2];

        // DRnPHI3
        B[i10] = 0.;
        B[i20] = 0.;

        // DRtPHI3
        B[i11] = 0.;
        B[i12] = 0.;
        B[i21] = 0.;
        B[i22] = 0.;

      }

    }

    // no convergence :
    //    result[i1] = fmax(0,result[i1]);
    //    result[i2] = fmax(0,result[i2]);

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
};



void zero3x3(double *a)
{
  a[_00] = 0.;
  a[_01] = 0.;
  a[_02] = 0.;
  a[_10] = 0.;
  a[_11] = 0.;
  a[_12] = 0.;
  a[_20] = 0.;
  a[_21] = 0.;
  a[_22] = 0.;
}


/* c = a*b */
/* simple 3x3 multiplication */
/* note: Strassen alg. and others better for size >> 3 */
void multm3x3(double *a, double *b, double *c)
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

void multv3x3(double *a, double *b, double *c)
{

  c[0] = a[_00] * b[0] + a[_01] * b[1] + a[_02] * b[2];
  c[1] = a[_10] * b[0] + a[_11] * b[1] + a[_12] * b[2];
  c[2] = a[_20] * b[0] + a[_21] * b[1] + a[_22] * b[2];

}

/* b += a */
void addm3x3(double *a, double *b)
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
  int k0 = i0 + n * j0;
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
  int k0 = i0 + n * j0;
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
  unsigned int erritermax = options->iparam[1];

  assert(itermax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  unsigned int problemSize2 = problemSize * problemSize;
  unsigned int _3problemSize = 3 * problemSize;

  void *buffer;

  if (!options->dWork)
    buffer = malloc((8 * problemSize +
                     problemSize2) * sizeof(double) +
                    problemSize * sizeof(int));
  else
    buffer = options->dWork;

  double *facWork = (double *) buffer; //malloc(problemSize*sizeof(double));
  double *A = facWork + problemSize; //malloc(3*problemSize*sizeof(double));
  double *B = A + _3problemSize; //malloc(3*problemSize*sizeof(double));
  double *rho = B + _3problemSize; //malloc(problemSize*sizeof(double));
  double *R = rho + problemSize;// malloc(problemSize*problemSize*sizeof(double));
  int *ipiv = (int *)(R + problemSize2);  // malloc(problemSize*sizeof(int));

  double w;
  int dgelsinfo[1];

  DGELS(problemSize, problemSize, 1, R, problemSize, facWork, problemSize, &w, -1, dgelsinfo);

  int LWORK = (int) w;

  double *WORK = (double *) malloc(w * sizeof(double));


  for (unsigned int i = 0; i < problemSize; ++i) rho[i] = 1.;

  info[0] = 1;

  // velocity <- M*reaction + qfree
  DGEMV(LA_NOTRANS, problemSize, problemSize, 1.,
        problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);
  DAXPY(problemSize, 1, problem->q, 1, velocity, 1);



  while (iter++ < itermax)
  {

    frictionContact3D_GlobalAlartCurnierFunction(problemSize, reaction, velocity, rho, problem->mu, facWork, A, B);

    /*    printm("fac",1,3,facWork);
          printm("A",3,3,A);
          printm("B",3,3,B);
    */

    /*    frictionContact3D_GlobalAlartCurnierFunctionO(problemSize, reaction, velocity, rho, problem->mu, facWork, A, B);
        printm("facO",1,3,facWork);
        printm("AO",3,3,A);
        printm("BO",3,3,B);
    */

    // AW + B
    double Wij[9], Aj[9], Bj[9], tmp[9];

    for (unsigned int jp3 = 0, jp9 = 0; jp3 < problemSize; jp3 += 3, jp9 += 9)
    {
      assert(jp9 < 3 * problemSize - 8);

      subextract3x3(3, jp9, 0, A, Aj);
      subextract3x3(3, jp9, 0, B, Bj);

      for (unsigned int ip3 = 0; ip3 < problemSize; ip3 += 3)
      {
        assert(ip3 < problemSize - 2);
        assert(jp3 < problemSize - 2);

        subextract3x3(problemSize, ip3, jp3, problem->M->matrix0, Wij);
        multm3x3(Aj, Wij, tmp);
        if (ip3 == jp3) addm3x3(Bj, tmp);
        subinsert3x3(tmp, problemSize, ip3, jp3, R);

      }

    }

    int fail;

    //DGESV(problemSize, 1, R, problemSize, ipiv, facWork, problemSize, fail );

    DGELS(problemSize, problemSize, 1, R, problemSize, facWork, problemSize, WORK, LWORK, fail);

    assert(fail >= 0);

    if (fail > 0)
      /*if (verbose>0)*/
      printf("GLOBALAC: warning DGESV fail with U(%d,%d) == 0.\n", fail, fail);


    DAXPY(problemSize, 1, facWork, 1, reaction, 1);

    // velocity <- M*reaction + qfree
    DGEMV(LA_NOTRANS, problemSize, problemSize, 1.,
          problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);
    DAXPY(problemSize, 1, problem->q, 1, velocity, 1);


    options->dparam[1] = INFINITY;

    if (!(iter % erritermax))
      FrictionContact3D_compute_error(problem, reaction, velocity,
                                      tolerance, options, &(options->dparam[1]));

    if (verbose > 0)
      printf("GLOBALAC: iteration %d : error=%g\n", iter, options->dparam[1]);

    if (options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }



  if (verbose > 0)
  {
    if (!info[0])
      printf("GLOBALAC: convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    else
      printf("GLOBALAC: no convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
  }

  if (!options->dWork)
  {
    assert(buffer);
    free(buffer);

    if (WORK)
      free(WORK);

  }
  else
  {
    assert(buffer == options->dWork);
  }


}

int frictionContact3D_GlobalAlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the GLOBALAC Solver\n");
  }

  /*  strcpy(options->solverName,"GLOBALAC");*/
  options->solverId = SICONOS_FRICTION_3D_GLOBALAC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *) malloc(options->iSize * sizeof(int));
  options->dparam = (double *) malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (unsigned int i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 200;
  options->iparam[1] = 1;
  options->dparam[0] = 1e-3;

  options->internalSolvers = NULL;

  return 0;
}
