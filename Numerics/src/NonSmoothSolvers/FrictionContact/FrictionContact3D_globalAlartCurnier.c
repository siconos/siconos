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
#include "op3x3.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_compute_error.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"


void frictionContact3D_gamma(double* x, double *r, double *p)
{
  assert(x);
  assert(r);
  // p optional

  double* x1 = ++x;
  double* x2 = ++x;

  //  double* r11 = r+_11;
  //double* r12 = r+_12;
  //double* r21 = r+_21;
  //double* r22 = r+_22;
  SET3X3(r);

  double invp_, p_, p3;

  if (p)
  {
    p_ = p[0];
  }
  else
  {
    p_ = hypot(*x1, *x2);
  }

  p3 = p_ * p_ * p_;
  invp_ = 1. / p_;
  *r11 = invp_ - *x1 * *x1 / p3;
  *r12 = - *x1 * *x2 / p3;
  *r21 = *r12;
  *r22 = invp_ - *x2 * *x2 / p3;

}

void frictionContact3D_localAlartCurnierFunctionGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
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

void frictionContact3D_localAlartCurnierFunction(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{
  double result[21]; //3 + 2 * 9

  SET3(reaction);
  SET3(velocity);
  SET3(rho);


  frictionContact3D_localAlartCurnierFunctionGenerated(
    *reaction0, *reaction1, *reaction2,
    *velocity0, *velocity1, *velocity2,
    mu,
    *rho0, *rho1, *rho2,
    result);

  cpy3(result, f);
  cpy3x3(result + 3, A);
  cpy3x3(result + 12, B);
}


void frictionContact3D_globalAlartCurnierFunction(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
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
  for (i = 0; i < problemSize; i += 3)
  {

    frictionContact3D_localAlartCurnierFunction(reaction,
        velocity,
        *mu,
        rho,
        result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;
    result += 3;
    A += 9;
    B += 9;

  }

}

void frictionContact3D_localAlartCurnierFunctionHandMade(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *ACF,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);
  assert(ACF);
  assert(A);
  assert(B);

  SET3X3(A);
  SET3X3(B);
  SET3(ACF);
  SET3(reaction);
  SET3(velocity);
  SET3(rho);



  *A01 = 0.;
  *A02 = 0.;
  *A10 = 0.;
  *A20 = 0.;

  *B01 = 0.;
  *B02 = 0.;


  double p0, p1, p2;

  p0 = *reaction0 - *rho0 * *velocity0;
  p1 = *reaction1 - *rho1 * *velocity1;
  p2 = *reaction2 - *rho2 * *velocity2;

  // note : rmu was reaction[i0]*mu[i]
  double rmu = mu * fmax(0, p0);

  double p = hypot(p1, p2);

  if (p0 >= 0)
  {
    *ACF0 -= *reaction0; // note : this is -PHI p425

    if (A && B)
    {
      // DUnPHI2
      *A00 = *rho0;
      *B00 = 0.;
    }


  }

  else
  {
    *ACF0 = -*reaction0;

    if (A && B)
    {
      // DUnPHI2
      *A00 = 0.;
      *B00 = 1.;
    }
  }


  if (p > rmu)
  {

    // outside disk
    if (A && B)
    {

      double rho1rmu, rho2rmu;
      rho1rmu = *rho1 * rmu;
      rho2rmu = *rho2 * rmu;

      frictionContact3D_gamma(ACF0, A00, &p);

      int positivepnv = *ACF0 > 0;

      // DRnPHI3
      *B10 = - mu * positivepnv * p1 / p; // not rmu
      *B20 = - mu * positivepnv * p2 / p; // not rmu

      // DRtPHI3
      *B11 = 1. - rmu* *A11;
      *B12 = - rmu* *A12;
      *B21 = *B12;
      *B22 = 1. - rmu* *A22;

      // DUtPHI3
      *A11 *= rho1rmu;
      *A12 *= rho1rmu;
      *A21 *= rho2rmu;
      *A22 *= rho2rmu;

    }

    assert(p > 0.);

    *ACF1 *= rmu / p;
    *ACF2 *= rmu / p;

  }

  else
  {
    // inside disk
    if (A && B)
    {

      // DUtPHI3
      *A11 = *rho1;
      *A12 = 0.;
      *A21 = 0.;
      *A22 = *rho2;

      // DRnPHI3
      *B10 = 0.;
      *B20 = 0.;

      // DRtPHI3
      *B11 = 0.;
      *B12 = 0.;
      *B21 = 0.;
      *B22 = 0.;

    }

  }

  *ACF1 -= *reaction1;
  *ACF2 -= *reaction2;
};


void frictionContact3D_globalAlartCurnierFunctionHandMade(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
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
  for (i = 0; i < problemSize; i += 3)
  {

    frictionContact3D_localAlartCurnierFunctionHandMade(reaction,
        velocity,
        *mu,
        rho,
        result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;
    result += 3;
    A += 9;
    B += 9;

  }

}



void frictionContact3D_globalAlartCurnier(
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

  DGELS(problemSize, problemSize,
        1, R, problemSize,
        facWork, problemSize, &w, -1, dgelsinfo);

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

    frictionContact3D_globalAlartCurnierFunction(problemSize,
        reaction, velocity,
        problem->mu, rho,
        facWork, A, B);

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

      extract3x3(3, jp9, 0, A, Aj);
      extract3x3(3, jp9, 0, B, Bj);

      for (unsigned int ip3 = 0; ip3 < problemSize; ip3 += 3)
      {
        assert(ip3 < problemSize - 2);
        assert(jp3 < problemSize - 2);

        extract3x3(problemSize, ip3, jp3, problem->M->matrix0, Wij);
        mm3x3(Aj, Wij, tmp);
        if (ip3 == jp3) add3x3(Bj, tmp);
        insert3x3(problemSize, ip3, jp3, R, tmp);

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

int frictionContact3D_globalAlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the GLOBALAC Solver\n");
  }

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
