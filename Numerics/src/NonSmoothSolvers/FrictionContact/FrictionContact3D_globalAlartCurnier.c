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

  double* r11 = r + _11;
  double* r12 = r + _12;
  double* r21 = r + _21;
  double* r22 = r + _22;

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


  double x0;
  double x1;
  double x2;
  double x3;
  double x4;
  double x5;
  double x6;
  double x7;
  double x8;
  double x9;
  double x10;
  double x11;
  double x12;
  double x13;
  double x14;
  double x15;
  double x16;
  double x17;
  double x18;
  double x19;
  double x20;
  double x21;
  double x22;
  double x23;
  double x24;
  x0 = rhon * un;
  x1 = rn - x0;
  x2 = x1 < 0;
  x3 = 0 <= x1;
  x4 = rhot1 * ut1;
  x5 = rt1 - x4;
  x6 = rhot2 * ut2;
  x7 = rt2 - x6;
  x8 = pow(x7, 2);
  x9 = pow(x5, 2);
  x10 = x8 + x9;
  x11 = mu * x1;
  x12 = pow(x10, (1.0 / 2.0));
  x13 = x12 <= x11;
  x14 = 1.0 / x12;
  x15 = x11 < x12;
  x16 = 0 < x11;
  x17 = x11 <= 0;
  x18 = (0, x11 <= 0);
  x19 = (0, x12 <= x11);
  x20 = pow(x14, 3);
  x21 = x4 - rt1;
  x22 = (1, x12 <= x11);
  x23 = x6 - rt2;
  x24 = x11 * x14;
  result[6] = 0;
  result[9] = 0;
  result[15] = 0;
  result[18] = 0;
  if ((0 < x11) && (x11 < x12))
  {
    result[1] = -rt1 + x24 * x5;
    result[4] = mu * rhon * x14 * x5;
    result[7] = rhot1 * x24 - rhot1 * x11 * x20 * x9;
    result[10] = -rhot2 * x11 * x20 * x5 * x7;
    result[13] = -mu * x14 * x5;
    result[16] = 1 - x24 - mu * x1 * x20 * x21 * x5;
    result[19] = -mu * x1 * x20 * x23 * x5;
    result[2] = -rt2 + x24 * x7;
    result[5] = mu * rhon * x14 * x7;
    result[8] = -rhot1 * x11 * x20 * x5 * x7;
    result[11] = rhot2 * x24 - rhot2 * x11 * x20 * x8;
    result[14] = -mu * x14 * x7;
    result[17] = -mu * x1 * x20 * x21 * x7;
    result[20] = 1 - x24 - mu * x1 * x20 * x23 * x7;
  };
  if (0 <= x1)
  {
    result[0] = -x0;
    result[3] = rhon;
    result[12] = 0;
  };
  if ((0 < x11) && (x12 <= x11))
  {
    result[1] = x5 - rt1;
    result[4] = 0;
    result[7] = rhot1;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[19] = 0;
    result[2] = x7 - rt2;
    result[5] = 0;
    result[8] = 0;
    result[11] = rhot2;
    result[14] = 0;
    result[17] = 0;
    result[20] = 0;
  };
  if (x11 <= 0)
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
  if (x1 < 0)
  {
    result[0] = -rn;
    result[3] = 0;
    result[12] = 1;
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
  *B01 = 0.;
  *B02 = 0.;


  double D0, D1, D2, muD0;

  D0 = *reaction0 - *rho0 * *velocity0;
  D1 = *reaction1 - *rho1 * *velocity1;
  D2 = *reaction2 - *rho2 * *velocity2;

  muD0 = mu * D0;

  double hypotD1D2 = hypot(D1, D2);


  if (muD0 <= 0.)
  {


    *ACF1 = *reaction1;
    *ACF2 = *reaction2;
    *A10 = 0.;
    *A11 = 0.;
    *A12 = 0.;
    *A20 = 0.;
    *A21 = 0.;
    *A22 = 0.;
    *B10 = 0.;
    *B11 = 1.;
    *B12 = 0.;
    *B20 = 0.;
    *B22 = 1.;

  };

  if (hypotD1D2 <= muD0)
  {
    *ACF1 = *reaction1 - D1;
    *ACF2 = *reaction2 - D2;
    *A10 = 0.;
    *A11 = *rho1;
    *A12 = 0.;
    *A20 = 0.;
    *A21 = 0.;
    *A22 = *rho2;

    *B10 = 0.;
    *B11 = 0.;
    *B12 = 0.;
  }

  if (D0 < 0.)
  {
    *ACF0 = *reaction0;
    *A00 = 0.;
    *B00 = 1.;
  }

  if (D0 >= 0.)
  {
    *ACF0 = *reaction0 - D0;
    *A00 = *rho0;
    *B00 = 0.;
  }

  if (muD0 < hypotD1D2)
  {
    double cubehypotD1D2 = hypotD1D2 * hypotD1D2 * hypotD1D2;
    double muD0rho1 = muD0* *rho1;

    *ACF1 = *reaction1 - muD0 * D1 / hypotD1D2;
    *ACF2 = *reaction2 - muD0 * D2 / hypotD1D2;

    *A10 = mu * D1* *rho0 / hypotD1D2;
    *A11 = -muD0rho1 * D1 * D1 / cubehypotD1D2 + muD0rho1 / hypotD1D2;
    *A12 = - D0 * D1 * D2 * mu**rho2 / cubehypotD1D2;

    *A20 = D2 * mu* *rho0 / hypotD1D2;
    *A21 = - D0 * D1 * D2 * mu* *rho1 / cubehypotD1D2;
    *A22 = - muD0* *rho2 * D2 * D2 / cubehypotD1D2 + muD0* *rho2 / hypotD1D2;

    *B10 = - mu * D1 / hypotD1D2;
    *B11 = 1 + muD0 * D1 * D1 / cubehypotD1D2 - muD0 / hypotD1D2;
    *B12 =  muD0 * D1 * D2 / cubehypotD1D2;

    *B20 = - mu * D2 / hypotD1D2;
    *B21 =  muD0 * D1 * D2 / cubehypotD1D2;
    *B22 = 1 + muD0 * D2 * D2 / cubehypotD1D2 - muD0 / hypotD1D2;

  }

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
