/* Siconos-Numerics, Copyright INRIA 2005-2012.
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


#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_compute_error.h"
#include "FischerBurmeisterGenerated.h"
#include "FrictionContactNonsmoothEqn.h"
#include "FrictionContact3D_localFischerBurmeister.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "SiconosLapack.h"


void frictionContact3D_FischerBurmeisterFunction(
  unsigned int problemSize,
  FischerBurmeisterFun3x3Ptr computeACFun3x3,
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

  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for (i = 0; i < problemSize; i += 3)
  {

    computeACFun3x3(reaction, velocity, *mu, rho, result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;

    if (result)
      result += 3;

    if (A)
      A += 9;

    if (B)
      B += 9;

  }

}



void frictionContact3D_localFischerBurmeister(
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

  FischerBurmeisterFun3x3Ptr computeACFun3x3;

  switch (options->iparam[10])
  {
  case 0:
  {

    computeACFun3x3 = &frictionContact3D_FischerBurmeisterFunctionGenerated;
    break;
  }
  }


  if (!problem->M->matrix0)
  {
    frictionContact3D_sparseLocalFischerBurmeister(
      problem,
      reaction,
      velocity,
      info,
      options);
    return;
  }

  assert(problem->M->matrix0);

  unsigned int problemSize = 3 * problem->numberOfContacts;

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[7];

  assert(itermax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  unsigned int problemSize2 = problemSize * problemSize;
  unsigned int _3problemSize = 3 * problemSize;

  void *buffer;

  if (!options->dWork)
  {
#ifndef NDEBUG
    buffer = malloc((14 * problemSize +
                     2 * problemSize2) * sizeof(double) +
                    problemSize * sizeof(int));
#else
    buffer = malloc((14 * problemSize +
                     problemSize2) * sizeof(double) +
                    problemSize * sizeof(int));
#endif
  }
  else
    buffer = options->dWork;

  double *F = (double *) buffer; //malloc(problemSize*sizeof(double));
  double *tmp1 = (double *) F + problemSize; //malloc(problemSize*sizeof(double));
  double *tmp2 = (double *) tmp1 + problemSize; //malloc(problemSize*sizeof(double));
  double *A = tmp2 + problemSize; //malloc(3*problemSize*sizeof(double));
  double *B = A + _3problemSize; //malloc(3*problemSize*sizeof(double));
  double *rho = B + _3problemSize; //malloc(problemSize*sizeof(double));
  double *AWpB = rho + problemSize;// malloc(problemSize*problemSize*sizeof(double));
  int *ipiv = (int *)(AWpB + problemSize2);  // malloc(problemSize*sizeof(int));
#ifndef NDEBUG
  double *AWpB_ = (double *) ipiv + problemSize;
#endif

  for (unsigned int i = 0; i < problemSize; ++i) rho[i] = 1.;

  info[0] = 1;

  // velocity <- M*reaction + qfree


  cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
  cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1.,
        problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);


  while (iter++ < itermax)
  {

    frictionContact3D_FischerBurmeisterFunction(
      problemSize,
      computeACFun3x3,
      reaction, velocity,
      problem->mu, rho,
      F, A, B);

    // AW + B
    computeAWpB(problemSize, A, problem->M->matrix0, B, AWpB);

    int info2 = 0;

    cblas_dcopy(problemSize, F, 1, tmp1, 1);
    cblas_dscal(problemSize, -1., tmp1, 1);

    if (options->iparam[2])
    {
      DGELS(LA_NOTRANS,problemSize, problemSize, 1, AWpB, problemSize,
            tmp1, problemSize, &info2);
    }
    else
    {

#ifndef NDEBUG
      cblas_dcopy(problemSize * problemSize, AWpB, 1, AWpB_, 1);
#endif

      // tmp1 <- sol (AWpB * tmp1 = -F)
      DGESV(problemSize, 1, AWpB, problemSize, ipiv,
            tmp1, problemSize, &info2);

#ifndef NDEBUG
      cblas_dcopy(problemSize, F, 1, tmp2, 1);
      cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1., AWpB_,
            problemSize, tmp1, 1, 1., tmp2, 1);
      assert(cblas_ddot(problemSize, tmp2, 1, tmp2, 1) < 1e-10);
#endif

    }

    assert(info2 >= 0);

    if (info2 > 0)
      /*if (verbose>0)*/
      printf("LOCALFB: warning DGESV failed with U(%d,%d) == 0.\n", info2, info2);

    // line search
    double alpha = 1;
    int info_ls = globalLineSearchGP(problemSize, computeACFun3x3, reaction, velocity, problem->mu, rho, F, A, B,
                                     problem->M->matrix0, problem->q, AWpB, tmp1, tmp2, &alpha, options->iparam[12]);

    if (!info_ls)
      cblas_daxpy(problemSize, alpha, tmp1, 1, reaction, 1);
    else
      cblas_daxpy(problemSize, 1, tmp1, 1., reaction, 1);


    // velocity <- M*reaction + qfree
    cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
    cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1.,
          problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);

    options->dparam[1] = INFINITY;

    if (!(iter % erritermax))
    {
      frictionContact3D_FischerBurmeisterFunction(
        problemSize,
        computeACFun3x3,
        reaction, velocity,
        problem->mu, rho,
        F, NULL, NULL);


      frictionContact3D_FischerBurmeister_compute_error(problem, reaction, velocity,
                                                        tolerance, options, &(options->dparam[1]));


      assert((cblas_dnrm2(problemSize, F, 1)
              / (1 + cblas_dnrm2(problemSize, problem->q, 1)))
             <= (10 * options->dparam[1] + 1e-15));

    }

    if (verbose > 0)
      printf("LOCALFB: iteration %d : error=%g\n", iter, options->dparam[1]);

    if (options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }



  if (verbose > 0)
  {
    if (!info[0])
      printf("LOCALFB: convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    else
    {
      printf("LOCALFB: no convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    }
  }

#ifdef DUMP_PROBLEM
  if (info[0])
  {
    static int file_counter = 0;
    char filename[64];
    printf("LOCALFB: dumping problem\n");
    sprintf(filename, "LOCALFB_failure%d.dat", file_counter++);
    FILE* file = fopen(filename, "w");
    frictionContact_printInFile(problem, file);
    fclose(file);
  }
#endif

  if (!options->dWork)
  {
    assert(buffer);
    free(buffer);
  }
  else
  {
    assert(buffer == options->dWork);
  }


}

int frictionContact3D_FischerBurmeister_compute_error(
    FrictionContactProblem* problem,
    double *z , double *w, double tolerance,
    SolverOptions * options, double * error)
{

  double *A = NULL;
  double *B = NULL;

  unsigned int problemSize = 3 * problem->numberOfContacts;

  double *rho = (double*) malloc(problemSize*sizeof(double));
  double *F = (double *) malloc(problemSize*sizeof(double));

  FischerBurmeisterFun3x3Ptr computeACFun3x3;

  switch (options->iparam[10])
  {
  case 0:
  {

    computeACFun3x3 = &frictionContact3D_FischerBurmeisterFunctionGenerated;
    break;
  }
  }

  frictionContact3D_FischerBurmeisterFunction(
    problemSize,
    computeACFun3x3,
    z, w,
    problem->mu, rho,
    F, A, B);

  *error=0.;
  for(unsigned int i=0; i<problemSize;
      i+=3)
  {
    *error += sqrt(F[i]*F[i] + F[i+1]*F[i+1] + F[i+2]*F[i+2]);
  }

  *error /= (problem->numberOfContacts + 1);

  free(F);
  free(rho);

  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - FrictionContact3D_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
  {
    return 0;
  }
}

int frictionContact3D_FischerBurmeister_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the LOCALFB Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_LOCALFB;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 14;
  options->dSize = 14;
  options->iparam = (int *) malloc(options->iSize * sizeof(int));
  options->dparam = (double *) malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;
  for (unsigned int i = 0; i < 14; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 200;
  options->iparam[1] = 1;
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;
  options->iparam[7] = 1;      /* erritermax */
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1;      /* default rho */

  options->iparam[8] = -1;     /* mpi com fortran */
  options->iparam[10] = 0;
  options->iparam[11] = 0;     /* 0 GoldsteinPrice line search, 1 FBLSA */
  options->iparam[12] = 100;   /* max iter line search */

  options->internalSolvers = NULL;

  return 0;
}


#ifdef WITH_MUMPS

void frictionContact3D_sparseLocalFischerBurmeisterInit(
  SolverOptions *options)
{
  frictionContactNonsmoothEqnInit(options);
}

typedef struct
{
  FischerBurmeisterFun3x3Ptr computeACFun3x3;
} FischerBurmeisterParams;

void nonsmoothEqnFischerBurmeisterFun(void* arg,
                                      unsigned int problemSize,
                                      double* reaction,
                                      double* velocity,
                                      double* mu,
                                      double* rho,
                                      double* result,
                                      double* A,
                                      double* B);
void nonsmoothEqnFischerBurmeisterFun(void* arg,
                                      unsigned int problemSize,
                                      double* reaction,
                                      double* velocity,
                                      double* mu,
                                      double* rho,
                                      double* result,
                                      double* A,
                                      double* B)
{
  FischerBurmeisterParams* acparams_p = (FischerBurmeisterParams *) arg;

  frictionContact3D_FischerBurmeisterFunction(problemSize,
                                         acparams_p->computeACFun3x3,
                                         reaction,
                                         velocity,
                                         mu,
                                         rho,
                                         result,
                                         A,
                                         B);
}




void frictionContact3D_sparseLocalFischerBurmeister(
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

  assert(!problem->M->matrix0);
  assert(problem->M->matrix1);

  assert(!options->iparam[4]); // only host

  FischerBurmeisterParams acparams;

  switch (options->iparam[10])
  {
  case 0:
  {
    acparams.computeACFun3x3 = &frictionContact3D_FischerBurmeisterFunctionGenerated;
    break;
  }
  }

  FrictionContactNonsmoothEqn equation;

  equation.problem = problem;
  equation.data = (void *) &acparams;
  equation.function = &nonsmoothEqnFischerBurmeisterFun;

  frictionContactNonsmoothEqnSolve(&equation, reaction, velocity, info,
                                   options);

}

#else /*WITH_MUMPS*/

void frictionContact3D_sparseLocalFischerBurmeisterInit(
  SolverOptions *options)
{
  fprintf(stderr, "The sparse global Fischer & Burmeister solver needs -DWITH_MUMPS for the compilation of Siconos/Numerics\n");
}

void frictionContact3D_sparseLocalFischerBurmeister(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
{
  fprintf(stderr, "The sparse global Fischer & Burmeister solver needs -DWITH_MUMPS for the compilation of Siconos/Numerics\n");
}
#endif
