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

#include "debug.h"
#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_compute_error.h"
#include "AlartCurnierGenerated.h"
#include "FrictionContact3D_nonsmooth_Newton_solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "SiconosLapack.h"


void frictionContact3D_AlartCurnierFunction(
  unsigned int problemSize,
  AlartCurnierFun3x3Ptr computeACFun3x3,
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

int frictionContact3D_AlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the LOCALAC Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_LOCALAC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 14;
  options->dSize = 14;
  options->iparam = (int *) malloc(options->iSize * sizeof(int));
  options->dparam = (double *) malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  null_SolverOptions(options);
  for (unsigned int i = 0; i < 14; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 200;
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;
  options->iparam[7] = 1;      /* erritermax */
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1;      /* default rho */

  options->iparam[8] = -1;     /* mpi com fortran */
  options->iparam[10] = 2;     /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */
  options->iparam[11] = 0;     /* 0 GoldsteinPrice line search, 1 FBLSA */
  options->iparam[12] = 100;   /* max iter line search */

#ifdef WITH_MUMPS
  options->iparam[13] = 1;
#else
  options->iparam[13] = 0;     /* Linear solver used at each Newton iteration. 0: cs_lusol, 1 mumps */
#endif

  options->internalSolvers = NULL;

#ifdef HAVE_MPI
  options->solverData = MPI_COMM_NULL;
#endif

  return 0;
}


typedef struct
{
  AlartCurnierFun3x3Ptr computeACFun3x3;
} AlartCurnierParams;

void nonsmoothEqnAlartCurnierFun(void* arg,
                                   unsigned int problemSize,
                                   double* reaction,
                                   double* velocity,
                                   double* mu,
                                   double* rho,
                                   double* result,
                                   double* A,
                                   double* B);
void nonsmoothEqnAlartCurnierFun(void* arg,
                                 unsigned int problemSize,
                                 double* reaction,
                                 double* velocity,
                                 double* mu,
                                 double* rho,
                                 double* result,
                                 double* A,
                                 double* B)
{
  AlartCurnierParams* acparams_p = (AlartCurnierParams *) arg;

  frictionContact3D_AlartCurnierFunction(problemSize,
                                         acparams_p->computeACFun3x3,
                                         reaction,
                                         velocity,
                                         mu,
                                         rho,
                                         result,
                                         A,
                                         B);
}




void frictionContact3D_nonsmooth_Newton_AlartCurnier(
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

  assert(!options->iparam[4]); // only host

  AlartCurnierParams acparams;

  switch (options->iparam[10])
  {
  case 0:
  {
    acparams.computeACFun3x3 = &computeAlartCurnierSTD;
    break;
  }
  case 1:
  {
    acparams.computeACFun3x3 = &computeAlartCurnierJeanMoreau;
    break;
  };
  case 2:
  {
    acparams.computeACFun3x3 = &frictionContact3D_AlartCurnierFunctionGenerated;
    break;
  }
  case 3:
  {
    acparams.computeACFun3x3 = &frictionContact3D_AlartCurnierJeanMoreauFunctionGenerated;
    break;
  }
  }

  FrictionContact_nonsmooth_Newton_solvers equation;

  equation.problem = problem;
  equation.data = (void *) &acparams;
  equation.function = &nonsmoothEqnAlartCurnierFun;

  frictionContact_nonsmooth_Newton_solvers_solve(&equation, reaction, velocity, info,
                                   options);

}
