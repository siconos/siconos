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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "Relay_Solvers.h"
#include "relay_cst.h"
#include <assert.h>
#include "SiconosBlas.h"
#include "misc.h"

void relay_pgs(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{


  double* M = problem->M->matrix0;
  double* q = problem->q;
  int n = problem -> size;
  double *a = problem->lb;
  double *b = problem->ub;

  assert(M);
  assert(q);
  assert(n);
  assert(a);
  assert(b);



  int itermax = options->iparam[0];
  double tol = options->dparam[0];


  int i;
  double zi;
  double * diag = (double*)malloc(n * sizeof(double));


  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(M[i * n + i]) < DBL_EPSILON)
    {
      if (verbose > 0)
      {
        printf("Numerics::lcp_pgs, error: vanishing diagonal term \n");
        printf(" The problem cannot be solved with this method \n");
      }

      *info = 2;
      free(diag);
      return;
    }
    else diag[i] = 1.0 / M[i * n + i];
  }

  /* Iterations*/
  int iter = 0;
  double err  = 1.;
  *info = 1;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Initialization of w with q */
    cblas_dcopy(n , q , 1 , w , 1);

    for (i = 0 ; i < n ; ++i)
    {
      z[i] = 0.0;
      zi = -(q[i] + cblas_ddot(n , &M[i] , n , z , 1)) * diag[i];
      z[i] = zi;
      if (zi < a[i]) z[i] = a[i];
      else if (zi > b[i]) z[i] = b[i];
    }
    /* **** Criterium convergence **** */
    *info = relay_compute_error(problem, z, w, tol, &err);

    if (verbose == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", z[i]);
      for (i = 0 ; i < n ; ++i) printf(" %g", w[i]);
      printf("\n");
    }
  }

  options->iparam[1] = iter;
  options->dparam[1] = err;


  if (err > tol)
  {
    printf("Siconos/Numerics: relay_pgs: No convergence of PGS after %d iterations\n" , iter);
    printf("The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: relay_pgs: Convergence of PGS after %d iterations\n" , iter);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }




  free(diag);



}
int relay_pgs_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PGS Solver\n");
  }
  /*  strcpy(options->solverName,"PGS");*/
  options->solverId = SICONOS_RELAY_PGS;
  options->numberOfInternalSolvers = 0;
  options->internalSolvers = NULL;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-6;
  options->dparam[1] = 1.0;

  return 0;
}

