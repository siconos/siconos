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
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"

#include "SiconosBlas.h"
#include "numerics_verbose.h"

/*\warning omega is not explicitely used. must be completed    */
void lcp_psor(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;

  double * q = problem->q;

  /* size of the LCP */
  int n = problem->size;
  int incx = 1, incy = 1;
  int incxn;
  int i, iter;

  double qs, err;
  double *ww, *diag;
  int itermax = options->iparam[0];
  double tol = options->dparam[0];
  double omega = options->dparam[2]; // Not yet used
  printf("Warning : omega %f is not used !!!!!\n", omega);

  incxn = n;

  /* Initialize output */

  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  /* Allocation */

  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  /* Check for non trivial case */

  qs = cblas_dnrm2(n , q , incx);

  if (verbose > 0) printf("\n ||q||= %g \n", qs);

  // Note FP : den never used ...
  //if (qs > DBL_EPSILON) den = 1.0 / qs;
  //else
  if (qs <= DBL_EPSILON)
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }

    free(ww);
    free(diag);

    *info = 0;
    return;
  }

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
  }

  cblas_dcopy(n , q , incx , w , incy);
  /* Intialization of w and z */
  /*if(initmethod == 0) {*/
  /* dcopy_( (integer *)&n , q , (integer *)&incx , w , (integer *)&incy );*/
  /*dcopy_( (integer *)&n , ww , (integer *)&incx , z , (integer *)&incy );*/
  /*}    else if(initmethod == 0) { */


  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(M[i * n + i]) < DBL_EPSILON)
    {

      if (verbose > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diag);
      free(ww);

      return;
    }
    else diag[i] = 1.0 / M[i * n + i];
  }

  /*start iterations*/

  iter = 0;
  err  = 1.;

  qs   = -1.0;

  cblas_dcopy(n , q , incx , w , incy);

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    cblas_dcopy(n , w , incx , ww , incy);   /* w --> ww */
    cblas_dcopy(n , q , incx , w , incy);    /* q --> w */

    for (i = 0 ; i < n ; ++i)
    {

      z[i] = 0.0;

      /*      zi = -( q[i] + ddot_( (integer *)&n , &M[i] , (integer *)&incxn , z , (integer *)&incy ))*diag[i]; */

      /*       if( zi < 0 ) z[i] = 0.0;  */
      /*       else z[i] = zi; */

      z[i] = fmax(0.0, -(q[i] + cblas_ddot(n , &M[i] , incxn , z , incy)) * diag[i]);

    }

    /* **** Criterium convergence **** */
    lcp_compute_error(problem, z, w, tol, &err);

    /* **** ********************* **** */
  }

  options->iparam[1] = iter;
  options->dparam[1] = err;

  if (err > tol)
  {
    printf("Siconos/Numerics: lcp_psor: No convergence of PSOR after %d iterations\n" , iter);
    printf("Siconos/Numerics: lcp_psor: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: lcp_psor: Convergence of PSOR after %d iterations\n" , iter);
      printf("Siconos/Numerics: lcp_psor: The residue is : %g \n", err);
    }
    *info = 0;
  }


  free(ww);
  free(diag);

  return;
}
int linearComplementarity_psor_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PSOR Solver\n");
  }


  /*  strcpy(options->solverName,"PSOR");*/
  options->solverId = SICONOS_LCP_PSOR;
  options->numberOfInternalSolvers = 0;
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
  options->dparam[2] = 0.1;


  return 0;
}
