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
#include "SiconosBlas.h"
#include <math.h>
#include <float.h>
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"

#include "numerics_verbose.h"
#define EPSDIAG DBL_EPSILON
void lcp_rpgs(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;

  double * q = problem->q;

  /* size of the LCP */
  int n = problem->size;

  int incx, incy;
  int i, iter;
  double qs, err, zi;
  double *diag;
  double Mii, ziprev;

  /*  double *buffer_errors;
    FILE *ficbuffer_errors;*/
  /* Recup input */

  int itermax = options->iparam[0];

  /*  buffer_errors = malloc( itermax*sizeof( double ) );*/

  double tol = options->dparam[0];
  double rho = options->dparam[2];
  // double omega = options->dparam[3]; // Not yet used


  /* Initialize output */

  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  /* Allocation */

  /*  ww   = ( double* )malloc( n*sizeof( double ) );*/
  /*  zprev = ( double* )malloc( n*sizeof( double ) );*/
  diag = (double*)malloc(n * sizeof(double));
  /*  diagprev = ( int* )malloc( n*sizeof( int ) );*/

  /*  qs = 0.;*/
  incx = 1;
  qs = cblas_dnrm2(n , q , incx);
  if (verbose > 0) printf("\n ||q||= %g \n", qs);
  // den = 1.0 / qs; // Note FP, den is never used.

  /* Initialization of z & w */
  /*
    for( i = 0 ; i < n ; i++ ) z[i] = 0;

    incx = 1;
    incy = 1;
    dcopy_( (integer *)&n , q , (integer *)&incx , w , (integer *)&incy );
  */
  /* Preparation of the diagonal of the inverse matrix */
  /*  maxdiag = 0.0;
    for( i = 0 ; i < n ; i++ ){
      if (M[i*(n+1)] > maxdiag) maxdiag = M[i*(n+1)];
    }

    if (maxdiag > EPSDIAG) rho = maxdiag;
    invrho = 1.0/rho;*/

  for (i = 0 ; i < n ; ++i)
  {
    Mii = M[i * (n + 1)];
    /* if(abs(Mii+rho)<EPSDIAG){   */ /* Version of Pascal */
    if (Mii < -EPSDIAG)
    {
      if (verbose > 0)
      {
        printf(" RPGS : Warning negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diag);
      /*      free(zprev);*/
      /*      free(ww);*/
      /*      free(diagprev);*/
      /*      free(buffer_errors);*/
      return;
    }
    else
    {
      diag[i] = 1.0 / (Mii + rho);
      /*        qs += pow(q[i]*diag[i],2);*/
      /*        if (Mii < EPSDIAG ){
                  diag[i] = invrho;
                  diagprev[i] = 1;
              }
              else {
                  diag[i] = 1.0/Mii;
                  diagprev[i] = 0;
              }
      */
    }
  }
  /*  den = 1.0/sqrt(qs);*/

  /*start iterations*/

  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = n;
    incy = 1;
    for (i = 0 ; i < n ; ++i)
    {
      ziprev = z[i];
      z[i] = 0.0;

      zi = -(q[i] - (rho * ziprev) + cblas_ddot(n , &M[i] , incx , z , incy)) * diag[i];

      if (zi > 0) z[i] = zi;

    }
    /* **** Criterium convergence **** */
    lcp_compute_error(problem, z, w, tol, &err);

    /*    buffer_errors[iter-1] = err;*/

    if (verbose == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", z[i]);
      for (i = 0 ; i < n ; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* **** ********************* **** */

  }

  options->iparam[1] = iter;
  options->dparam[1] = err;

  if (verbose > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of RPGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of RPGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }

  /*  if (iter == itermax)
    {
          ficbuffer_errors = fopen("errors.dat","w");
          if (ficbuffer_errors == NULL)
          {
              printf("Impossible d'ouvrir errors.dat !\n");
          } else {
              for(i = 0 ; i < itermax ; i++)
              {
                  fprintf(ficbuffer_errors,"%g\n",buffer_errors[i]);
              }
              fclose(ficbuffer_errors);
          }
    }
  */
  /*  free(ww);*/
  free(diag);
  /*  free(zprev);*/
  /*   free(diagprev);*/
  /*  free(buffer_errors);*/

  return;
}
int linearComplementarity_rpgs_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the RPGS Solver\n");
  }


  /*  strcpy(options->solverName,"RPGS");*/
  options->solverId = SICONOS_LCP_RPGS;
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
  options->dparam[1] = 1.0;


  return 0;
}
