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


#include "VariationalInequality.h"
#include "SolverOptions.h"
#include "VariationalInequality_computeError.h"
#include "SiconosLapack.h"
#include "SiconosSets.h"

#include <math.h>
#include <assert.h>
#include <float.h>

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
#include "numerics_verbose.h"

int variationalInequality_computeError(
  VariationalInequality* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options, double * error)
{

  assert(problem);
  assert(z);
  assert(w);
  assert(error);
  
  int incx = 1;
  int n = problem->size;

  *error = 0.;
  if (!options->dWork)
  {
    options->dWork = (double*)calloc(2*n,sizeof(double));
  }
  double *ztmp =  options->dWork;
  double *wtmp =  &(options->dWork[n]);

  
  if (!problem->istheNormVIset)
  {
    for (int i=0;i<n;i++)
    {
      ztmp[i]=0.0 ;
    }
    problem->F(problem,n,ztmp,w);
    problem->normVI= cblas_dnrm2(n , w , 1);
    DEBUG_PRINTF("problem->normVI = %12.8e\n", problem->normVI);
    problem->istheNormVIset=1;
  }

  double norm_q =problem->normVI;
  DEBUG_PRINTF("norm_q = %12.8e\n", norm_q);

  problem->F(problem,n,z,w);
  
  cblas_dcopy(n , z , 1 , ztmp, 1);
  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;

  problem->ProjectionOnX(problem,ztmp,wtmp);

  cblas_daxpy(n, -1.0, z , 1, wtmp , 1) ;
  *error = cblas_dnrm2(n , wtmp , incx);

  /* Computes error */
  if (fabs(norm_q) > DBL_EPSILON)
    *error /= norm_q;

  DEBUG_PRINTF("error = %e\n",*error);
  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - variationalInequality_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
    return 0;
}

int variationalInequality_compute_error_box(
  VariationalInequality* problem,
  double* x, double* F, double tolerance, double * error)
{
  assert(problem);
  assert(x);
  assert(F);
  assert(error);
  assert(problem->set);

  double* lb = ((box_constraints*) problem->set)->lb;
  double* ub = ((box_constraints*) problem->set)->ub;
  double diff;
  double err = 0;

  // compute componentwise \Pi_box(x-F(x)) - x
  for (int i = 0; i < problem->size; ++i)
  {
    diff = x[i] - F[i];
    if (diff < lb[i])
    {
      diff = lb[i] - x[i];
    }
    else if (diff > ub[i])
    {
      diff = ub[i] - x[i];
    }
    else
    {
      diff = F[i]; // should be -F, but we square it anyway ...
    }
    err += diff*diff;
  }
  error[0] = sqrt(err);

  if (error[0] > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - variationalInequality_compute_error: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}

/*
int variationalInequality_computeError_wait(
  VariationalInequality* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options, double * error)
{

  assert(problem);
  assert(z);
  assert(w);
  assert(error);
  
  int incx = 1;

  int n = problem->size;

  *error = 0.;
  double *ztmp = (double*)malloc(n* sizeof(double));
  cblas_dcopy(n , z , 1 , ztmp, 1);

  problem->F(problem,ztmp,w);
  double norm_q = cblas_dnrm2(n , w , incx);
  
  
  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;
  
  problem->ProjectionOnX(problem,ztmp,w);
  
  cblas_dcopy(n , z , 1 , ztmp, 1);

  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;

  *error = cblas_dnrm2(n , ztmp , incx);
  free(ztmp);
  
  problem->F(problem,z,w);


  
  // Computes error
  *error = *error / (norm_q + 1.0);
  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - variotionalInequality_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
    return 0;
}
*/
