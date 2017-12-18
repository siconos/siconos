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


#include "SecondOrderConeLinearComplementarityProblem.h"
#include "SolverOptions.h"
#include "soclcp_compute_error.h"
#include "soclcp_projection.h"
#include "projectionOnCone.h"
#include "projectionOnCylinder.h"
#include "SiconosLapack.h"
#include "numerics_verbose.h"
#include <math.h>
#include <assert.h>

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
void soclcp_unitary_compute_and_add_error(double *z , double *w, unsigned int dim, double mu, double * error,
                                          double * worktmp)
{

  double rho = 1.0;
  for (unsigned int i =0; i < dim; ++i)
  {
    worktmp[i] = z[i] - rho * w[i];
  }
  /* printf("mu = %f\n", mu); */
  /* for (int i=0; i < dim; i++ ) */
  /* { */
  /*   printf("-- worktmp[%i]=%e\t\t\t",i,(worktmp)[i]); */
  /*   printf("z[%i]=%e\n",i,(z)[i]); */
  /* } */
  projectionOnSecondOrderCone(worktmp, mu, dim);
  /* for (int i=0; i < dim; i++ ) */
  /* { */
  /*   printf("-- worktmp[%i]=%e\t\t\n",i,(worktmp)[i]); */
  /* } */

  for (unsigned int i =0; i < dim; ++i)
  {
    worktmp[i] = z[i] -  worktmp[i];
    *error +=  worktmp[i] * worktmp[i];
  }
}
int soclcp_compute_error(
  SecondOrderConeLinearComplementarityProblem* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options, double * error)
{
  assert(problem);
  assert(z);
  assert(w);
  assert(error);

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->nc;
  double *mu = problem->tau;
  int n = problem->n;
  
  cblas_dcopy(n , problem->q , incx , w , incy); // w <-q
  // Compute the current velocity
  NM_gemv(1.0, problem->M, z, 1.0, w);

  /* for (int i=0; i < n ; i++ ) */
  /* { */
  /*   printf("w[%i]=%e\t\t\t",i,w[i]); */
  /*   printf("z[%i]=%e\n",i,z[i]); */
  /* } */
  /* printf("\n"); */
  
  *error = 0.;

  int ic;
  int dim;
  unsigned int dim_max=0;
  for (int i =0; i <nc; i++)
  {
    dim_max=max(dim_max,problem->coneIndex[i+1]-problem->coneIndex[i]);
  }
  double *worktmp = (double *)calloc(dim_max,sizeof(double));

  for(ic = 0 ; ic < nc ; ic++)
  {
    dim = problem->coneIndex[ic+1]- problem->coneIndex[ic];
    soclcp_unitary_compute_and_add_error(z + problem->coneIndex[ic],
                                         w + problem->coneIndex[ic],
                                         dim , mu[ic], error, worktmp);
    /* for (int i=0; i < dim; i++ ) */
    /* { */
    /*   printf("-- w[%i]=%e\t\t\t",i,(w + problem->coneIndex[ic])[i]); */
    /*   printf("z[%i]=%e\n",i,(z + problem->coneIndex[ic])[i]); */
    /* } */
  }
  free(worktmp);
  *error = sqrt(*error);

  /* Computes error */
  double norm_q = cblas_dnrm2(n , problem->q , incx);
  DEBUG_PRINTF("norm_q = %12.8e\n", norm_q);
  *error = *error / (norm_q + 1.0);
  
  if(*error > tolerance)
  {
    if(verbose > 1)
      printf(" Numerics - soclcp_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
    return 0;
}



int soclcp_compute_error_v(SecondOrderConeLinearComplementarityProblem* problem, double *z , double *w, double tolerance, SolverOptions *options, double * error)
{
  /* Checks inputs */
  if(problem == NULL || z == NULL || w == NULL)
    numerics_error("soclcp_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->nc;
  int n = problem->n;
  double *mu = problem->tau;

  double invmu = 0.0;
  cblas_dcopy(n , problem->q , incx , z , incy); // z <-q

  // Compute the current reaction
  NM_gemv(1.0, problem->M, w, 1.0, z);

  *error = 0.;
  double rho = 1.0;
  for(int ic = 0 ; ic < nc ; ic++)
  {
    int dim = problem->coneIndex[ic+1]-problem->coneIndex[ic];
    double * worktmp = (double *)malloc(dim*sizeof(double)) ;
    int nic = problem->coneIndex[ic];
    for (int i=0; i < dim; i++)
    {
      worktmp[i] = w[nic+i] - rho * z[nic+i];
    }
    invmu = 1.0 / mu[ic];
    projectionOnSecondOrderCone(worktmp, invmu, dim);
    for (int i=0; i < dim; i++)
    {
      worktmp[i] = w[nic+i] - worktmp[i];
      *error +=  worktmp[i] * worktmp[i];
    }
    free(worktmp);
  }
  *error = sqrt(*error);

  /* Computes error */
  double norm_q = cblas_dnrm2(n , problem->q , incx);
  *error = *error / (norm_q + 1.0);
  if(*error > tolerance)
  {
    /*      if (verbose > 0) printf(" Numerics - soclcp_compute_error_velocity failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}
