/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include <assert.h>                        // for assert
#include <float.h>                         // for DBL_EPSILON
#include <math.h>                          // for fabs, sqrt, INFINITY
#include <stdio.h>                         // for NULL, fprintf, printf, stderr
#include <stdlib.h>                        // for free, malloc, calloc
#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "Friction_cst.h"                  // for SICONOS_FRICTION_3D_IPARAM...
#include "LCP_Solvers.h"                   // for lcp_nsgs_SBM_buildLocalPro...
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix, RawNumeric...
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "SparseBlockMatrix.h"             // for SparseBlockStructuredMatrix
#include "fc2d_Solvers.h"                  // for fc2d_sparse_nsgs, fc2d_spa...
#include "fc2d_compute_error.h"            // for fc2d_compute_error
#include "numerics_verbose.h"              // for numerics_printf, verbose
#include "SiconosBlas.h"                         // for cblas_dnrm2


#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

static
void accumulateLightErrorSum(double *light_error_sum, double localreaction[3],
                             double *oldreaction)
{
  *light_error_sum += ((oldreaction[0] - localreaction[0])*(oldreaction[0] - localreaction[0]) +
                       (oldreaction[1] - localreaction[1])*(oldreaction[1] - localreaction[1])) ;
}
static
double calculateLightError(double light_error_sum, unsigned int nc, double *reaction)
{
  double error = sqrt(light_error_sum);
  double norm_r = cblas_dnrm2(nc*2, reaction, 1);
  if(fabs(norm_r) > DBL_EPSILON)
    error /= norm_r;
  return error;
}
static
int determine_convergence(double error, double tolerance, int iter,
                          SolverOptions *options)
{
  int hasNotConverged = 1;
  if(error < tolerance)
  {
    hasNotConverged = 0;
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e < %7.3e\n", iter, error, tolerance);
  }
  else
  {
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e > %7.3e\n", iter, error, tolerance);
  }
  return hasNotConverged;
}

static
double calculateFullErrorFinal(FrictionContactProblem *problem, SolverOptions *options,
                               double *reaction, double *velocity, double tolerance,
                               double norm_q)
{
  double absolute_error;
  /* (*computeError)(problem, reaction , velocity, tolerance, */
  /*                 options, norm_q, &absolute_error); */

  fc2d_compute_error(problem, reaction, velocity, tolerance, norm_q, &absolute_error);


  if(verbose > 0)
  {
    if(absolute_error > options->dparam[SICONOS_DPARAM_TOL])
    {
      numerics_printf("-- FC2D - NSGS - Warning absolute "
                      "Residual = %14.7e is larger than required precision = %14.7e",
                      absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
    else
    {
      numerics_printf("-- FC2D - NSGS - absolute "
                      "Residual = %14.7e is smaller than required precision = %14.7e",
                      absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
  }
  return absolute_error;
}



static
int determine_convergence_with_full_final(FrictionContactProblem *problem, SolverOptions *options,
    double *reaction, double *velocity,
    double *tolerance, double norm_q, double error,
    int iter)
{
  int hasNotConverged = 1;
  if(error < *tolerance)
  {
    hasNotConverged = 0;
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e < %7.3e", iter, error, *tolerance);

    double absolute_error = calculateFullErrorFinal(problem, options,
                            reaction, velocity,
                            options->dparam[SICONOS_DPARAM_TOL],
                            norm_q);
    if(absolute_error > options->dparam[SICONOS_DPARAM_TOL])
    {
      *tolerance = error/absolute_error*options->dparam[SICONOS_DPARAM_TOL];
      assert(*tolerance > 0.0 && "tolerance has to be positive");
      /* if (*tolerance < DBL_EPSILON) */
      /* { */
      /*   numerics_warning("determine_convergence_with_full_fina", "We try to set a very smal tolerance"); */
      /*   *tolerance = DBL_EPSILON; */
      /* } */
      numerics_printf("-- FC2D - NSGS - We modify the required incremental precision to reach accuracy to %e", *tolerance);
      hasNotConverged = 1;
    }
    else
    {
      numerics_printf("-- FC2D - NSGS - The incremental precision is sufficient to reach accuracy to %e", *tolerance);
    }




  }
  else
  {
    numerics_printf("-- FC2D - NSGS - Iteration %i "
                    "Residual = %14.7e > %7.3e", iter, error, *tolerance);
  }
  return hasNotConverged;
}


static int fc2dLocalSolve(double *W, double *q, double mu, double *P, double *U);

int fc2dLocalSolve(double *W, double *q, double mu, double *P, double *U)
{
  double D, muPn;

  /* | Wnn Wnt |
     | Wtn Wtt | */

#define Wnn W[0]
#define Wtn W[1]
#define Wnt W[2]
#define Wtt W[3]

  if(q[0] > 0)
  {
    P[0] = 0;
    P[1] = 0;
  }
  else
  {
    /* solve WP + q = 0  */
    D = Wnn * Wtt - Wnt * Wtn;
    if(D < DBL_EPSILON) return(1);

    P[0] = - (Wtt * q[0] - Wnt * q[1]) / D;
    P[1] = - (-Wtn * q[0] + Wnn * q[1]) / D;

    muPn = mu * P[0];

    if(fabs(P[1]) > muPn)
      /* outside cone */
    {

      if(P[1] + muPn < 0)
      {

        P[0] = - q[0] / (Wnn - mu * Wnt);
        P[1] = - mu * P[0];
      }
      else
      {

        P[0] = - q[0] / (Wnn + mu * Wnt);
        P[1] = mu * P[0];

      }
    }
  }


#undef Wnn
#undef Wnt
#undef Wtn
#undef Wtt

  return(0);
}


void fc2d_sparse_nsgs(FrictionContactProblem* problem, double *z, double *w,
                      int *info, SolverOptions* options)
{
  /* Notes:
     - we suppose that the trivial solution case has been checked before,
     and that all inputs differs from NULL since this function is
     supposed to be called from lcp_driver_global().

     - Input matrix M of the problem is supposed to be sparse-block
     with no null row (ie no rows with all blocks equal to null)
  */

  assert(problem->M->matrix1);

  /*
    The options for the global "block" solver are defined in
    options[0].
   */

  /* Global Solver parameters*/
  int * iparam = options->iparam;
  double * dparam = options->dparam;


  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  /* Matrix M/vector q of the LCP */
  SparseBlockStructuredMatrix* blmat = problem->M->matrix1;
  double * q = problem->q;

  int nc = problem->numberOfContacts;
  double norm_q = cblas_dnrm2(nc*2, problem->q, 1);


  assert(blmat->nbblocks >= 1);

  /* Local problem initialization */

  LinearComplementarityProblem * local_problem = (LinearComplementarityProblem *)
      malloc(sizeof(*local_problem));
  local_problem->M = (NumericsMatrix *)malloc(sizeof(*local_problem->M));
  local_problem->M->storageType = 0; // dense storage
  local_problem->M->matrix0 = NULL;
  local_problem->M->matrix1 = NULL;
  local_problem->M->matrix2 = NULL;
  local_problem->M->internalData = NULL;

  /* Memory allocation for q. Size of q = blsizemax, size of the
     largest square-block in blmat */

  int blsizemax = blmat->blocksize0[0];
  int k;
  for(unsigned int i = 1 ; i < blmat->blocknumber0 ; i++)
  {
    k = blmat->blocksize0[i] - blmat->blocksize0[i - 1];
    if(k > blsizemax) blsizemax = k;
  }
  local_problem->q = (double*)malloc(blsizemax * sizeof(double));
  double localreaction[2];

  /* Current row (of blocks) number */
  unsigned int rowNumber;

  /*****  Gauss-Seidel iterations *****/
  int iter = 0; /* Current iteration number */
  double error = INFINITY; /* Current error */
  int hasNotConverged = 1;

  int pos = 0;

  /* Output from local solver */
  int infoLocal = -1;

  while((iter < itermax) && hasNotConverged)
  {
    ++iter;

    double light_error_sum = 0.0;
    /* Loop over the rows of blocks in blmat */
    for(pos = 0, rowNumber = 0; rowNumber < blmat->blocknumber0; ++rowNumber, ++pos, ++pos)
    {
      /* Local problem formalization */
      lcp_nsgs_SBM_buildLocalProblem(rowNumber, blmat, local_problem, q, z);


      localreaction[0] = z[pos];
      localreaction[1] = z[pos+1];


      /* Solve local problem */
      infoLocal = fc2dLocalSolve(local_problem->M->matrix0,
                                 local_problem->q,
                                 problem->mu[rowNumber],
                                 localreaction, &w[pos]);


      if(infoLocal)
      {
        free(local_problem->q);
        free(local_problem->M);
        free(local_problem);
        /* Number of GS iterations */
        options[0].iparam[SICONOS_IPARAM_ITER_DONE] = iter;
        fprintf(stderr, "fc2d_nsgs error: local LCP solver failed at global iteration %d.\n for block-row number %d. Output info equal to %d.\n", iter, rowNumber, infoLocal);

        *info = infoLocal;
        return;
      }
      if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT ||
          iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL
        )
        accumulateLightErrorSum(&light_error_sum, localreaction, &z[pos]);

      z[pos]   = localreaction[0];
      z[pos+1] = localreaction[1];

    }
    /* error evaluation */
    if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
    {
      error = calculateLightError(light_error_sum, nc, z);
      hasNotConverged = determine_convergence(error, tolerance, iter, options);
    }
    else if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
    {
      error = calculateLightError(light_error_sum, nc, z);
      hasNotConverged = determine_convergence_with_full_final(problem,  options,
                        z, w,
                        &tolerance, norm_q, error,
                        iter);

    }

    /* if (erriter >= erritermax) */
    /* { */
    /*   erriter = 0; */
    /*   fc2d_compute_error(problem, z, w, tolerance, &error); */
    /*   hasNotConverged = error > tolerance  ; */
    /* } */
  }
  /* Full criterium */
  if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
  {
    error = calculateFullErrorFinal(problem, options, z, w,
                                    tolerance, norm_q);

    hasNotConverged = determine_convergence(error,  dparam[SICONOS_DPARAM_TOL], iter, options);


  }


  // numerics_printf("Siconos Numerics : problem size=%d, nb iterations=%d, error=%g\n",
  //          blmat->blocknumber0,
  //          iter,
  //          error);

  *info = hasNotConverged;

  /* Number of GS iterations */
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /* Resulting error */
  dparam[SICONOS_DPARAM_RESIDU] = error;

  free(local_problem->q);
  free(local_problem->M);
  free(local_problem);
}

// options setup is done through fc2d_nsgs_set_default.

