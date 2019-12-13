
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <assert.h>                                    // for assert
#include <math.h>                                      // for sqrt
#include <stdio.h>                                     // for fprintf, stderr
#include <stdlib.h>                                    // for malloc, free
#include "ConvexQP.h"                                  // for ConvexQP
#include "ConvexQP_Solvers.h"                          // for convexQP_ADMM_...
#include "ConvexQP_cst.h"                              // for SICONOS_CONVEX...
#include "Friction_cst.h"                              // for SICONOS_FRICTI...
#include "GlobalFrictionContactProblem.h"              // for GlobalFriction...
#include "GlobalFrictionContactProblem_as_ConvexQP.h"  // for GlobalFriction...
#include "NumericsFwd.h"                               // for ConvexQP, Solv...
#include "NumericsMatrix.h"                            // for NM_clear, NM_tr...
#include "SiconosBlas.h"                               // for cblas_dcopy
#include "SolverOptions.h"                             // for SolverOptions
#include "debug.h"                                     // for DEBUG_EXPR
#include "gfc3d_Solvers.h"                             // for gfc3d_set_inte...
#include "gfc3d_compute_error.h"                       // for gfc3d_compute_...
#include "numerics_verbose.h"                          // for numerics_print...


/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalSolverPtr)(ConvexQP* ,
                                  double *, double *, double *, double *,
                                  int* , SolverOptions* );
void gfc3d_ACLMFixedPoint(GlobalFrictionContactProblem* restrict problem, double* restrict reaction, double* restrict velocity,
                          double* restrict globalVelocity, int* restrict info, SolverOptions* restrict options)
{

  /* verbose=1; */
  
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(n , problem->q , 1);
  double norm_b = cblas_dnrm2(m , problem->b , 1);



  if (options->numberOfInternalSolvers < 1)
  {
    numerics_error("gfc3d_ACLMFixedpoint", "The ACLM Fixed Point method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }

  SolverOptions * internalsolver_options = options->internalSolvers[0];

  if (verbose > 0)
  {
    solver_options_print(options);
  }


  /*****  Fixed Point Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  internalSolverPtr internalsolver;

  ConvexQP *cqp = (ConvexQP *)malloc(sizeof(ConvexQP));
  cqp->size = n;
  cqp->m = m;
  cqp->M = problem->M;

  cqp->q = (double *) malloc(n* sizeof(double));
  for(int i = 0; i < n ; i++)
  {
    cqp->q[i] = -problem->q[i];
  }
  DEBUG_EXPR(NM_display(problem->H));

  cqp->A = NM_transpose(problem->H);
  DEBUG_EXPR(NM_display(cqp->A));
  cqp->b = (double *) malloc(m* sizeof(double));
  cqp->ProjectionOnC = &Projection_ConvexQP_GFC3D_DualCone;
  cqp->normConvexQP= norm_q;
  cqp->istheNormConvexQPset=1;
  double * w = (double *) malloc(n* sizeof(double));

  GlobalFrictionContactProblem_as_ConvexQP *gfc3d_as_cqp= (GlobalFrictionContactProblem_as_ConvexQP*)malloc(sizeof(GlobalFrictionContactProblem_as_ConvexQP));
  cqp->env = gfc3d_as_cqp ;
  gfc3d_as_cqp->cqp = cqp;
  gfc3d_as_cqp->gfc3d = problem;
  gfc3d_as_cqp->options = options;
  if (internalsolver_options->solverId == SICONOS_CONVEXQP_ADMM  )
  {
    numerics_printf_verbose(1," ========================== set ADMM solver internal ConveQP problem ==========================\n");
    internalsolver = &convexQP_ADMM;
    convexQP_ADMM_init(cqp, options->internalSolvers[0]);
  }
  else
  {
    fprintf(stderr, "Numerics, gfc3d_ACLMFixedPoint failed. Unknown internal solver.\n");
    exit(EXIT_FAILURE);
  }

  double normUT;
  int cumul_iter =0;
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    // internal solver for the regularized problem

    /* Compute the value of the initial value of b */
    cblas_dcopy(m,problem->b,1,cqp->b,1);
    for (int ic = 0 ; ic < nc ; ic++)
    {
      normUT = sqrt(velocity[ic*3+1] * velocity[ic*3+1] + velocity[ic*3+2] * velocity[ic*3+2]);
      cqp->b[3*ic] += problem->mu[ic]*normUT;
    }

    gfc3d_set_internalsolver_tolerance(problem,options,internalsolver_options, error);

    (*internalsolver)(cqp, globalVelocity, w,  reaction , velocity , info , internalsolver_options);

    cumul_iter +=  internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE];
    /* **** Criterium convergence **** */

    gfc3d_compute_error(problem, reaction , velocity, globalVelocity, tolerance, options,
                        norm_q, norm_b, &error);

    numerics_printf_verbose(1,"---- GFC3D - ACLMFP - Iteration %i Residual = %14.7e", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }

  numerics_printf_verbose(1,"---- GFC3D - ACLMFP - # Iteration %i Final Residual = %14.7e", iter, error);
  numerics_printf_verbose(1,"---- GFC3D - ACLMFP - #              internal iteration = %i", cumul_iter);

  NM_clear(cqp->A);
  free(cqp->b);
  free(cqp->q);
  free(w);
  free(cqp);


 if (internalsolver_options->solverId == SICONOS_CONVEXQP_ADMM  )
  {
    convexQP_ADMM_free(cqp, options->internalSolvers[0]);
  }

  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

}



void gfc3d_aclmfp_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] =2.0;

  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_CONVEXQP_ADMM);
  options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] =1000;
}
