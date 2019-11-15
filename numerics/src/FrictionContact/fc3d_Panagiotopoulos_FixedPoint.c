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

#include <math.h>                                // for fmax
#include <stdio.h>                               // for printf, NULL
#include <stdlib.h>                              // for malloc, calloc, free
#include "ConvexQP.h"                            // for ConvexQP
#include "ConvexQP_Solvers.h"                    // for convexQP_VI_solver_s...
#include "ConvexQP_cst.h"                        // for SICONOS_CONVEXQP_PG
#include "FrictionContactProblem.h"              // for SplittedFrictionCont...
#include "FrictionContactProblem_as_ConvexQP.h"  // for FrictionContactProbl...
#include "Friction_cst.h"                        // for SICONOS_FRICTION_3D_...
#include "LCP_Solvers.h"                         // for linearComplementarit...
#include "LinearComplementarityProblem.h"        // for LinearComplementarit...
#include "NumericsFwd.h"                         // for SolverOptions, ConvexQP
#include "NumericsMatrix.h"                      // for NM_gemv
#include "SolverOptions.h"                       // for SolverOptions, SICON...
#include "fc3d_Solvers.h"                        // for fc3d_set_internalsol...
#include "fc3d_compute_error.h"                  // for fc3d_compute_error
#include "lcp_cst.h"                             // for SICONOS_LCP_CONVEXQP_PG
#include "numerics_verbose.h"                    // for numerics_error, verbose
#include "SiconosBlas.h"                               // for cblas_dcopy, cblas_d...

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*normalInternalSolverPtr)(LinearComplementarityProblem*, double*, double*, int *, SolverOptions *);
typedef void (*tangentInternalSolverPtr)(ConvexQP*, double*, double*, int *, SolverOptions *);


void fc3d_Panagiotopoulos_FixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{

  /* verbose=1; */
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);


  SolverOptions ** internalsolver_options = options->internalSolvers;

  if (verbose) solver_options_print(options);

  /*****  Fixed Point Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  normalInternalSolverPtr internalsolver_normal;
  tangentInternalSolverPtr internalsolver_tangent;
  options->dWork = (double *) calloc(nc, sizeof(double));
  options->dWorkSize = nc;
  double * mu = options->dWork;
  // Warning : same dwork for current and internal solver !!
  internalsolver_options[0]->dWork = options->dWork;

  double * r_n = (double *) malloc(nc * sizeof(double));

  double * r_t = (double *) malloc(2* nc * sizeof(double));
  for (int contact = 0 ; contact < nc; contact ++)
  {
    r_n[contact] = reaction[contact*3];
    r_t[2*contact] = reaction[contact*3+1];
    r_t[2*contact+1] = reaction[contact*3+2];
  }

  SplittedFrictionContactProblem * splitted_problem = (SplittedFrictionContactProblem *)malloc(sizeof(SplittedFrictionContactProblem));

  createSplittedFrictionContactProblem(problem, splitted_problem);

  LinearComplementarityProblem* normal_lcp_problem;
  ConvexQP * tangent_cqp;

  if (options->numberOfInternalSolvers !=2)
    numerics_error("fc3d_Panagiotopoulos_FixedPoint", " the solver requires 2 internal solver");

  if (internalsolver_options[0]->solverId == SICONOS_LCP_PGS||
      internalsolver_options[0]->solverId == SICONOS_LCP_CONVEXQP_PG)
  {
 
    normal_lcp_problem = (LinearComplementarityProblem*)malloc(sizeof(LinearComplementarityProblem));
    normal_lcp_problem->size = nc;
    normal_lcp_problem->M = splitted_problem->M_nn;

    /* for (int contact = 0 ; contact < nc; contact ++) */
    /* { */
    /*   problem->mu[contact] =0.0; */
    /* } */
    /* splitted_problem->M_nn->matrix0 =(double *)malloc(splitted_problem->M_nn->size0*splitted_problem->M_nn->size1*sizeof(double)); */
    /* SBM_to_dense(splitted_problem->M_nn->matrix1, splitted_problem->M_nn->matrix0); */
    /* splitted_problem->M_nn->storageType=NM_DENSE; */

    normal_lcp_problem->q = (double *)malloc(nc*sizeof(double));
    
  }
  else
  {
    numerics_error("fc3d_Panagiotopoulos_FixedPoint", "Unknown internal solver for the normal part.");
  }
 if (internalsolver_options[1]->solverId == SICONOS_CONVEXQP_PG ||
     internalsolver_options[1]->solverId == SICONOS_CONVEXQP_VI_FPP||
     internalsolver_options[1]->solverId == SICONOS_CONVEXQP_VI_EG)
  {
    tangent_cqp = (ConvexQP *)malloc(sizeof(ConvexQP));
    tangent_cqp->M = splitted_problem->M_tt;
    tangent_cqp->q = (double *) malloc(2* nc * sizeof(double));
    tangent_cqp->ProjectionOnC = &Projection_ConvexQP_FC3D_Disk;
    tangent_cqp->A=NULL;
    tangent_cqp->b= NULL;
    FrictionContactProblem_as_ConvexQP *fc3d_as_cqp= (FrictionContactProblem_as_ConvexQP*)malloc(sizeof(FrictionContactProblem_as_ConvexQP));
    tangent_cqp->env = fc3d_as_cqp ;
    tangent_cqp->size = nc*2;

    /*set the norm of the VI to the norm of problem->q  */
    double norm_q_t = cblas_dnrm2(nc*2 , splitted_problem->q_t , 1);
    tangent_cqp->normConvexQP= norm_q_t;
    tangent_cqp->istheNormConvexQPset=1;

    fc3d_as_cqp->cqp = tangent_cqp;
    fc3d_as_cqp->fc3d = problem;
    fc3d_as_cqp->options = options;
  }
  else
  {
    numerics_error("fc3d_Panagiotopoulos_FixedPoint", "Unknown internal solver for the tangent part.");
  }

 if (internalsolver_options[0]->solverId == SICONOS_LCP_PGS)
 {
   if (verbose > 0)
     printf(" ========================== Call LCP_PGS solver for Friction-Contact 3D problem ==========================\n");
   internalsolver_normal = &lcp_pgs;
 }  
 else if (internalsolver_options[0]->solverId == SICONOS_LCP_CONVEXQP_PG)
 {
   if (verbose > 0)
     printf(" ========================== Call LCP_CONVEX_QP solver for Friction-Contact 3D problem ==========================\n");
   internalsolver_normal = &lcp_ConvexQP_ProjectedGradient;
 }
 else
 {
   numerics_error("fc3d_Panagiotopoulos_FixedPoint", "Unknown internal solver for the normal part.");
 }
   

 
 if (internalsolver_options[1]->solverId == SICONOS_CONVEXQP_PG)
 {
   if (verbose > 0)
     printf(" ========================== Call SICONOS_CONVEX_QP solver for Friction-Contact 3D problem ==========================\n");
   internalsolver_tangent = &convexQP_ProjectedGradient;
 }
 else if  (internalsolver_options[1]->solverId == SICONOS_CONVEXQP_VI_FPP ||
           internalsolver_options[1]->solverId == SICONOS_CONVEXQP_VI_EG )
 {
   if (verbose > 0)
     printf(" ========================== Call SICONOS_CONVEX_VI_FPP solver for Friction-Contact 3D problem ==========================\n");
   internalsolver_tangent = &convexQP_VI_solver;
 }

  int cumul_internal=0;
  //verbose=1;
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;


    fc3d_set_internalsolver_tolerance(problem,options, internalsolver_options[0], error);

    /* ----------------- */
    /* normal resolution */
    /* ----------------- */

    /* compute the rhs of the normal problem */
    cblas_dcopy(nc , splitted_problem->q_n , 1 , normal_lcp_problem->q, 1);
    NM_gemv(1.0, splitted_problem->M_nt, r_t, 1.0, normal_lcp_problem->q);

    (*internalsolver_normal)(normal_lcp_problem, r_n , velocity , info , internalsolver_options[0]);
    cumul_internal += internalsolver_options[0]->iparam[SICONOS_IPARAM_ITER_DONE];

    for (int contact = 0 ; contact < nc; contact ++)
    {
      reaction[contact*3]= r_n[contact];
    }

    fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q,  &error);

    if (error < tolerance)
    {
      hasNotConverged = 0;
    }
    else
    {

      /* ------------------ */
      /* tangent resolution */
      /* ------------------ */

      fc3d_set_internalsolver_tolerance(problem,options,internalsolver_options[1], error);
      /* compute the rhs of the tangent problem */
      cblas_dcopy(2*nc , splitted_problem->q_t, 1 , tangent_cqp->q, 1);
      NM_gemv(1.0, splitted_problem->M_tn, r_n, 1.0, tangent_cqp->q);

      /* Compute the value of the initial value friction threshold*/
      for (int ic = 0 ; ic < nc ; ic++) mu[ic] = fmax(0.0, problem->mu[ic] *  reaction [ic * 3]);

      /* if (verbose>0) */
      /*   printf("norm of mu = %10.5e \n", cblas_dnrm2(nc , mu , 1)); */
      fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q,  &error);

      (*internalsolver_tangent)(tangent_cqp, r_t , velocity , info , internalsolver_options[1]);
      cumul_internal += internalsolver_options[1]->iparam[SICONOS_IPARAM_ITER_DONE];

      for (int contact = 0 ; contact < nc; contact ++)
      {
        reaction[contact*3+1]= r_t[2*contact];
        reaction[contact*3+2]= r_t[2*contact+1];
      }


      /* **** Criterium convergence **** */
      fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q,  &error);


      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                                 reaction, velocity, error, NULL);
      }

      if (error < tolerance) hasNotConverged = 0;
    }



    *info = hasNotConverged;

    if (verbose > 0)
    {
      if (hasNotConverged)
      {
        printf("--------------- FC3D - PFP - Iteration %i error = %14.7e > %10.5e\n", iter, error, tolerance);
      }
      else
      {
        printf("--------------- FC3D - PFP - Iteration %i error = %14.7e < %10.5e\n", iter, error, tolerance);
        printf("--------------- FC3D - PFP - #              Internal iteration = %i\n", cumul_internal);
      }
    }
  }
  internalsolver_options[0]->dWork = NULL;
  mu = NULL;
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

}



void fc3d_pfp_set_options(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] =10.0;

  // Two internal solvers
  assert(options->numberOfInternalSolvers == 2);

  options->internalSolvers[0] = solver_options_create(SICONOS_LCP_PGS);
  options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] =1000;

  options->internalSolvers[1] = solver_options_create(SICONOS_CONVEXQP_VI_FPP);
  options->internalSolvers[1]->iparam[SICONOS_IPARAM_MAX_ITER] =1000;
}
