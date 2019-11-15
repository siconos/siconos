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

#include <math.h>                           // for fabs
#include <stdio.h>                          // for printf
#include <stdlib.h>                         // for free, malloc
#include "ConvexQP.h"                       // for ConvexQP
#include "ConvexQP_Solvers.h"               // for convexQP_VI_solver, conve...
#include "ConvexQP_as_VI.h"                 // for ConvexQP_as_VI, Function_...
#include "ConvexQP_computeError.h"          // for convexQP_compute_error_re...
#include "ConvexQP_cst.h"                   // for SICONOS_CONVEXQP_VI_FPP
#include "NumericsFwd.h"                    // for SolverOptions, Variationa...
#include "NSSTools.h"                    // for min
#include "SolverOptions.h"                  // for SolverOptions, solver_opt...
#include "VI_cst.h"                         // for SICONOS_VI_EG, SICONOS_VI...
#include "VariationalInequality.h"          // for VariationalInequality
#include "VariationalInequality_Solvers.h"  // for variationalInequality_set...
#include "numerics_verbose.h"               // for numerics_error, verbose
#include "SiconosBlas.h"                          // for cblas_dnrm2

const char* const   SICONOS_CONVEXQP_VI_FPP_STR = "CONVEXQP VI FPP";
const char* const   SICONOS_CONVEXQP_VI_EG_STR = "CONVEXQP VI EG";



void convexQP_VI_solver(ConvexQP* problem, double *z, double *w, int* info, SolverOptions* options);


void convexQP_VI_solver(ConvexQP* problem, double *z, double *w, int* info, SolverOptions* options)
{
  NumericsMatrix* A = problem->A;
  if (A)
  {
    numerics_error("ConvexQP_VI_Solver", "This solver does not support a specific matrix A different from the identity");
  }
  /* Dimension of the problem */
  int n = problem->size;

  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  //vi.self = &vi;
  vi->F = &Function_VI_CQP;
  vi->ProjectionOnX = &Projection_VI_CQP;

  int iter=0;
  double error=1e24;

  ConvexQP_as_VI *convexQP_as_vi= (ConvexQP_as_VI*)malloc(sizeof(ConvexQP_as_VI));
  vi->env =convexQP_as_vi ;
  vi->size =  n;

  /*set the norm of the VI to the norm of problem->q  */
  double norm_q = cblas_dnrm2(n, problem->q , 1);
  vi->normVI= norm_q;
  vi->istheNormVIset=1;

  convexQP_as_vi->vi = vi;
  convexQP_as_vi->cqp = problem;

  if (options->solverId==SICONOS_CONVEXQP_VI_FPP)
  {
    variationalInequality_FixedPointProjection(vi, z, w , info , options);
  }
  else if (options->solverId==SICONOS_CONVEXQP_VI_EG)
  {
    variationalInequality_ExtraGradient(vi, z, w , info , options);
  }


  /* **** Criterium convergence **** */
  convexQP_compute_error_reduced(problem, z , w, options->dparam[SICONOS_DPARAM_TOL], options, norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  if (verbose > 0)
  {

    if (options->solverId==SICONOS_CONVEXQP_VI_FPP)
    {
      printf("--------------- CONVEXQP - VI solver (VI_FPP) - #Iteration %i Final Residual = %14.7e\n",
             options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
    }
    else if (options->solverId==SICONOS_CONVEXQP_VI_EG)
    {
      printf("--------------- CONVEXQP - VI solver (VI_EG) - #Iteration %i Final Residual = %14.7e\n",
                          options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
    }
    
  }
  free(vi);

  free(convexQP_as_vi);



}

