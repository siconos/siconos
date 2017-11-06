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

#include "ConvexQP.h"
#include "ConvexQP_as_VI.h"
#include "VariationalInequality_Solvers.h"
#include "ConvexQP_Solvers.h"
#include "ConvexQP_computeError.h"
#include "SiconosCompat.h"

#include "SolverOptions.h"
#include "numerics_verbose.h"

  
#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

  SolverOptions * visolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));


  if (options->solverId==SICONOS_CONVEXQP_VI_FPP)
  {
    variationalInequality_setDefaultSolverOptions(visolver_options, SICONOS_VI_FPP);
  }
  else if (options->solverId==SICONOS_CONVEXQP_VI_EG)
  {
    variationalInequality_setDefaultSolverOptions(visolver_options, SICONOS_VI_EG);
  }
    
  int isize = options->iSize;
  int dsize = options->dSize;
  int vi_isize = visolver_options->iSize;
  int vi_dsize = visolver_options->dSize;
  if (isize != vi_isize )
  {
    printf("Warning: options->iSize in convexQP_VI_solver is not consitent with options->iSize in VI solver\n");
  }
  if (dsize != vi_dsize )
  {
    printf("Warning: options->iSize in convexQP_VI_solver is not consitent with options->iSize in VI solver\n");
  }
  int i;
  for (i = 0; i < min(isize,vi_isize); i++)
  {
    if (options->iparam[i] != 0 )
      visolver_options->iparam[i] = options->iparam[i] ;
  }
  for (i = 0; i <  min(dsize,vi_dsize); i++)
  {
    if (fabs(options->dparam[i]) >= 1e-24 )
      visolver_options->dparam[i] = options->dparam[i] ;
  }

  
  if (options->solverId==SICONOS_CONVEXQP_VI_FPP)
  {
    variationalInequality_FixedPointProjection(vi, z, w , info , visolver_options);
  }
  else if (options->solverId==SICONOS_CONVEXQP_VI_EG)
  {
    variationalInequality_ExtraGradient(vi, z, w , info , visolver_options);
  }


  /* **** Criterium convergence **** */
  convexQP_computeError(problem, z , w, options->dparam[0], options, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  error = visolver_options->dparam[SICONOS_DPARAM_RESIDU];
  iter = visolver_options->iparam[SICONOS_IPARAM_ITER_DONE];

  options->dparam[SICONOS_DPARAM_RESIDU] = error;
  options->dparam[3] = visolver_options->dparam[SICONOS_VI_EG_DPARAM_RHO];
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;


  if (verbose > 0)
  {

    if (options->solverId==SICONOS_CONVEXQP_VI_FPP)
    {
      printf("--------------- CONVEXQP - VI solver (VI_FPP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
    }
    else if (options->solverId==SICONOS_CONVEXQP_VI_EG)
    {
      printf("--------------- CONVEXQP - VI solver (VI_EG) - #Iteration %i Final Residual = %14.7e\n", iter, error);
    }
    
  }
  free(vi);

  solver_options_delete(visolver_options);
  free(visolver_options);
  free(convexQP_as_vi);



}


int convexQP_VI_solver_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ConvexQP_VI_solver Solver\n");
  }
  variationalInequality_FixedPointProjection_setDefaultSolverOptions(options);

  options->solverId = SICONOS_CONVEXQP_VI_FPP;

  return 0;
}
