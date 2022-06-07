/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include <stdio.h>                               // for printf, NULL
#include <stdlib.h>                              // for free, malloc
#include "ConvexQP.h"                            // for ConvexQP
#include "ConvexQP_Solvers.h"                    // for convexQP_ProjectedGr...
#include "ConvexQP_cst.h"                        // for SICONOS_CONVEXQP_PG
#include "FrictionContactProblem.h"              // for FrictionContactProblem
#include "FrictionContactProblem_as_ConvexQP.h"  // for FrictionContactProbl...
#include "NumericsFwd.h"                         // for ConvexQP, SolverOptions
#include "SiconosBlas.h"                         // for cblas_dnrm2
#include "SolverOptions.h"                       // for SolverOptions, SICON...
#include "fc3d_Solvers.h"                        // for fc3d_ConvexQP_Projec...
#include "fc3d_compute_error.h"                  // for fc3d_Tresca_compute_...
#include "numerics_verbose.h"                    // for verbose

void fc3d_ConvexQP_ProjectedGradient_Cylinder(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* Dimension of the problem */
  int n = 3 * nc;


  ConvexQP *cqp = (ConvexQP *)malloc(sizeof(ConvexQP));
  cqp->size=n;
  cqp->M = problem->M;
  cqp->q = problem->q;
  cqp->A = NULL; /* The A matrix is the identity and b is equal to zero */
  cqp->b = NULL;

  cqp->ProjectionOnC = &Projection_ConvexQP_FC3D_Cylinder;

  double error=1e24;

  FrictionContactProblem_as_ConvexQP *fc3d_as_cqp= (FrictionContactProblem_as_ConvexQP*)malloc(sizeof(FrictionContactProblem_as_ConvexQP));
  cqp->env = fc3d_as_cqp ;
  cqp->size = n;

  /*set the norm of the ConvexQP to the norm of problem->q  */
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);
  cqp->normConvexQP= norm_q;
  cqp->istheNormConvexQPset=1;

  fc3d_as_cqp->cqp = cqp;
  fc3d_as_cqp->fc3d = problem;
  fc3d_as_cqp->options = options;
  /* frictionContact_display(fc3d_as_cqp->fc3d); */
  // options->solverId = SICONOS_CONVEXQP_PG;

  // Warning : a new solver options is required here, because dWork
  // is used in convexQP_compute_error_reduced, ProjectionOnC ...
  //
  SolverOptions * cqpsolver_options = solver_options_create(SICONOS_CONVEXQP_PG);
  cqpsolver_options->dparam[SICONOS_DPARAM_TOL] = options->dparam[SICONOS_DPARAM_TOL];
  cqpsolver_options->iparam[SICONOS_IPARAM_MAX_ITER] = options->iparam[SICONOS_IPARAM_MAX_ITER];
  //cqpsolver_options->dWork =  options->dWork;
  convexQP_ProjectedGradient(cqp, reaction, velocity, info, cqpsolver_options);
  //options->solverId = SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER;

  /* **** Criterium convergence **** */
  // Warning: the function below uses options->dWork
  fc3d_Tresca_compute_error(problem, reaction, velocity, options->dparam[SICONOS_DPARAM_TOL], options, norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */
  if(verbose > 0)
  {
    printf("--------------- FC3D - ConvexQP Fixed Point Projection (ConvexQP_FPP) - #Iteration %i Final Residual = %14.7e\n", options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  free(cqp);
  free(fc3d_as_cqp);



}
