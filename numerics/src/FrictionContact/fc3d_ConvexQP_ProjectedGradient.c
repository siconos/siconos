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

#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FrictionContactProblem_as_ConvexQP.h"
#include "ConvexQP_Solvers.h"
#include "SiconosCompat.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"

#include "SolverOptions.h"
#include "numerics_verbose.h"





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

  int iter=0;
  double error=1e24;

  FrictionContactProblem_as_ConvexQP *fc3d_as_cqp= (FrictionContactProblem_as_ConvexQP*)malloc(sizeof(FrictionContactProblem_as_ConvexQP));
  cqp->env = fc3d_as_cqp ;
  cqp->size = n;

  /*set the norm of the ConvexQP to the norm of problem->q  */
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);
  cqp->normConvexQP= norm_q;
  cqp->istheNormConvexQPset=1;

  fc3d_as_cqp->cqp = cqp;
  fc3d_as_cqp->fc3d = problem;
  fc3d_as_cqp->options = options;
  /* frictionContact_display(fc3d_as_cqp->fc3d); */

  SolverOptions * cqpsolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));

  convexQP_ProjectedGradient_setDefaultSolverOptions(cqpsolver_options);

  int isize = options->iSize;
  int dsize = options->dSize;
  int cqp_isize = cqpsolver_options->iSize;
  int cqp_dsize = cqpsolver_options->dSize;
  if (isize != cqp_isize )
  {
    printf("Warning: options->iSize in fc3d_ConvexQP_FixedPointProjection is not consitent with options->iSize in ConvexQP_FPP\n");
  }
  if (dsize != cqp_dsize )
  {
    printf("Warning: options->iSize in fc3d_ConvexQP_FixedPointProjection is not consitent with options->iSize in ConvexQP_FPP\n");
  }
  int i;
  for (i = 0; i < min(isize,cqp_isize); i++)
  {
    if (options->iparam[i] != 0 )
      cqpsolver_options->iparam[i] = options->iparam[i] ;
  }
  for (i = 0; i <  min(dsize,cqp_dsize); i++)
  {
    if (fabs(options->dparam[i]) >= 1e-24 )
      cqpsolver_options->dparam[i] = options->dparam[i] ;
  }

  convexQP_ProjectedGradient(cqp, reaction, velocity , info , cqpsolver_options);

  /* **** Criterium convergence **** */

  fc3d_Tresca_compute_error(problem, reaction , velocity, options->dparam[0], options, norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  error = cqpsolver_options->dparam[1];
  iter = cqpsolver_options->iparam[7];

  options->dparam[SICONOS_DPARAM_RESIDU] = error;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;


  if (verbose > 0)
  {
    printf("--------------- FC3D - ConvexQP Fixed Point Projection (ConvexQP_FPP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  free(cqp);

  solver_options_delete(cqpsolver_options);
  free(cqpsolver_options);
  free(fc3d_as_cqp);



}


int fc3d_ConvexQP_ProjectedGradient_Cylinder_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the FixedPointProjection Cylinder Solver\n");
  }
  convexQP_ProjectedGradient_setDefaultSolverOptions(options);

  options->solverId = SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder;
  
  return 0;
}


