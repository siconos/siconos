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

#include "FrictionContactProblem_as_VI.h"
#include "VariationalInequality_Solvers.h"
#include "SiconosCompat.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"

#include "SolverOptions.h"
#include "numerics_verbose.h"



void fc3d_VI_FixedPointProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* Dimension of the problem */
  int n = 3 * nc;


  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  //vi.self = &vi;
  vi->F = &Function_VI_FC3D;
  vi->ProjectionOnX = &Projection_VI_FC3D;

  int iter=0;
  double error=1e24;

  FrictionContactProblem_as_VI *fc3d_as_vi= (FrictionContactProblem_as_VI*)malloc(sizeof(FrictionContactProblem_as_VI));
  vi->env =fc3d_as_vi ;
  vi->size =  n;

  /*set the norm of the VI to the norm of problem->q  */
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);
  vi->normVI= norm_q;
  vi->istheNormVIset=1;

  fc3d_as_vi->vi = vi;
  fc3d_as_vi->fc3d = problem;
  /* frictionContact_display(fc3d_as_vi->fc3d); */

  SolverOptions * visolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));
  variationalInequality_setDefaultSolverOptions(visolver_options, SICONOS_VI_FPP);

  int isize = options->iSize;
  int dsize = options->dSize;
  int vi_isize = visolver_options->iSize;
  int vi_dsize = visolver_options->dSize;
  if (isize != vi_isize )
  {
    printf("Warning: options->iSize in fc3d_VI_FixedPointProjection is not consitent with options->iSize in VI_FPP\n");
  }
  if (dsize != vi_dsize )
  {
    printf("Warning: options->iSize in fc3d_VI_FixedPointProjection is not consitent with options->iSize in VI_FPP\n");
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




  variationalInequality_FixedPointProjection(vi, reaction, velocity , info , visolver_options);



  /* **** Criterium convergence **** */
  fc3d_compute_error(problem, reaction , velocity, options->dparam[0], options, norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  error = visolver_options->dparam[1];
  iter = visolver_options->iparam[7];

  options->dparam[SICONOS_DPARAM_RESIDU] = error;
  options->dparam[3] = visolver_options->dparam[SICONOS_VI_EG_DPARAM_RHO];
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;


  if (verbose > 0)
  {
    printf("--------------- FC3D - VI Fixed Point Projection (VI_FPP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  free(vi);

  solver_options_delete(visolver_options);
  free(visolver_options);
  free(fc3d_as_vi);



}


int fc3d_VI_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the FixedPointProjection Solver\n");
  }
  variationalInequality_FixedPointProjection_setDefaultSolverOptions(options);

  options->solverId = SICONOS_FRICTION_3D_VI_FPP;

  return 0;
}

void fc3d_VI_FixedPointProjection_Cylinder(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* Dimension of the problem */
  int n = 3 * nc;


  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  //vi.self = &vi;
  vi->F = &Function_VI_FC3D_Cylinder;
  vi->ProjectionOnX = &Projection_VI_FC3D_Cylinder;

  int iter=0;
  double error=1e24;

  FrictionContactProblem_as_VI *fc3d_as_vi= (FrictionContactProblem_as_VI*)malloc(sizeof(FrictionContactProblem_as_VI));
  vi->env = fc3d_as_vi ;
  vi->size = n;

  /*set the norm of the VI to the norm of problem->q  */
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);
  vi->normVI= norm_q;
  vi->istheNormVIset=1;

  fc3d_as_vi->vi = vi;
  fc3d_as_vi->fc3d = problem;
  fc3d_as_vi->options = options;
  /* frictionContact_display(fc3d_as_vi->fc3d); */

  SolverOptions * visolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));
  variationalInequality_setDefaultSolverOptions(visolver_options,
                                                SICONOS_VI_FPP);

  int isize = options->iSize;
  int dsize = options->dSize;
  int vi_isize = visolver_options->iSize;
  int vi_dsize = visolver_options->dSize;
  if (isize != vi_isize )
  {
    printf("Warning: options->iSize in fc3d_VI_FixedPointProjection is not consitent with options->iSize in VI_FPP\n");
  }
  if (dsize != vi_dsize )
  {
    printf("Warning: options->iSize in fc3d_VI_FixedPointProjection is not consitent with options->iSize in VI_FPP\n");
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

  variationalInequality_FixedPointProjection(vi, reaction, velocity , info , visolver_options);

  /* **** Criterium convergence **** */

  fc3d_Tresca_compute_error(problem, reaction , velocity, options->dparam[0], options, norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  error = visolver_options->dparam[1];
  iter = visolver_options->iparam[7];

  options->dparam[SICONOS_DPARAM_RESIDU] = error;
  options->dparam[3] =  visolver_options->dparam[SICONOS_VI_EG_DPARAM_RHO];
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;


  if (verbose > 0)
  {
    printf("--------------- FC3D - VI Fixed Point Projection (VI_FPP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  free(vi);

  solver_options_delete(visolver_options);
  free(visolver_options);
  free(fc3d_as_vi);



}


int fc3d_VI_FixedPointProjection_Cylinder_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the FixedPointProjection Cylinder Solver\n");
  }
  variationalInequality_FixedPointProjection_setDefaultSolverOptions(options);

  options->solverId = SICONOS_FRICTION_3D_VI_FPP_Cylinder;
  
  return 0;
}


