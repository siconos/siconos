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

#include <stdio.h>                          // for printf
#include <stdlib.h>                         // for free, malloc
#include "FrictionContactProblem.h"         // for FrictionContactProblem
#include "FrictionContactProblem_as_VI.h"   // for FrictionContactProblem_as_VI
#include "NumericsFwd.h"                    // for VariationalInequality
#include "SiconosBlas.h"                    // for cblas_dnrm2
#include "SolverOptions.h"                  // for SolverOptions, SICONOS_DP...
#include "VI_cst.h"                         // for SICONOS_VI_FPP
#include "VariationalInequality.h"          // for VariationalInequality
#include "VariationalInequality_Solvers.h"  // for variationalInequality_Fix...
#include "fc3d_Solvers.h"                   // for fc3d_VI_FixedPointProjection
#include "fc3d_compute_error.h"             // for fc3d_Tresca_compute_error
#include "numerics_verbose.h"               // for verbose

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

  double error=1e24;

  FrictionContactProblem_as_VI *fc3d_as_vi= (FrictionContactProblem_as_VI*)malloc(sizeof(FrictionContactProblem_as_VI));
  vi->env =fc3d_as_vi ;
  vi->size =  n;

  /*set the norm of the VI to the norm of problem->q  */
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);
  vi->normVI= norm_q;
  vi->istheNormVIset=1;

  fc3d_as_vi->vi = vi;
  fc3d_as_vi->fc3d = problem;
  /* frictionContact_display(fc3d_as_vi->fc3d); */
  variationalInequality_FixedPointProjection(vi, reaction, velocity, info, options);

  /* **** Criterium convergence **** */
  fc3d_compute_error(problem, reaction, velocity, options->dparam[SICONOS_DPARAM_TOL], options, norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  if(verbose > 0)
  {
    printf("--------------- FC3D - VI Fixed Point Projection (VI_FPP) - #Iteration %i Final Residual = %14.7e\n",
           options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  free(vi);
  free(fc3d_as_vi);



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

  double error=1e24;

  FrictionContactProblem_as_VI *fc3d_as_vi= (FrictionContactProblem_as_VI*)malloc(sizeof(FrictionContactProblem_as_VI));
  vi->env = fc3d_as_vi ;
  vi->size = n;

  /*set the norm of the VI to the norm of problem->q  */
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);
  vi->normVI= norm_q;
  vi->istheNormVIset=1;

  fc3d_as_vi->vi = vi;
  fc3d_as_vi->fc3d = problem;
  fc3d_as_vi->options = options;
  /* frictionContact_display(fc3d_as_vi->fc3d); */
  options->solverId = SICONOS_VI_FPP;
  variationalInequality_FixedPointProjection(vi, reaction, velocity, info, options);

  /* **** Criterium convergence **** */

  fc3d_Tresca_compute_error(problem, reaction, velocity, options->dparam[SICONOS_DPARAM_RESIDU], options, norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  if(verbose > 0)
  {
    printf("--------------- FC3D - VI Fixed Point Projection (VI_FPP) - #Iteration %i Final Residual = %14.7e\n",
           options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  free(vi);
  free(fc3d_as_vi);
}



