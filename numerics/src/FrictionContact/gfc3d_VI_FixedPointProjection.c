/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include <stdio.h>                               // for printf
#include <stdlib.h>                              // for free, malloc, calloc
#include <string.h>                              // for memcpy
#include "GlobalFrictionContactProblem.h"        // for GlobalFrictionContac...
#include "GlobalFrictionContactProblem_as_VI.h"  // for GlobalFrictionContac...
#include "NumericsFwd.h"                         // for VariationalInequality
#include "NumericsMatrix.h"                      // for NumericsMatrix
#include "SiconosBlas.h"                         // for cblas_dnrm2
#include "SolverOptions.h"                       // for SolverOptions, SICON...
#include "VariationalInequality.h"               // for VariationalInequality
#include "VariationalInequality_Solvers.h"       // for variationalInequalit...
#include "siconos_debug.h"                               // for DEBUG_EXPR, DEBUG_BEGIN
#include "gfc3d_Solvers.h"                       // for gfc3d_VI_FixedPointP...
#include "gfc3d_compute_error.h"                 // for gfc3d_compute_error
#include "numerics_verbose.h"                    // for verbose

#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

void gfc3d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem,
                                   double *reaction, double *velocity,
                                   double *globalVelocity,  int* info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem, ... \n");
  DEBUG_EXPR(verbose=1;);
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* Dimension of the problem */
  int m = 3 * nc;
  int n = problem->M->size0;

  DEBUG_EXPR(NM_vector_display(reaction, m););
  DEBUG_EXPR(NM_vector_display(velocity, m););
  DEBUG_EXPR(NM_vector_display(globalVelocity, n););

  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  //vi.self = &vi;
  vi->F = &Function_VI_GFC3D;
  vi->ProjectionOnX = &Projection_VI_GFC3D;

  double error=1e24;

  GlobalFrictionContactProblem_as_VI *gfc3d_as_vi= (GlobalFrictionContactProblem_as_VI*)malloc(sizeof(GlobalFrictionContactProblem_as_VI));
  vi->env =  gfc3d_as_vi ;
  vi->size = n+m;


  /*Set the norm of the VI to the norm of problem->q  */
  vi->normVI= cblas_dnrm2(n, problem->q, 1);
  vi->istheNormVIset=1;

  gfc3d_as_vi->vi = vi;
  gfc3d_as_vi->gfc3d = problem;
  /* frictionContact_display(fc3d_as_vi->fc3d); */

  DEBUG_EXPR(NV_display(reaction,m););
  DEBUG_EXPR(NV_display(globalVelocity,n););

  double * z = (double*)malloc((n+m)*sizeof(double));
  double * Fz = (double*)calloc((n+m),sizeof(double));

  memcpy(z, globalVelocity, n * sizeof(double));
  memcpy(&z[n], reaction, m * sizeof(double));
  DEBUG_EXPR(NV_display(z,m+n););
  DEBUG_EXPR(NV_display(Fz,m+n););
  // Warning : tolerance might be updated by VI determine_convergence function !
  variationalInequality_FixedPointProjection(vi, z, Fz, info, options);

  memcpy(globalVelocity, z,  n * sizeof(double));
  memcpy(reaction, &z[n], m * sizeof(double));
  memcpy(velocity, &Fz[m],  m * sizeof(double))  ;

  /* for (int k =0 ;  k< n; k++) */
  /* { */
  /*   globalVelocity[k] = z[k]; */
  /* } */
  /* for (int k =0 ;  k< m; k++) */
  /* { */
  /*   reaction[k] = z[k+n]; */
  /* } */
  /* for (int k =0 ;  k< m; k++) */
  /* { */
  /*   velocity[k] = Fz[k+n]; */
  /* } */

  free(z);
  free(Fz);

  /* **** Criterium convergence **** */
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  double norm_b = cblas_dnrm2(m, problem->b, 1);

  gfc3d_compute_error(problem, reaction, velocity, globalVelocity, options->dparam[SICONOS_DPARAM_TOL],
                      options, norm_q, norm_b, &error);

  DEBUG_EXPR(NM_vector_display(reaction,m));
  DEBUG_EXPR(NM_vector_display(velocity,m));
  DEBUG_EXPR(NM_vector_display(globalVelocity,n));

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  if(verbose > 0)
  {
    printf("--------------- GFC3D - VI Fixed Point Projection (VI_FPP) - #Iteration %i Final Residual = %14.7e\n",
           options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  free(vi);
  free(gfc3d_as_vi);

  DEBUG_END("gfc3d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem, ... \n")


}
