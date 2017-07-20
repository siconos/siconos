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
#include <string.h>

#include "GlobalFrictionContactProblem_as_VI.h"
#include "GlobalFrictionContactProblem.h"
#include "VariationalInequality_Solvers.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
#include "NumericsMatrix.h"

#include "SolverOptions.h"
#include "numerics_verbose.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
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

  int iter=0;
  double error=1e24;

  GlobalFrictionContactProblem_as_VI *gfc3d_as_vi= (GlobalFrictionContactProblem_as_VI*)malloc(sizeof(GlobalFrictionContactProblem_as_VI));
  vi->env =  gfc3d_as_vi ;
  vi->size = n+m;


  /*Set the norm of the VI to the norm of problem->q  */
  vi->normVI= cblas_dnrm2(n , problem->q , 1);
  vi->istheNormVIset=1;

  gfc3d_as_vi->vi = vi;
  gfc3d_as_vi->gfc3d = problem;
  /* frictionContact_display(fc3d_as_vi->fc3d); */

  SolverOptions * visolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));
  variationalInequality_setDefaultSolverOptions(visolver_options,
                                                SICONOS_VI_EG);
  DEBUG_EXPR(NV_display(reaction,m););
  DEBUG_EXPR(NV_display(globalVelocity,n););
  int isize = options->iSize;
  int dsize = options->dSize;
  int vi_isize = visolver_options->iSize;
  int vi_dsize = visolver_options->dSize;

  if (isize != vi_isize )
  {
    printf("size prolem in gfc3d_VI_FixedPointProjection\n");
  }
  if (dsize != vi_dsize )
  {
    printf("size prolem in gfc3d_VI_FixedPointProjection\n");
  }
  int i;
  for (i = 0; i < isize; i++)
  {
    visolver_options->iparam[i] = options->iparam[i] ;
  }
  for (i = 0; i < dsize; i++)
  {
    visolver_options->dparam[i] = options->dparam[i] ;
  }


  double * z = (double*)malloc((n+m)*sizeof(double));
  double * Fz = (double*)calloc((n+m),sizeof(double));

  memcpy(z, globalVelocity, n * sizeof(double));
  memcpy(&z[n], reaction, m * sizeof(double));
  DEBUG_EXPR(NV_display(z,m+n););
  DEBUG_EXPR(NV_display(Fz,m+n););
  variationalInequality_FixedPointProjection(vi, z, Fz , info , visolver_options);

  memcpy(globalVelocity, z,  n * sizeof(double)  );
  memcpy(reaction, &z[n], m * sizeof(double)  );
  memcpy(velocity, &Fz[m],  m * sizeof(double)  )  ;
  
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
  double norm_q = cblas_dnrm2(n , problem->q , 1);
  gfc3d_compute_error(problem, reaction , velocity, globalVelocity, options->dparam[0], norm_q, &error);

  DEBUG_EXPR(NM_vector_display(reaction,m));
  DEBUG_EXPR(NM_vector_display(velocity,m));
  DEBUG_EXPR(NM_vector_display(globalVelocity,n));

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  error = visolver_options->dparam[1];
  iter = visolver_options->iparam[7];

  options->dparam[1] = error;
  options->iparam[7] = iter;


  if (verbose > 0)
  {
    printf("----------------------------------- FC3D - VI Extra Gradient (VI_EG) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  free(vi);

  solver_options_delete(visolver_options);
  free(visolver_options);
  visolver_options=NULL;
  free(gfc3d_as_vi);

  DEBUG_END("gfc3d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem, ... \n")


}


int gfc3d_VI_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the FixedPointProjection Solver\n");
  }

  variationalInequality_FixedPointProjection_setDefaultSolverOptions(options);
  options->solverId = SICONOS_GLOBAL_FRICTION_3D_VI_FPP;
  return 0;
}
