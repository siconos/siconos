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
#include <assert.h>                              // for assert
#include <stdio.h>                               // for printf, fclose, fopen
#include <stdlib.h>                              // for malloc, free, exit
#include "FrictionContactProblem.h"              // for FrictionContactProblem
#include "GlobalFrictionContactProblem.h"        // for GlobalFrictionContac...
#include "NumericsFwd.h"                         // for NumericsMatrix, Fric...
#include "NumericsMatrix.h"                      // for NumericsMatrix, NM_gemv
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_d...
#include "fc3d_Solvers.h"                        // for fc3d_DeSaxceFixedPoint
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"  // for fc3d_nonsmooth_Newto...
#include "gfc3d_Solvers.h"                       // for gfc3d_DeSaxceFixedPo...
#include "numerics_verbose.h"                    // for verbose, numerics_pr...
#include "gfc3d_compute_error.h"
#include "SolverOptions.h"                       // for SICONOS_DPARAM_TOL

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"                                // for DEBUG_EXPR, DEBUG_P...
#include "gfc3d_ipm.h"                           // for primalResidual, dualResidual ...
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void  gfc3d_nsgs_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  char *blk_num_name = NULL;
  if (options->solverId == SICONOS_GLOBAL_FRICTION_3D_NSGS_SEP_WR)
  {
    blk_num_name = (char *)malloc(10*sizeof(char));
    char *blk_num_ptr = options->solverData;
    strcpy(blk_num_name, blk_num_ptr);
    free(options->solverData); options->solverData = NULL;
  }

  char *str = (char *) malloc(200);
  if (problem->name)
  {
    strcpy( str, problem->name );
  }
  else
  {
    strcpy( str, "foo_" );
  }
  const char * separators = "/";
  char *strToken = strtok( str, separators );
  for(int i=0; i<5; i++)
  {
    if(strToken != NULL) strToken = strtok ( NULL, separators );
  }
  strToken = strtok ( strToken, "." );
  strcat(strToken, blk_num_name);
  FILE *fileName = fopen("problem_name.res", "w");
  fprintf(fileName, "%s", strToken);
  fclose(fileName);
  free(str);
  free(blk_num_name);

  /* verbose=1; */
  DEBUG_BEGIN("gfc3d_nsgs_wr\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    // call nsgs solver for the local problem
    long clk_tck = CLOCKS_PER_SEC;
    clock_t t1 = clock();
    fc3d_nsgs(localproblem, reaction, velocity, info, options);
    clock_t t2 = clock();
    printf("\nTIME = %10.4f\n",(double)(t2-t1)/(double)clk_tck);


    // // printf("\n\n Compute v of NSGS: \n");
    // FILE *sol_file = fopen("sol_data.res", "r");
    // for (int i=0; i < problem->numberOfContacts*3; i++)
    // {
    //   fscanf(sol_file, "%lf ", reaction+i);
    // }
    // fscanf(sol_file, "\n");
    // fclose(sol_file);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);


    // printf("\n\nNSGS v = "); printBlockVec(globalVelocity, problem->M->size0, 3, 0);
    // printf("\n\n");

    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;
    gfc3d_compute_error(problem,  reaction, velocity, globalVelocity,  options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
  DEBUG_END("gfc3d_nsgs_wr\n");
}


void  gfc3d_admm_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_admm_wr\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_admm(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
  DEBUG_END("gfc3d_admm_wr\n");
}

void  gfc3d_nonsmooth_Newton_AlartCurnier_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););


    numerics_printf("gfc3d_nonsmooth_Newton_AlartCurnier_wr - Call to the fc3d solver ...\n");

    fc3d_nonsmooth_Newton_AlartCurnier(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }

  DEBUG_END("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n")


}

void  gfc3d_nonsmooth_Newton_AlartCurnier_new_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("gfc3d_nonsmooth_Newton_AlartCurnier_new_wr(...)\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););


    numerics_printf("gfc3d_nonsmooth_Newton_AlartCurnier_new_wr - Call to the fc3d solver ...\n");

    fc3d_nonsmooth_Newton_AlartCurnier_new(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }

  DEBUG_END("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n")


}
void  gfc3d_nsgs_velocity_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_nsgs_velocity(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}

void  gfc3d_proximal_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););
    DEBUG_EXPR(frictionContact_display(localproblem););
    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_proximal(localproblem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}

void  gfc3d_DeSaxceFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_DeSaxceFixedPoint(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}

void  gfc3d_TrescaFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_TrescaFixedPoint(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }

}


void gfc3d_ipm_snm_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  verbose = 1;
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_IPM_SNM(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}

void gfc3d_ipm_wr(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int *info, SolverOptions* options)
{
  // verbose = 1;
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ... this make take a while");
    FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_IPM(localproblem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    frictionContactProblem_free(localproblem);
  }
  else
  {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
}


