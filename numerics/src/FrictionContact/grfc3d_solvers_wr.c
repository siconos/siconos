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
#include <assert.h>                              // for assert
#include <stdio.h>                               // for printf, fclose, fopen
#include <stdlib.h>                              // for malloc, free, exit
#include "GlobalRollingFrictionContactProblem.h" // for GlobalRollingFrictionContac...
#include "GlobalFrictionContactProblem.h"        // for GlobalFrictionContac...
#include "RollingFrictionContactProblem.h"        // for GlobalFrictionContac...
#include "NumericsFwd.h"                         // for NumericsMatrix, Fric...
#include "NumericsMatrix.h"                      // for NumericsMatrix, NM_gemv
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_d...
#include "rolling_fc_Solvers.h"                        // for fc3d_DeSaxceFixedPoint
#include "numerics_verbose.h"                    // for verbose, numerics_pr...
//#include "gfc3d_compute_error.h"
#include "SolverOptions.h"                       // for SICONOS_DPARAM_TOL

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"                                // for DEBUG_EXPR, DEBUG_P...


#pragma GCC diagnostic ignored "-Wmissing-prototypes"


void  grfc3d_nsgs_wr(GlobalRollingFrictionContactProblem* problem,
                     double *reaction, double *velocity, double* globalVelocity,
                     int *info, SolverOptions* options)
{
  /* verbose=1; */
  DEBUG_BEGIN("grfc3d_nsgs_wr\n");
  NumericsMatrix *H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1/5);
  if(H->size1 > 0)
  {
    // Reformulation
    numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
    RollingFrictionContactProblem* localproblem = globalRollingFrictionContact_reformulation_RollingFrictionContact(problem);
    DEBUG_EXPR(rollingFrictionContact_display(localproblem););

    if(verbose)
    {
      printf("Call to the rfc3d solver ...\n");
    }
    // call nsgs solver for the local problem
    rolling_fc3d_nsgs(localproblem, reaction, velocity, info, options);

    globalRollingFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;
    //gfc3d_compute_error(problem,  reaction, velocity, globalVelocity,  options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);

    rollingFrictionContactProblem_free(localproblem);
  }
  else
  {
    globalRollingFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0 ;
  }
  DEBUG_END("grfc3d_nsgs_wr\n");
}
