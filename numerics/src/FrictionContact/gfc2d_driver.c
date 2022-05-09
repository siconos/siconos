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
#include <assert.h>                        // for assert
#include <stdio.h>                         // for NULL, fprintf, stderr
#include <stdlib.h>                        // for exit, EXIT_FAILURE
#include "Friction_cst.h"                  // for SICONOS_GLOBAL_FRICTION_3D...
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for gfc3d_driver
#include "NumericsFwd.h"                   // for SolverOptions, GlobalFrict...
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "siconos_debug.h"                         // for DEBUG_EXPR
#include "fc2d_Solvers.h"                  // for fc2d_nsgs

#include "numerics_verbose.h"              // for numerics_printf_verbose
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_dscal

#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
#endif

//const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_WR_STR = "GFC3D_NSGS_WR";



/* static int gfc3d_balancing_check_drift(GlobalFrictionContactProblem* balanced_problem, */
/*                                        GlobalFrictionContactProblem* problem, */
/*                                        double *reaction, double *velocity, */
/*                                        double* globalVelocity,  SolverOptions* options) */
/* { */
/*   if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0) */
/*   { */
/*     size_t nc = problem->numberOfContacts; */
/*     size_t n = problem->M->size0; */
/*     size_t m = 3 * nc; */

/*     double norm_b = cblas_dnrm2(m, balanced_problem->b, 1); */
/*     double norm_q = cblas_dnrm2(n, balanced_problem->q, 1); */
/*     double error_balancing = 1e24; */
/*     double tolerance = options->dparam[SICONOS_DPARAM_TOL]; */
/*     gfc3d_compute_error(balanced_problem,  reaction, velocity, globalVelocity, */
/*                         tolerance, options, */
/*                         norm_q, norm_b,  &error_balancing); */

/*     /\* Come back to original variables *\/ */
/*     gfc3d_balancing_back_to_original_variables(balanced_problem, options, */
/*                                                reaction, velocity, globalVelocity); */

/*     norm_b = cblas_dnrm2(m, problem->b, 1); */
/*     norm_q = cblas_dnrm2(n, problem->q, 1); */
/*     double error =0.0; */
/*     gfc3d_compute_error(problem,  reaction, velocity, globalVelocity, */
/*                         tolerance, options, */
/*                         norm_q, norm_b,  &error); */

/*     numerics_printf_verbose(0,"error with balancing = %8.4e", error_balancing); */
/*     numerics_printf_verbose(0,"error with original = %8.4e", error); */
/*   } */
/*   //else continue */

/*   return 0; */

/* } */

//#define DUMP_PROBLEM
#ifdef DUMP_PROBLEM
static int fccounter = 0;
#endif
//#define DUMP_PROBLEM_IF_INFO
#ifdef DUMP_PROBLEM_IF_INFO
static int fccounter = 0;
#endif

int gfc2d_driver(GlobalFrictionContactProblem* problem, double *reaction, double *velocity,
                 double* globalVelocity,  SolverOptions* options)
{
  assert(options->isSet);
  DEBUG_EXPR(NV_display(globalVelocity,problem_ori->M->size0););

#ifdef DUMP_PROBLEM
  char fname[256];
  int ncc = problem->numberOfContacts;
  sprintf(fname, "gfc2d_granularflowonwall_%.5d_%.5d.dat", ncc, fccounter++);
  printf("Dump %s file\n", fname);

  FILE * foutput  =  fopen(fname, "w");
  globalFrictionContact_printInFile(problem, foutput);
  fclose(foutput);
#endif

  /* verbose=1; */

  if(verbose > 0)
    solver_options_print(options);

  /* Solver name */
  /*  const char* const  name = options->solverName;*/

  int info = -1 ;

  if(problem->dimension != 2)
    numerics_error("gfc2d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* if there is no contact, we compute directly the global velocity as M^{-1}q */
  int m = problem->H->size1;
  if(m ==0)
  {
    numerics_printf_verbose(1,"---- GFC2D - DRIVER . No contact case. Direct computation of global velocity");
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    return 0;
  }

  /* Non Smooth Gauss Seidel (NSGS) */
  switch(options->solverId)
  {
  case SICONOS_FRICTION_2D_NSGS:
  {
    numerics_printf_verbose(1," ========================== Call NSGS solver with reformulation into Friction-Contact 3D problem ==========================\n");
    /* verbose=1; */
    // We compute only if the local problem has contacts
    DEBUG_PRINTF("Number of contacts = %i \n", H->size1/3);
    if(problem->H->size1 > 0)
    {
      // Reformulation
      numerics_printf_verbose(1,"Reformulation info a reduced problem onto local variables ...\n");
      FrictionContactProblem* localproblem = globalFrictionContact_reformulation_FrictionContact(problem);
      DEBUG_EXPR(frictionContact_display(localproblem););

      // call solver for the local problem
      info = fc2d_driver(localproblem, reaction, velocity, options);

      globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

      
      /* /\* Number of contacts *\/ */
      /* int nc = problem->numberOfContacts; */
      /* /\* Dimension of the problem *\/ */
      /* int m = 3 * nc; */
      /* int n = problem->M->size0; */
      /* double norm_q = cblas_dnrm2(n, problem->q, 1); */
      /* double norm_b = cblas_dnrm2(m, problem->b, 1); */
      /* double error; */
      /* gfc3d_compute_error(problem,  reaction, velocity, globalVelocity,  options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error); */


      frictionContactProblem_free(localproblem);
    }
    else
    {
      globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
      info = 0 ;
    }

    break;

  }
  default:
  {
    fprintf(stderr, "Numerics, gfc3d_driver failed. Unknown solver %d.\n", options->solverId);
    exit(EXIT_FAILURE);

  }
  }

  return info;

}

/* int gfc3d_checkTrivialCaseGlobal(int n, double* q, double* velocity, double* reaction, double * globalVelocity, SolverOptions* options) */
/* { */
/*   /\* norm of vector q *\/ */
/*   /\*   double qs = cblas_dnrm2( n , q , 1 ); *\/ */
/*   /\*   int i; *\/ */
/*   int info = -1; */
/*   /\*   if( qs <= DBL_EPSILON )  *\/ */
/*   /\*     { *\/ */
/*   /\*       // q norm equal to zero (less than DBL_EPSILON) *\/ */
/*   /\*       // -> trivial solution: reaction = 0 and velocity = q *\/ */
/*   /\*       for( i = 0 ; i < n ; ++i ) *\/ */
/*   /\*  { *\/ */
/*   /\*    velocity[i] = q[i]; *\/ */
/*   /\*    reaction[i] = 0.; *\/ */
/*   /\*  } *\/ */
/*   /\*       iparam[2] = 0; *\/ */
/*   /\*       iparam[4]= 0; *\/ */
/*   /\*       dparam[1] = 0.0; *\/ */
/*   /\*       dparam[3] = 0.0; *\/ */
/*   /\*       info = 0; *\/ */
/*   /\*       if(iparam[1]>0) *\/ */
/*   /\*  printf("fc3d driver, norm(q) = 0, trivial solution reaction = 0, velocity = q.\n"); *\/ */
/*   /\*     } *\/ */
/*   return info; */
/* } */
