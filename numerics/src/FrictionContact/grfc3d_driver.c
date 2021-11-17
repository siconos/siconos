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
#include <assert.h>                        // for assert
#include <stdio.h>                         // for NULL, fprintf, stderr
#include <stdlib.h>                        // for exit, EXIT_FAILURE
#include "Friction_cst.h"                  // for SICONOS_GLOBAL_FRICTION_3D...
#include "GlobalRollingFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for g_rolling_fc3d_driver
#include "NumericsFwd.h"                   // for SolverOptions, GlobalFrict...
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "NumericsMatrix.h"                // 

#include "global_rolling_fc_Solvers.h"     // for grfc3d...
#include "numerics_verbose.h"              // for numerics_printf_verbose
//#include "gfc3d_compute_error.h"
//#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_dscal

#include "siconos_debug.h"                         // for DEBUG_EXPR
#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
#endif

const char* const SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR_STR = "GFC3D_NSGS_WR";

#ifdef WITH_FCLIB
#include "string.h"                  // for strcpy, strcat
#include "fclib_interface.h"         // for frictionContact_fclib_write, fri...
#endif

//#define FCLIB_OUTPUT

#ifdef  FCLIB_OUTPUT
#ifdef WITH_FCLIB
#include "string.h"                  // for strcpy, strcat
#include "fclib_interface.h"         // for frictionContact_fclib_write, fri...
#endif
static int fccounter = -1;
#endif

int g_rolling_fc3d_driver(GlobalRollingFrictionContactProblem* problem, double *reaction, double *velocity,
                  double* globalVelocity,  SolverOptions* options)
{
#ifdef FCLIB_OUTPUT
#ifdef WITH_FCLIB
  fccounter ++;
  int freq_output=1;
  int nc = problem->numberOfContacts;
  if(nc >0)
  {
    if(fccounter % freq_output == 0)
    {
      char fname[256];
      sprintf(fname, "GRFC3D-%.5d-%.5d.hdf5",  (int)nc, fccounter);
      printf("Dump GRFC3D-%.5d-%.5d.hdf5.\n",  (int)nc, fccounter);
      /* printf("ndof = %i.\n", ndof); */
      /* FILE * foutput  =  fopen(fname, "w"); */
      int n = 100;
      char * title = (char *)malloc(n * sizeof(char));
      strcpy(title, "GRFC3 dump in hdf5");
      char * description = (char *)malloc(n * sizeof(char));
      strcpy(description, "Rewriting in hdf5 through siconos of  ");
      strcat(description, fname);
      strcat(description, " in FCLIB format");
      char * mathInfo = (char *)malloc(n * sizeof(char));
      strcpy(mathInfo,  "unknown");
      globalRollingFrictionContact_fclib_write(problem,
                                               title,
                                               description,
                                               mathInfo,
                                               fname);
    }
    /* fclose(foutput); */
  }
#else
  printf("Fclib is not available ...\n");
#endif
#endif




  assert(options->isSet);
  DEBUG_EXPR(NV_display(globalVelocity,problem_ori->M->size0););
  if(verbose > 0)
    solver_options_print(options);

  /* Solver name */
  /*  const char* const  name = options->solverName;*/

  int info = -1 ;

  if(problem->dimension != 5)
    numerics_error("grfc3d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* if there is no contact, we compute directly the global velocity as M^{-1}q */
  int m = problem->H->size1;
  if(m ==0)
  {
    numerics_printf_verbose(1,"---- GFC3D - DRIVER . No contact case. Direct computation of global velocity");
    globalRollingFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    return 0;
  }

  /* Non Smooth Gauss Seidel (NSGS) with reduction*/
  switch(options->solverId)
  {
  case SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGS_WR solver with reformulation into Rolling Friction-Contact 3D problem ==========================\n");
    grfc3d_nsgs_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;
  }
   case SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM:
  {
     grfc3d_IPM(problem, reaction, velocity,
                globalVelocity, &info, options);
     break;

  }
  default:
  {
    fprintf(stderr, "Numerics, grfc3d_driver failed. Unknown solver %d.\n", options->solverId);
    exit(EXIT_FAILURE);

  }
  }

  return info;

}

int grfc3d_checkTrivialCaseGlobal(int n, double* q, double* velocity, double* reaction, double * globalVelocity, SolverOptions* options)
{
  /* norm of vector q */
  /*   double qs = cblas_dnrm2( n , q , 1 ); */
  /*   int i; */
  int info = -1;
  /*   if( qs <= DBL_EPSILON )  */
  /*     { */
  /*       // q norm equal to zero (less than DBL_EPSILON) */
  /*       // -> trivial solution: reaction = 0 and velocity = q */
  /*       for( i = 0 ; i < n ; ++i ) */
  /*  { */
  /*    velocity[i] = q[i]; */
  /*    reaction[i] = 0.; */
  /*  } */
  /*       iparam[2] = 0; */
  /*       iparam[4]= 0; */
  /*       dparam[1] = 0.0; */
  /*       dparam[3] = 0.0; */
  /*       info = 0; */
  /*       if(iparam[1]>0) */
  /*  printf("fc3d driver, norm(q) = 0, trivial solution reaction = 0, velocity = q.\n"); */
  /*     } */
  return info;
}
