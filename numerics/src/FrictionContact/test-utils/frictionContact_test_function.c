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
#include <math.h>                        // for isfinite
#include <stdio.h>                       // for printf, fclose, fopen, FILE
#include <stdlib.h>                      // for calloc, free, rand, srand
#include "FrictionContactProblem.h"      // for frictionContactProblem_free
#include "NonSmoothDrivers.h"            // for fc2d_driver, fc3d_driver
#include "NumericsFwd.h"                 // for SolverOptions, FrictionConta...
#include "SiconosConfig.h"               // for WITH_FCLIB
#include "SolverOptions.h"               // for SolverOptions, solver_option...
#include "frictionContact_test_utils.h"  // for frictionContact_test_function
#include "test_utils.h"                  // for TestCase
#include <time.h>
#include "SiconosConfig.h" // for WITH_FCLIB, HAVE_GAMS_C_API // IWYU pragma: keep
#include "SiconosLapack.h"

// avoid a conflict with old csparse.h in case fclib includes it
#define _CS_H

#if defined(WITH_FCLIB)
#include <fclib_interface.h>             // for globalFrictionContact_fclib_...
#include <time.h>                        // for time
#include "fc3d_solvers_wr.h"             // for fc3d_reformulation_global_pr...
#endif

// --- Extra setup for options when the solver belongs to GAMS family ---
#ifdef HAVE_GAMS_C_API
#include "GAMSlink.h"                    // for SN_GAMSparams
void frictionContact_test_gams_opts(SolverOptions * options)
{
  int solverId = options->solverId;
  if(solverId == SICONOS_FRICTION_3D_GAMS_PATH ||
      solverId == SICONOS_FRICTION_3D_GAMS_LCP_PATH ||
      solverId == SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH||
      solverId == SICONOS_FRICTION_3D_GAMS_PATH ||
      solverId == SICONOS_FRICTION_3D_GAMS_LCP_PATH ||
      solverId == SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH ||
      solverId == SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH ||
      solverId == SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI
    )
  {
    assert(options->solverParameters);
    SN_GAMSparams* GP = (SN_GAMSparams*)options->solverParameters;
    GP->model_dir = strdup(GAMS_MODELS_SOURCE_DIR);
    GP->filename = current->filename;

    if(solverId == SICONOS_FRICTION_3D_GAMS_PATHVI ||
        solverId == SICONOS_FRICTION_3D_GAMS_LCP_PATHVI ||
        solverId == SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI)
    {
      add_GAMS_opt_str(GP, "avi_start", "ray_first", GAMS_OPT_SOLVER);
      add_GAMS_opt_str(GP, "ratio_tester", "expand", GAMS_OPT_SOLVER);
      add_GAMS_opt_double(GP, "expand_eps", 0., GAMS_OPT_SOLVER);
      add_GAMS_opt_bool(GP, "ratio_tester_tfirst", false, GAMS_OPT_SOLVER);
      //    add_GAMS_opt_int(GP, "scheduler_decompose", 1, GAMS_OPT_SOLVER);
      //    add_GAMS_opt_str(GP, "lemke_factorization_method", "minos_blu", GAMS_OPT_SOLVER);
    }
    else if(solverId == SICONOS_FRICTION_3D_GAMS_PATH ||
            solverId == SICONOS_FRICTION_3D_GAMS_LCP_PATH ||
            solverId == SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH
           )
    {
      add_GAMS_opt_int(GP, "linear_model_perturb", 0, GAMS_OPT_SOLVER);
      add_GAMS_opt_double(GP, "proximal_perturbation", 0., GAMS_OPT_SOLVER);
      add_GAMS_opt_double(GP, "proximal_initial_maximum", 0., GAMS_OPT_SOLVER);
      add_GAMS_opt_str(GP, "crash_method", "none", GAMS_OPT_SOLVER);
      add_GAMS_opt_int(GP, "crash_perturb", 0, GAMS_OPT_SOLVER);
      add_GAMS_opt_int(GP, "restart_limit", 0, GAMS_OPT_SOLVER);
      //    add_GAMS_opt_str(GP, "lemke_start", "first", GAMS_OPT_SOLVER);
      //    add_GAMS_opt_int(GP, "output_linear_model", 1, GAMS_OPT_SOLVER);
      //    add_GAMS_opt_int(GP, "output_minor_iterations_frequency", 1, GAMS_OPT_SOLVER);
      //    add_GAMS_opt_int(GP, "output_linear_model", 1, GAMS_OPT_SOLVER);
    }
    add_GAMS_opt_int(GP, "minor_iteration_limit", 100000, GAMS_OPT_SOLVER);
    add_GAMS_opt_int(GP, "major_iteration_limit", 20, GAMS_OPT_SOLVER);
    add_GAMS_opt_double(GP, "expand_delta", 1e-10, GAMS_OPT_SOLVER);
  }
  else
  {
    // nothing, just pass
  }
}
#endif

int frictionContact_test_function(TestCase* current)
{
  int k, info = -1 ;

  FrictionContactProblem * problem  = frictionContact_new_from_filename(current->filename);
  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = frictionContact_printInFile(problem, foutput);

#ifdef WITH_FCLIB
  int global_hdf5_output =0;

  if(global_hdf5_output)
  {
    /* get the current calendar time */
    int stime;
    long ltime;
    ltime = time(NULL);
    stime = (unsigned) ltime/2;
    srand(stime);
    char filename[100];
    sprintf(filename,"gfc3d_%d.hdf5",rand());

    GlobalFrictionContactProblem * gfc3d = fc3d_reformulation_global_problem(problem);

    char * title = "--";
    char * description = "--";
    char * mathInfo = "--";

    globalFrictionContact_fclib_write(gfc3d,
                                      title,
                                      description,
                                      mathInfo,
                                      filename);

  }
#endif

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  //int dim = problem->numberOfContacts;

  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));

  long clk_tck = CLOCKS_PER_SEC;

  solver_options_print(current->options);

// --- Extra setup for options when the solver belongs to GAMS family ---
#ifdef HAVE_GAMS_C_API
  frictionContact_test_gams_opts(current->options);
#endif

  clock_t t1 = clock();

  if(dim == 2)
  {
    info = fc2d_driver(problem,
                       reaction, velocity,
                       current->options);
  }
  else if(dim == 3)
  {
    info = fc3d_driver(problem,
                       reaction, velocity,
                       current->options);
  }
  else
    info = 1;

  for(int k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }


  clock_t t2 = clock();

  int print_size = 10;

  printf("Norm velocity:  %12.8e\n", cblas_dnrm2(NC*dim, velocity, 1));
  printf("Norm reaction:  %12.8e\n", cblas_dnrm2(NC*dim, reaction, 1));

  if(dim * NC >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for(k = 0 ; k < print_size; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
    }
  }
  else
  {
    for(k = 0 ; k < dim * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
    }
  }
  printf(" ..... \n");


  if(!info)
    printf("test: success\n");
  else
    printf("test: failure\n");

  printf("\nsumry: %d  %9.2e  %5i  %10.4f", info, current->options->dparam[SICONOS_DPARAM_RESIDU], current->options->iparam[SICONOS_IPARAM_ITER_DONE], (double)(t2-t1)/(double)clk_tck);
  printf("%3i %5i     %s\n\n", dim, NC, current->filename);

  free(reaction);
  free(velocity);
  frictionContactProblem_free(problem);
  fclose(foutput);

  return info;
}


