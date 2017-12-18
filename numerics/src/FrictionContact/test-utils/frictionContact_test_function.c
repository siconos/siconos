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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SparseMatrix_internal.h"

// avoid a conflict with old csparse.h in case fclib includes it
#define _CS_H

#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"
#include "fc3d_Solvers.h"
#include "gfc3d_Solvers.h"
#include "Friction_cst.h"
#if defined(WITH_FCLIB)
#include <fclib.h>
#include <fclib_interface.h>
#endif
#include "numerics_verbose.h"
#include "SiconosCompat.h"

void frictionContact_test_gams_opts(SN_GAMSparams* GP, int solverId)
{
#ifdef HAVE_GAMS_C_API
  if (solverId == SICONOS_FRICTION_3D_GAMS_PATHVI ||
      solverId == SICONOS_FRICTION_3D_GAMS_LCP_PATHVI ||
      solverId == SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI
      )
  {
    add_GAMS_opt_str(GP, "avi_start", "ray_first", GAMS_OPT_SOLVER);
    add_GAMS_opt_str(GP, "ratio_tester", "expand", GAMS_OPT_SOLVER);
    add_GAMS_opt_double(GP, "expand_eps", 0., GAMS_OPT_SOLVER);
    add_GAMS_opt_bool(GP, "ratio_tester_tfirst", false, GAMS_OPT_SOLVER);
//    add_GAMS_opt_int(GP, "scheduler_decompose", 1, GAMS_OPT_SOLVER);
//    add_GAMS_opt_str(GP, "lemke_factorization_method", "minos_blu", GAMS_OPT_SOLVER);
  }
  else if (solverId == SICONOS_FRICTION_3D_GAMS_PATH ||
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
  else
  {
    fprintf(stderr, "frictionContact_test_gams_opts :: ERROR unknown solverId = %d e.g. solver named %s", solverId, solver_options_id_to_name(solverId));
  }
  add_GAMS_opt_int(GP, "minor_iteration_limit", 100000, GAMS_OPT_SOLVER);
  add_GAMS_opt_int(GP, "major_iteration_limit", 20, GAMS_OPT_SOLVER);
  add_GAMS_opt_double(GP, "expand_delta", 1e-10, GAMS_OPT_SOLVER);
#endif
}

int frictionContact_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  /* numerics_set_verbose(1); */
  info = frictionContact_newFromFile(problem, f);

  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = frictionContact_printInFile(problem, foutput);

  solver_options_print(options);
  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  //int dim = problem->numberOfContacts;
  
  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));

  if (dim == 2)
  {
    info = fc2d_driver(problem,
		       reaction , velocity,
		       options);
  }
  else if (dim == 3)
  {
    info = fc3d_driver(problem,
		       reaction , velocity,
		       options);
  }
  else
  {
    info = 1;
  }
  printf("\n");

  int print_size =10;

  if  (dim * NC >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for (k = 0 ; k < print_size; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf(" ..... \n");
  }
  else
  {
    for (k = 0 ; k < dim * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf("\n");
  }

  /* for (k = 0 ; k < dim * NC; k++) */
  /* { */
  /*   printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]); */
  /* } */
  /* printf("\n"); */

  if (!info)
  {
    printf("test successful, residual = %g\n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful, residual = %g, info = %d, nb iter = %d\n", options->dparam[1], info, options->iparam[1] ? options->iparam[1] : options->iparam[7]);
  }
  free(reaction);
  free(velocity);

  freeFrictionContactProblem(problem);
  fclose(foutput);

  return info;

}


#if defined(WITH_FCLIB)
int frictionContact_test_function_hdf5(const char * path, SolverOptions * options)
{

  int k, info = -1 ;
  /* FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem)); */
  /* info = frictionContact_newFromFile(problem, f); */
  
  FrictionContactProblem* problem = frictionContact_fclib_read(path);
  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = frictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  //int dim = problem->numberOfContacts;
  
  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));

  if (dim == 2)
  {
    info = fc2d_driver(problem,
		       reaction , velocity,
		       options);
  }
  else if (dim == 3)
  {
    info = fc3d_driver(problem,
		       reaction , velocity,
		       options);
  }
  else
  {
    info = 1;
  }
  printf("\n");

  int print_size =10;

  if  (dim * NC >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for (k = 0 ; k < print_size; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf(" ..... \n");
  }
  else
  {
    for (k = 0 ; k < dim * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf("\n");
  }

  for (k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }
  /* for (k = 0 ; k < dim * NC; k++) */
  /* { */
  /*   printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]); */
  /* } */
  /* printf("\n"); */

  if (!info)
  {
    printf("test successful, residual = %g\n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful, residual = %g\n", options->dparam[1]);
  }
  free(reaction);
  free(velocity);

  freeFrictionContactProblem(problem);
  fclose(foutput);

  return info;

}


#endif
