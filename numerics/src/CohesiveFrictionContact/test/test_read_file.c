/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*!
 * \file test_read_file.c
 * \brief Test that reads a cohesive friction contact problem from file and solves it
 *
 * This test reads a problem from data/cohesive_test_2x2.dat that was written
 * using cohesiveFrictionContact_printInFile, solves it with NSGS solver,
 * and prints the results.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CohesiveFrictionContactProblem.h"
#include "CohesiveFrictionContact_options.h"
#include "NonSmoothGaussSeidel_options.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "NumericsVerbose.h"
#include "SolverOptions.h"
#include "cohesive_friction_3d_driver.h"

#include "numerics_verbose.h"

/* Helper function to print a vector */
static void print_vector(const char* name, double* vec, int size) {
  printf("%s = [", name);
  for (int i = 0; i < size; i++) {
    printf(" %10.6e", vec[i]);
    if (i < size - 1) printf(",");
  }
  printf(" ]\n");
}

static int test_nsgs_on_filename(const char* filename) {
  printf("Reading problem from file: %s\n", filename);
  CohesiveFrictionContactProblem* problem = cohesiveFrictionContact_newFromFilename((char*)filename);
  if (!problem) {
    fprintf(stderr, "ERROR: Failed to read problem from file: %s\n", filename);
    return 1;
  }

  int dim = problem->dimension;
  int nc = problem->numberOfContacts;
  int ncoh = problem->numberOfCohesivePoints;
  int M_size = (nc + ncoh) * dim;

  printf("\nProblem successfully loaded!\n");
  printf("  Dimension: %d\n", dim);
  printf("  Number of contacts: %d\n", nc);
  printf("  Number of cohesive points: %d\n", ncoh);
  printf("  Total problem size: %d\n", M_size);

  /* Display problem details */
  printf("\n");
  cohesiveFrictionContact_display(problem);

  /* Allocate solution vectors */
  double* reaction = (double*)calloc(M_size, sizeof(double));
  double* velocity = (double*)calloc(M_size, sizeof(double));

  /* Create solver options */
  printf("\nSetting up solver...\n");
  SolverOptions* options = solver_options_create(SICONOS_COHESIVE_FRICTION_3D_NSGS);
  if (!options) {
    fprintf(stderr, "ERROR: Failed to create solver options\n");
    cohesiveFrictionContactProblem_free(problem);
    free(reaction);
    free(velocity);
    return 1;
  }

  /* Set solver parameters */
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->iparam[SICONOS_NSGS_FREEZING_CONTACT] = 0;
  options->iparam[SICONOS_NSGS_SHUFFLE] = SICONOS_NSGS_SHUFFLE_TRUE;
  options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
    SICONOS_NSGS_ERROR_EVALUATION_FULL;  
  
  printf("  Solver: NSGS (Non-Smooth Gauss-Seidel)\n");
  printf("  Tolerance: %.2e\n", options->dparam[SICONOS_DPARAM_TOL]);
  printf("  Max iterations: %d\n", options->iparam[SICONOS_IPARAM_MAX_ITER]);

  /* Solve the problem */
  printf("\nSolving...\n");
  numerics_set_verbose(1);
  /* double alpha =100.; */
  /* for (unsigned int k = 0; k < problem->numberOfCohesivePoints; k++) { */
  /*   problem->c_n[k] = alpha*problem->c_n[k] ; */
  /*   problem->c_t[k] = alpha*problem->c_t[k] ; */
  /*   }     */
  /* for (unsigned int k = 0; k < problem->numberOfContacts; k++) { */
  /*   //problem->mu[k] = 1.0; */
  /* }     */

  
  int info = cohesive_friction_3d_driver(problem, reaction, velocity, options);


  

  /* Print results */
  printf("\n=================================================================\n");
  printf("=== Results =====================================================\n");
  printf("=================================================================\n");


    for (int i = 0; i < nc; i++) {
      printf("Contact %d reaction (r):\n", i);
      printf("  r_n  = %12.6e (normal)\n", reaction[0 + dim * i]);
      printf("  r_t1 = %12.6e (tangent 1)\n", reaction[1 + dim * i]);
      printf("  r_t2 = %12.6e (tangent 2)\n", reaction[2 + dim * i]);
      printf("\nContact %d velocity (v):\n", i);
      printf("  v_n  = %12.6e (normal)\n", velocity[0 + dim * i]);
      printf("  v_t1 = %12.6e (tangent 1)\n", velocity[1 + dim * i]);
      printf("  v_t2 = %12.6e (tangent 2)\n", velocity[2 + dim * i]);
      
      /* Verify friction cone condition */
      double r_n = reaction[dim * i];
      double r_t_norm = sqrt(reaction[1 + dim * i] * reaction[1 + dim * i] + 
                             reaction[2 + dim * i] * reaction[2 + dim * i]);
      double mu_r_n = problem->mu[i] * r_n;

      printf("\nVerification:\n");
      printf("  Normal reaction r_n = %.6e\n", r_n);
      printf("  Tangential reaction ||r_t|| = %.6e\n", r_t_norm);
      printf("  Friction bound mu*r_n = %.6e\n", mu_r_n);
      printf("  Friction cone condition: ||r_t|| <= mu*r_n : %s\n",
             (r_t_norm <= mu_r_n + 1e-10) ? "SATISFIED" : "VIOLATED");

      /* Check complementarity */
      double u_n = velocity[dim * i];
      double comp = r_n * u_n;
      printf("  Complementarity r_n * v_n = %.6e (should be ~0)\n", comp);
      printf("\n");
    }
    
    for (int k = 0; k < ncoh; k++) {
      printf("Cohesive point %d reaction (r_coh):\n", k);
      printf("  r_n  = %12.6e (normal)\n", reaction[dim * nc + 0 + dim * k]);
      printf("  r_t1 = %12.6e (tangent 1)\n", reaction[dim * nc + 1 + dim * k]);
      printf("  r_t2 = %12.6e (tangent 2)\n", reaction[dim * nc + 2 + dim * k]);

      printf("\nCohesive point %d displacement rate (u):\n", k);
      printf("  u_n  = %12.6e (normal)\n", velocity[dim * nc + 0 + dim * k]);
      printf("  u_t1 = %12.6e (tangent 1)\n", velocity[dim * nc + 1 + dim * k]);
      printf("  u_t2 = %12.6e (tangent 2)\n", velocity[dim * nc + 2 + dim * k]);
      printf("\n");
    }
  if (info == 0) {
    printf("Solver converged successfully!\n\n");    
  } else {
    printf("Solver failed with error code: %d\n", info);
  }

  printf("Iterations: %d\n", options->iparam[SICONOS_IPARAM_ITER_DONE]);
  printf("Final residual: %.6e\n", options->dparam[SICONOS_DPARAM_RESIDU]);

  /* Cleanup */
  solver_options_delete(options);
  cohesiveFrictionContactProblem_free(problem);
  free(reaction);
  free(velocity);

  return info;
}    


int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  printf("=================================================================\n");
  printf("=== Cohesive Friction 3D Test - Read from File =================\n");
  printf("=================================================================\n\n");

  int info=-1;

  
  /* Read problem from file */
  /* const char* filename_0 = "data/cohesive_test_2x2.dat"; */
  /* info = test_nsgs_on_filename(filename_0); */
  
  /* const char* filename_1 = "data/sphere_2x2_mu0.dat"; */
  /* info += test_nsgs_on_filename(filename_1); */
  
  /* const char* filename_2 = "data/sphere_2x2.dat"; */
  /* info += test_nsgs_on_filename(filename_2); */

  const char* filename_3 = "data/diamond_sphere_sphere_pack.dat";
  info += test_nsgs_on_filename(filename_3);

  printf("\n=================================================================\n");
  printf("=== Test %s ======================================\n",
         (info == 0) ? "PASSED" : "FAILED");
  printf("=================================================================\n");

  return info;
}
