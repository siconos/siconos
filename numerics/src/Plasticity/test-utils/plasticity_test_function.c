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
#include <math.h>    // for isfinite
#include <stdio.h>   // for printf, fclose, fopen, FILE
#include <stdlib.h>  // for calloc, free
#include <time.h>

#include "NonSmoothDrivers.h"     // for plasticity_2d_driver
#include "NumericsFwd.h"          // for SolverOptions
#include "PlasticityProblem.h"    // for plasticity2D_new_from_filename, etc
#include "SiconosBlas.h"
#include "SiconosConfig.h"
#include "SolverOptions.h"  // for SolverOptions
#include "test_utils.h"     // for TestCase
#include "plasticity_test_utils.h"  // for plasticity_test_function

int plasticity_test_function(TestCase* current) {
  int info = -1;

  PlasticityProblem* problem = plasticity2D_new_from_filename(current->filename);
  if (!problem) {
    printf("Failed to load problem from %s\n", current->filename);
    return 1;
  }

  int NC = problem->numberOfCones;
  int dim = problem->dimension;

  double* stress = (double*)calloc(dim * NC, sizeof(double));
  double* plastic_strain_rate = (double*)calloc(dim * NC, sizeof(double));

  long clk_tck = CLOCKS_PER_SEC;

  solver_options_print(current->options);

  clock_t t1 = clock();
  info = plasticity_2d_driver(problem, stress, plastic_strain_rate, current->options);

  int print_size = 10;

  for (int k = 0; k < dim * NC; ++k) {
    info = info == 0 ? !(isfinite(plastic_strain_rate[k]) && isfinite(stress[k])) : info;
  }

  clock_t t2 = clock();

  printf("Norm plastic_strain_rate:  %12.8e\n", cblas_dnrm2(NC * dim, plastic_strain_rate, 1));
  printf("Norm stress:  %12.8e\n", cblas_dnrm2(NC * dim, stress, 1));

  if (dim * NC >= print_size) {
    printf("First values (%i)\n", print_size);
    for (int k = 0; k < print_size; k++) {
      printf("plastic_strain_rate[%i] = %12.8e \t \t stress[%i] = %12.8e\n", 
             k, plastic_strain_rate[k], k, stress[k]);
    }
  } else {
    for (int k = 0; k < dim * NC; k++) {
      printf("plastic_strain_rate[%i] = %12.8e \t \t stress[%i] = %12.8e\n", 
             k, plastic_strain_rate[k], k, stress[k]);
    }
  }
  printf(" ..... \n");

  if (!info)
    printf("test success, residual = %9.2e, info = %d, nb iter = %i\n",
           current->options->dparam[SICONOS_DPARAM_RESIDU], info,
           current->options->iparam[SICONOS_IPARAM_ITER_DONE]);
  else
    printf("test failure, residual = %9.2e, info = %d, nb iter = %i\n",
           current->options->dparam[SICONOS_DPARAM_RESIDU], info,
           current->options->iparam[SICONOS_IPARAM_ITER_DONE]);

  printf("\nsumry: %d  %9.2e  %5i  %10.4f", info, 
         current->options->dparam[SICONOS_DPARAM_RESIDU],
         current->options->iparam[SICONOS_IPARAM_ITER_DONE], 
         (double)(t2 - t1) / (double)clk_tck);
  printf("%3i %5i     %s\n\n", dim, NC, current->filename);

  free(stress);
  free(plastic_strain_rate);
  plasticity2DProblem_free(problem);

  return info;
}
