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
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "MohrCoulomb2DProblem.h"
#include "NonSmoothDrivers.h"
#include "mc2d_solvers.h"
#include "NumericsMatrix.h"
#include "Plasticity_cst.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"

#include "fc3d_Solvers.h"
#include "frictionContact_test_utils.h"
#include "FrictionContactProblem.h"
#include "Friction_cst.h"

#include "numerics_verbose.h"

static int test_unit(char* filename, SolverOptions* options) {
  numerics_set_verbose(1);

  MohrCoulomb2DProblem* problem = mohrCoulomb2D_new_from_filename(filename);
  // mohrCoulomb2D_display(problem);

  double stress[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double plastic_strain_rate[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  int info = mc2d_driver(problem, stress, plastic_strain_rate, options);
  int NC = problem->numberOfCones;
  int dim = problem->dimension;
  int print_size = 10;
  printf("Norm velocity:  %12.8e\n", cblas_dnrm2(NC * dim, plastic_strain_rate, 1));
  printf("Norm stress:  %12.8e\n", cblas_dnrm2(NC * dim, stress, 1));

  if (dim * NC >= print_size) {
    printf("First values (%i)\n", print_size);
    for (int k = 0; k < print_size; k++) {
      printf("plastic_strain_rate[%i] = %12.8e \t \t stress[%i] = %12.8e\n", k,
             plastic_strain_rate[k], k, stress[k]);
    }
  } else {
    for (int k = 0; k < dim * NC; k++) {
      printf("plastic_strain_rate[%i] = %12.8e \t \t stress[%i] = %12.8e\n", k,
             plastic_strain_rate[k], k, stress[k]);
    }
  }

  mohrCoulomb2DProblem_free(problem);

  /* FrictionContactProblem* FC = frictionContactProblem_new_with_data(3, 3, W, q, mu); */
  /* double r[9] = {0.}; */
  /* double u[9] = {0.}; */

  /* FC->M = NULL; */
  /* FC->q = NULL; */
  /* FC->mu = NULL; */
  /* frictionContactProblem_free(FC); */

  /* FrictionContactProblem* FCdense = frictionContactProblem_new_with_data(3, 3, tmpM, q, mu);
   */
  /* double rdense[9] = {0.}; */
  /* double udense[9] = {0.}; */

  /* FCdense->M = NULL; */
  /* FCdense->q = NULL; */
  /* FCdense->mu = NULL; */
  /* frictionContactProblem_free(FCdense); */
  printf(" ..... \n\n\n\n");
  return info;
}

int main(void) {
  int total_info = 0;

  SolverOptions* options = solver_options_create(MOHR_COULOMB_2D_NSGS);
  /* total_info +=  test_unit("./data/mc2d_example1.dat", options); */
  /* total_info +=  test_unit("./data/mc2d_example1_mu0.dat", options); */


  options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
  
  solver_options_update_internal(
      options, 0, MOHR_COULOMB_2D_ONECONTACT_ProjectionOnConeWithLocalIteration);

  options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
  options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  total_info += test_unit("./data/mc2d_example1.dat", options);
  total_info += test_unit("./data/mc2d_example1_mu0.dat", options);

  return total_info;
}
