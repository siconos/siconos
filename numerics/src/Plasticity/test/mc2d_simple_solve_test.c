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

#include "FrictionContactProblem.h"
#include "Friction_cst.h"
#include "MohrCoulomb2DProblem.h"
#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "Plasticity_cst.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#include "frictionContact_test_utils.h"
#include "mc2d_solvers.h"
#include "numerics_verbose.h"

static int test_unit(char* filename, SolverOptions* options) {


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
  options->dparam[SICONOS_DPARAM_TOL] = 1e-16;


  numerics_set_verbose(0);
  printf("#######\ntest with default options\n");
  total_info += test_unit("./data/mc2d_example1.dat", options);
  total_info += test_unit("./data/mc2d_example1_mu0.dat", options);
  
  numerics_set_verbose(0);
  printf("#######\n test with pure Newton local solver \n");
  solver_options_update_internal(options, 0, MOHR_COULOMB_2D_ONECONE_NSN);
  /* parameters for hybrid solvers */
  options->internalSolvers[0]->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] =
      PLASTICITY_NSN_HYBRID_STRATEGY_NO;
  options->internalSolvers[0]->iparam[PLASTICITY_NSN_FORMULATION] = PLASTICITY_NSN_FORMULATION_ALARTCURNIER_STD;
  options->internalSolvers[0]->iparam[PLASTICITY_NSN_FORMULATION] = PLASTICITY_NSN_FORMULATION_NATURALMAP;
  options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
  total_info += test_unit("./data/mc2d_example1.dat", options);
  /* total_info += test_unit("./data/mc2d_example1_mu0.dat", options); */



  numerics_set_verbose(0);
  printf("#######\ntest with projection on Cone with local iteration solver \n");
  solver_options_update_internal(
      options, 0, MOHR_COULOMB_2D_ONECONE_ProjectionOnConeWithLocalIteration);

  options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
  options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  total_info += test_unit("./data/mc2d_example1.dat", options);
  total_info += test_unit("./data/mc2d_example1_mu0.dat", options);

  return total_info;
}
