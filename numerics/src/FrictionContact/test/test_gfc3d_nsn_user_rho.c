/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#include <assert.h>  // for assert
#include <math.h>    // for fabs
#include <stdlib.h>  // for calloc, free, malloc
#include <string.h>  // for memset

#include "FrictionContactProblem.h"  // for FrictionContactProblem, frictionContactProblem_free
#include "FrictionContact_options.h"              // for SICONOS_GLOBAL_FRICTION_3D_NSN_AC, ...
#include "GlobalFrictionContactProblem.h"         // for GlobalFrictionContactProblem, ...
#include "SolverOptions.h"                        // for solver_options_create, ...
#include "fc3d_compute_rho.h"                     // for fc3d_compute_rho
#include "gfc3d_nonsmooth_Newton_AlartCurnier.h"  // for gfc3d_nonsmooth_Newton_AlartCurnier

int main() {
  const char* filename = "./data/GFC3D_OneContact.dat";

  /* Load a 3D global friction contact problem from a data file. */
  GlobalFrictionContactProblem* problem = globalFrictionContact_new_from_filename(filename);
  assert(problem);
  assert(problem->dimension == 3);
  assert(problem->numberOfContacts > 0);

  const int nc = problem->numberOfContacts;
  const int dim = 3 * nc;
  const int n = problem->M->size0;

  /* Work on a copy because the reformulation routine may factorize M in place. */
  GlobalFrictionContactProblem* problem_copy = globalFrictionContact_copy(problem);
  assert(problem_copy);

  /* Reformulate the GFC3D problem as an FC3D problem. The resulting local matrix
   * is W = H^T M^{-1} H and the local vector is q = H^T M^{-1} q_global + b. */
  FrictionContactProblem* local_problem =
      globalFrictionContact_reduce_to_FrictionContact(problem_copy);
  assert(local_problem);
  assert(local_problem->numberOfContacts == nc);
  assert(local_problem->dimension == 3);

  /* Compute a per-contact rho vector from the reformulated problem. */
  double* rho = (double*)calloc((size_t)dim, sizeof(double));
  int rho_info = fc3d_compute_rho(
      local_problem, SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM, 0.0, rho);
  assert(rho_info == 0);
  for (int i = 0; i < dim; ++i) assert(rho[i] > 0.0);

  /* Release the reformulated problem and the copy; the original problem is kept
   * unmodified for the GFC3D solver. */
  frictionContactProblem_free(local_problem);
  globalFrictionContact_free(problem_copy);

  /* Allocate output arrays for the GFC3D solver. */
  double* reaction = (double*)calloc((size_t)dim, sizeof(double));
  double* velocity = (double*)calloc((size_t)dim, sizeof(double));
  double* globalVelocity = (double*)calloc((size_t)n, sizeof(double));

  /* Set up the NSN Alart--Curnier solver and request the user-provided rho
   * strategy. The rho array is heap-allocated so that solver_options_delete can
   * release it together with the other solver data. */
  SolverOptions* options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSN_AC);
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD;
  options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_USER_PROVIDED;
  options->solverData = rho;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;

  int info = -1;
  gfc3d_nonsmooth_Newton_AlartCurnier(problem, reaction, velocity, globalVelocity, &info,
                                      options);
  assert(info == 0);

  /* solver_options_delete frees options->solverData, i.e. the rho vector. */
  solver_options_delete(options);
  globalFrictionContact_free(problem);

  free(reaction);
  free(velocity);
  free(globalVelocity);

  return 0;
}
