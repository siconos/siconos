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
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "CohesiveFrictionContactProblem.h"
#include "CohesiveFrictionContact_options.h"
#include "NumericsVector.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"
#include "SiconosBlas.h"
#include "cohesive_friction_3d_projection.h"
#include "cohesive_friction_3d_nsgs.h"
#include "numerics_verbose.h"
#include "numerics_errors.h"
#include "solver_registry.h"

void cohesive_friction_3d_projection_initialize(CohesiveFrictionContactProblem* main_problem) {}

void cohesive_friction_3d_projection_free(CohesiveFrictionContactProblem* main_problem, SolverOptions * options) {}


int cohesive_friction_3d_projection_solve(CohesiveFrictionContactProblem* localproblem, double* reaction,
					  SolverOptions* options) {
  /*  /\* Builds local problem for the current contact *\/ */
  /*   fc3d_projection_update(contact, reaction); */

  /* double* MLocal = localproblem->M->matrix0; */
  /* double* qLocal = localproblem->q; */
  /* double mu_i = localproblem->mu[0]; */
  /* /\* int nLocal = 3; *\/ */

  /* /\* this part is critical for the success of the projection *\/ */
  /* /\*double an = 1./(MLocal[0]);*\/ */
  /* /\*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; *\/ */
  /* /\*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + */
  /*  * MLocal[1*nLocal+2]; *\/ */
  /* /\*   double beta = alpha*alpha - 4*det; *\/ */
  /* /\*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); *\/ */

  /* // double an = 1./(MLocal[0]+mu_i); */
  /* double an = 1. / (MLocal[0]); */

  /* /\* int incx = 1, incy = 1; *\/ */
  /* double worktmp[3]; */
  /* double normUT; */
  /* /\* cblas_dcopy_msan(nLocal , qLocal, incx , worktmp , incy); *\/ */
  /* /\* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, */
  /*  * incx, 1.0, worktmp, incy); *\/ */

  /* for (int i = 0; i < 3; i++) */
  /*   worktmp[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i] + */
  /*                MLocal[i + 1 * 3] * reaction[1] + +MLocal[i + 2 * 3] * reaction[2]; */

  /* normUT = sqrt(worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]); */
  /* reaction[0] -= an * (worktmp[0] + mu_i * normUT); */
  /* reaction[1] -= an * worktmp[1]; */
  /* reaction[2] -= an * worktmp[2]; */

  // projectionOnCone(reaction, mu_i);

  reaction[0] = -localproblem->c_n[0];
  reaction[1] = 0.0;
  reaction[2] = 0.0;

 
  return 0;
}

/* Projection on Cone */
static void cohesive_friction_3d_proj_set_default(SolverOptions* options) { /* No specific defaults */ }

static int cohesive_friction_3d_proj_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int cohesive_friction_3d_proj_solve_wrap(void* problem, double* reaction, double* velocity,
                                     SolverOptions* options) {
  (void)velocity;
  return cohesive_friction_3d_projection_solve((CohesiveFrictionContactProblem*)problem, reaction, options);
}

REGISTER_SOLVER(SICONOS_COHESIVE_FRICTION_3D_PROJECTION,
                "SICONOS_COHESIVE_FRICTION_3D_PROJECTION", "Projection solver (local solver) for cohesion",
                cohesive_friction_3d_proj_init_wrap,
                cohesive_friction_3d_proj_solve_wrap, NULL, NULL,
                cohesive_friction_3d_proj_set_default,
                1000, 1e-14, 1)
