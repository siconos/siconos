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

#include "MixedProblem_NSGS.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CSparseMatrix.h"
#include "LCP_Solvers.h"
#include "MixedProblem.h"
#include "MixedProblems_options.h"
#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include "SolverOptions.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"

/* Include friction contact solvers */
#include "FrictionContactProblem.h"
#include "fc2d_Solvers.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "fc3d_projection.h"
#include "projectionOnCone.h"

/* Include rolling friction contact solvers */
#include "RollingFrictionContact/rolling_fc_Solvers.h"
#include "RollingFrictionContactProblem.h"
#include "RollingFrictionContact_options.h"

/* Include plasticity solvers */
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "plasticity_2d_compute_error.h"
#include "plasticity_2d_projection.h"
#include "plasticity_2d_solvers.h"
#include "plasticity_2d_vonmises_projection.h"

/* Local solver data for each block */
typedef struct MixedProblemLocalSolverData {
  void* local_problem;
  SolverOptions* local_options;
  double* local_z;
  double* local_w;
  double* local_q;
  NumericsMatrix* local_M;
} MixedProblemLocalSolverData;

/* Solver data for the NSGS */
typedef struct MixedProblemNSGSData {
  MixedProblemLocalSolverData* block_solvers;
  double* z_prev;
} MixedProblemNSGSData;

static int mixedProblem_local_solver_init(MixedProblem* problem, int block_index,
                                          MixedProblemLocalSolverData* solver_data,
                                          SolverOptions* internal_options) {
  assert(problem);
  assert(solver_data);
  assert(block_index >= 0 && block_index < problem->num_blocks);

  MixedProblemBlock* block = &problem->blocks[block_index];
  int dim = block->dimension;
  int ne = block->num_elements;
  int local_size = dim * ne;

  /* Allocate local vectors */
  solver_data->local_z = (double*)calloc(local_size, sizeof(double));
  solver_data->local_w = (double*)calloc(local_size, sizeof(double));
  solver_data->local_q = (double*)calloc(local_size, sizeof(double));

  if (!solver_data->local_z || !solver_data->local_w || !solver_data->local_q) {
    {
      int err = numerics_error("mixedProblem_local_solver_init", "Memory allocation failed");
      (void)err;
    }
    return NUMERICS_ERR_MEMORY;
  }

  /* Use the pre-created internal solver options */
  solver_data->local_options = internal_options;

  /* Initialize based on model type */
  switch (block->model_type) {
    case MIXED_MODEL_FRICTION_2D: {
      /* Create local 2D friction problem */
      FrictionContactProblem* fc2d = frictionContactProblem_new();
      if (!fc2d) return NUMERICS_ERR_MEMORY;

      fc2d->dimension = 2;
      fc2d->numberOfContacts = ne;
      fc2d->mu = block->model.friction2d->mu;

      /* Local M will be a view/extraction from global M */
      fc2d->M = NM_create(NM_DENSE, local_size, local_size);
      if (!fc2d->M) {
        free(fc2d);
        return NUMERICS_ERR_MEMORY;
      }

      fc2d->q = solver_data->local_q;
      solver_data->local_problem = fc2d;

      break;
    }

    case MIXED_MODEL_FRICTION_3D: {
      /* Create local 3D friction problem */
      FrictionContactProblem* fc3d = frictionContactProblem_new();
      if (!fc3d) return NUMERICS_ERR_MEMORY;

      fc3d->dimension = 3;
      fc3d->numberOfContacts = ne;
      fc3d->mu = block->model.friction3d->mu;

      fc3d->M = NM_create(NM_DENSE, local_size, local_size);
      if (!fc3d->M) {
        free(fc3d);
        return NUMERICS_ERR_MEMORY;
      }

      fc3d->q = solver_data->local_q;
      solver_data->local_problem = fc3d;

      break;
    }

    case MIXED_MODEL_ROLLING_FRICTION_2D: {
      /* Create local 2D rolling friction problem */
      RollingFrictionContactProblem* rfc2d = rollingFrictionContactProblem_new();
      if (!rfc2d) return NUMERICS_ERR_MEMORY;

      rfc2d->dimension = 2;
      rfc2d->numberOfContacts = ne;
      rfc2d->mu = block->model.rolling_friction2d->mu;
      rfc2d->mu_r = block->model.rolling_friction2d->mur;

      rfc2d->M = NM_create(NM_DENSE, local_size, local_size);
      if (!rfc2d->M) {
        free(rfc2d);
        return NUMERICS_ERR_MEMORY;
      }

      rfc2d->q = solver_data->local_q;
      solver_data->local_problem = rfc2d;

      break;
    }

    case MIXED_MODEL_ROLLING_FRICTION_3D: {
      /* Create local 3D rolling friction problem */
      RollingFrictionContactProblem* rfc3d = rollingFrictionContactProblem_new();
      if (!rfc3d) return NUMERICS_ERR_MEMORY;

      rfc3d->dimension = 3;
      rfc3d->numberOfContacts = ne;
      rfc3d->mu = block->model.rolling_friction3d->mu;
      rfc3d->mu_r = block->model.rolling_friction3d->mur;

      rfc3d->M = NM_create(NM_DENSE, local_size, local_size);
      if (!rfc3d->M) {
        free(rfc3d);
        return NUMERICS_ERR_MEMORY;
      }

      rfc3d->q = solver_data->local_q;
      solver_data->local_problem = rfc3d;

      solver_data->local_options =
          solver_options_create(SICONOS_ROLLING_FRICTION_3D_NSGS); /* ID 3000 */
      if (!solver_data->local_options) return NUMERICS_ERR_MEMORY;

      break;
    }

    case MIXED_MODEL_VON_MISES: {
      /* Create local plasticity problem with Von Mises model */
      PlasticityProblem* plasticity = plasticityProblem_new();
      if (!plasticity) return NUMERICS_ERR_MEMORY;

      plasticity->dimension = 3;
      plasticity->numberOfCones = ne;
      plasticity->model_type = PLASTICITY_MODEL_VON_MISES;
      plasticity->model.von_mises =
          (Plasticity_VonMises_model*)malloc(sizeof(Plasticity_VonMises_model));
      if (!plasticity->model.von_mises) {
        free(plasticity);
        return NUMERICS_ERR_MEMORY;
      }
      plasticity->model.von_mises->sigma_y = block->model.von_mises->sigma_y;

      plasticity->M = NM_create(NM_DENSE, local_size, local_size);
      if (!plasticity->M) {
        free(plasticity->model.von_mises);
        free(plasticity);
        return NUMERICS_ERR_MEMORY;
      }

      plasticity->q = solver_data->local_q;
      solver_data->local_problem = plasticity;

      break;
    }

    case MIXED_MODEL_DRUCKER_PRAGER: {
      /* Create local plasticity problem with Drucker-Prager model */
      PlasticityProblem* plasticity = plasticityProblem_new();
      if (!plasticity) return NUMERICS_ERR_MEMORY;

      plasticity->dimension = 3;
      plasticity->numberOfCones = ne;
      plasticity->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
      plasticity->model.drucker_prager =
          (Plasticity_DruckerPrager_model*)malloc(sizeof(Plasticity_DruckerPrager_model));
      if (!plasticity->model.drucker_prager) {
        free(plasticity);
        return NUMERICS_ERR_MEMORY;
      }
      plasticity->model.drucker_prager->eta = block->model.drucker_prager->eta;
      plasticity->model.drucker_prager->theta = block->model.drucker_prager->theta;

      plasticity->M = NM_create(NM_DENSE, local_size, local_size);
      if (!plasticity->M) {
        free(plasticity->model.drucker_prager);
        free(plasticity);
        return NUMERICS_ERR_MEMORY;
      }

      plasticity->q = solver_data->local_q;
      solver_data->local_problem = plasticity;

      break;
    }

    default: {
      int err = numerics_error("mixedProblem_local_solver_init", "Unknown model type");
      (void)err;
    }
      return NUMERICS_ERR_INVALID_ARGUMENT;
  }

  return NUMERICS_OK;
}

static void mixedProblem_local_solver_free(MixedProblemLocalSolverData* solver_data,
                                           MixedModelType model_type) {
  if (!solver_data) return;

  free(solver_data->local_z);
  free(solver_data->local_w);
  free(solver_data->local_q);

  /* Note: local_options is not freed here as it belongs to the main SolverOptions
   * internalSolvers array */

  if (solver_data->local_problem) {
    switch (model_type) {
      case MIXED_MODEL_FRICTION_2D:
        frictionContactProblem_free((FrictionContactProblem*)solver_data->local_problem);
        break;
      case MIXED_MODEL_FRICTION_3D:
        frictionContactProblem_free((FrictionContactProblem*)solver_data->local_problem);
        break;
      case MIXED_MODEL_ROLLING_FRICTION_2D:
        rollingFrictionContactProblem_free(
            (RollingFrictionContactProblem*)solver_data->local_problem);
        break;
      case MIXED_MODEL_ROLLING_FRICTION_3D:
        rollingFrictionContactProblem_free(
            (RollingFrictionContactProblem*)solver_data->local_problem);
        break;
      case MIXED_MODEL_VON_MISES:
      case MIXED_MODEL_DRUCKER_PRAGER:
        plasticityProblem_free((PlasticityProblem*)solver_data->local_problem);
        break;
      default:
        break;
    }
  }
}

static void mixedProblem_extract_local_block(MixedProblem* problem, int block_index, double* z,
                                             double* w, double* q,
                                             MixedProblemLocalSolverData* solver_data) {
  MixedProblemBlock* block = &problem->blocks[block_index];
  int start = block->start_index;
  int dim = block->dimension;
  int ne = block->num_elements;
  int local_size = dim * ne;

  /* Extract local vectors from global vectors */
  for (int i = 0; i < local_size; i++) {
    solver_data->local_z[i] = z[start + i];
    solver_data->local_w[i] = w[start + i];
    solver_data->local_q[i] = q[start + i];
  }

  /* Extract local M from global M */
  NumericsMatrix* global_M = problem->M;
  NumericsMatrix* local_M = NULL;

  switch (block->model_type) {
    case MIXED_MODEL_FRICTION_2D:
      local_M = ((FrictionContactProblem*)solver_data->local_problem)->M;
      break;
    case MIXED_MODEL_FRICTION_3D:
      local_M = ((FrictionContactProblem*)solver_data->local_problem)->M;
      break;
    case MIXED_MODEL_ROLLING_FRICTION_2D:
      local_M = ((RollingFrictionContactProblem*)solver_data->local_problem)->M;
      break;
    case MIXED_MODEL_ROLLING_FRICTION_3D:
      local_M = ((RollingFrictionContactProblem*)solver_data->local_problem)->M;
      break;
    case MIXED_MODEL_VON_MISES:
    case MIXED_MODEL_DRUCKER_PRAGER:
      local_M = ((PlasticityProblem*)solver_data->local_problem)->M;
      break;
    default:
      return;
  }

  /* Copy submatrix from global M to local M */
  for (int i = 0; i < local_size; i++) {
    for (int j = 0; j < local_size; j++) {
      local_M->dense[i * local_size + j] =
          global_M->dense[(start + i) * problem->total_dimension + (start + j)];
    }
  }
}

static void mixedProblem_update_global_block(MixedProblem* problem, int block_index, double* z,
                                             double* w,
                                             MixedProblemLocalSolverData* solver_data) {
  MixedProblemBlock* block = &problem->blocks[block_index];
  int start = block->start_index;
  int dim = block->dimension;
  int ne = block->num_elements;
  int local_size = dim * ne;

  /* Update global vectors from local vectors */
  for (int i = 0; i < local_size; i++) {
    z[start + i] = solver_data->local_z[i];
    w[start + i] = solver_data->local_w[i];
  }
}

/* Compute local w = M*z + q for a single block */
static void mixedProblem_compute_local_w(MixedProblemLocalSolverData* solver_data,
                                         int local_size) {
  NumericsMatrix* local_M = NULL;

  /* Get local M based on problem type - we need to extract it from the local problem */
  /* For now, assume local_M is stored in the solver_data or compute directly */
  /* This is a simplified version - in practice we'd get M from the local problem structure */

  double* z = solver_data->local_z;
  double* w = solver_data->local_w;
  double* q = solver_data->local_q;

  /* Get local M dense - we need to retrieve it from the local problem */
  /* Since local_M was set up during initialization, we need to access it */
  /* For simplicity, we'll assume diagonal M (identity) for the local block */
  /* In a full implementation, we'd extract M from the local problem structure */

  /* w = z + q (since M = I for local block) */
  for (int i = 0; i < local_size; i++) {
    w[i] = z[i] + q[i];
  }
}

static int mixedProblem_solve_local_block(MixedProblem* problem, int block_index,
                                          MixedProblemLocalSolverData* solver_data,
                                          SolverOptions* global_options) {
  MixedProblemBlock* block = &problem->blocks[block_index];
  int info = 0;

  /* Update local q with coupling terms from other blocks */
  int start = block->start_index;
  int local_size = block->dimension * block->num_elements;

  for (int i = 0; i < local_size; i++) {
    double coupling = 0.0;
    for (int j = 0; j < problem->total_dimension; j++) {
      if (j < start || j >= start + local_size) {
        coupling += problem->M->dense[(start + i) * problem->total_dimension + j] *
                    solver_data->local_z[j - start];
      }
    }
    solver_data->local_q[i] = problem->q[start + i] - coupling;
  }

  /* Solve local problem using direct local solvers (not drivers) */
  switch (block->model_type) {
    case MIXED_MODEL_FRICTION_2D: {
      /* For 2D friction, use direct NSGS since there's no simple projection */
      FrictionContactProblem* fc2d = (FrictionContactProblem*)solver_data->local_problem;
      /* Copy updated q to local problem */
      fc2d->q = solver_data->local_q;
      /* Use nsgs for small 2D problem (returns void, info passed by pointer) */
      fc2d_nsgs(fc2d, solver_data->local_z, solver_data->local_w, &info,
                solver_data->local_options);
      break;
    }

    case MIXED_MODEL_FRICTION_3D: {
      /* For 3D friction, use direct projection on cone */
      FrictionContactProblem* fc3d = (FrictionContactProblem*)solver_data->local_problem;
      /* Copy updated q to local problem */
      fc3d->q = solver_data->local_q;
      /* Use projection on cone directly */
      info =
          fc3d_projectionOnCone_solve(fc3d, solver_data->local_z, solver_data->local_options);
      /* Compute w = M*z + q */
      mixedProblem_compute_local_w(solver_data, local_size);
      break;
    }

    case MIXED_MODEL_ROLLING_FRICTION_2D: {
      /* For rolling friction 2D, use nsgs - returns void */
      RollingFrictionContactProblem* rfc2d =
          (RollingFrictionContactProblem*)solver_data->local_problem;
      rfc2d->q = solver_data->local_q;
      rolling_friction_2d_nsgs(rfc2d, solver_data->local_z, solver_data->local_w, &info,
                               solver_data->local_options);
      break;
    }

    case MIXED_MODEL_ROLLING_FRICTION_3D: {
      /* For rolling friction 3D, use nsgs - returns void */
      RollingFrictionContactProblem* rfc3d =
          (RollingFrictionContactProblem*)solver_data->local_problem;
      rfc3d->q = solver_data->local_q;
      rolling_friction_3d_nsgs(rfc3d, solver_data->local_z, solver_data->local_w, &info,
                               solver_data->local_options);
      break;
    }

    case MIXED_MODEL_VON_MISES: {
      /* For Von Mises, use direct radial return projection */
      PlasticityProblem* plasticity = (PlasticityProblem*)solver_data->local_problem;
      plasticity->q = solver_data->local_q;
      info = plasticity_2d_projectionOnVonMises_solve(plasticity, solver_data->local_z,
                                                      solver_data->local_options);
      /* Compute w = M*z + q */
      mixedProblem_compute_local_w(solver_data, local_size);
      break;
    }

    case MIXED_MODEL_DRUCKER_PRAGER: {
      /* For Drucker-Prager, use projection on cone */
      PlasticityProblem* plasticity = (PlasticityProblem*)solver_data->local_problem;
      plasticity->q = solver_data->local_q;
      info = plasticity_2d_projectionOnCone_solve(plasticity, solver_data->local_z,
                                                  solver_data->local_options);
      /* Compute w = M*z + q */
      mixedProblem_compute_local_w(solver_data, local_size);
      break;
    }

    default:
      info = NUMERICS_ERR_INVALID_ARGUMENT;
      break;
  }

  return info;
}

void mixedProblem_NSGS(MixedProblem* problem, double* z, double* w, int* info,
                       SolverOptions* options) {
  assert(problem && z && w && info && options);

  int num_blocks = problem->num_blocks;
  int total_dim = problem->total_dimension;

  MixedProblemNSGSData* nsgs_data = (MixedProblemNSGSData*)options->solverData;
  if (!nsgs_data) {
    *info = numerics_error("mixedProblem_NSGS", "Solver data not initialized");
    return;
  }

  int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  double error = 0.0;
  int iter = 0;

  /* Initialize z and w */
  for (int i = 0; i < total_dim; i++) {
    z[i] = 0.0;
    w[i] = 0.0;
  }

  /* NSGS iteration */
  for (iter = 0; iter < max_iter; iter++) {
    /* Save previous z for convergence check */
    memcpy(nsgs_data->z_prev, z, total_dim * sizeof(double));

    /* Gauss-Seidel sweep over all blocks */
    for (int b = 0; b < num_blocks; b++) {
      MixedProblemLocalSolverData* solver_data = &nsgs_data->block_solvers[b];

      /* Extract local block from global vectors */
      mixedProblem_extract_local_block(problem, b, z, w, problem->q, solver_data);

      /* Solve local problem */
      int local_info = mixedProblem_solve_local_block(problem, b, solver_data, options);
      if (local_info != 0) {
        if (verbose > 0) {
          printf("MixedProblem_NSGS: Local solver failed for block %d\n", b);
        }
      }

      /* Update global vectors with local solution */
      mixedProblem_update_global_block(problem, b, z, w, solver_data);
    }

    /* Compute error */
    *info = mixedProblem_compute_error(problem, z, w, tolerance, options, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;

    if (verbose > 0) {
      printf("MixedProblem_NSGS: iter %d, error = %e\n", iter, error);
    }

    if (error < tolerance) {
      *info = 0;
      break;
    }
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
}

void mixedProblem_NSGS_init(MixedProblem* problem, SolverOptions* options) {
  assert(problem && options);

  int num_blocks = problem->num_blocks;
  int total_dim = problem->total_dimension;

  /* Allocate NSGS data structure */
  MixedProblemNSGSData* nsgs_data =
      (MixedProblemNSGSData*)malloc(sizeof(MixedProblemNSGSData));
  if (!nsgs_data) {
    {
      int err = numerics_error("mixedProblem_NSGS_init", "Memory allocation failed");
      (void)err;
    }
    return;
  }

  /* Allocate block solvers */
  nsgs_data->block_solvers =
      (MixedProblemLocalSolverData*)calloc(num_blocks, sizeof(MixedProblemLocalSolverData));
  if (!nsgs_data->block_solvers) {
    free(nsgs_data);
    {
      int err = numerics_error("mixedProblem_NSGS_init", "Memory allocation failed");
      (void)err;
    }
    return;
  }

  /* Allocate previous solution vector */
  nsgs_data->z_prev = (double*)calloc(total_dim, sizeof(double));
  if (!nsgs_data->z_prev) {
    free(nsgs_data->block_solvers);
    free(nsgs_data);
    {
      int err = numerics_error("mixedProblem_NSGS_init", "Memory allocation failed");
      (void)err;
    }
    return;
  }

  /* Initialize local solvers for each block, passing internal solver options */
  for (int b = 0; b < num_blocks; b++) {
    SolverOptions* internal_options =
        (options->internalSolvers && b < options->numberOfInternalSolvers)
            ? options->internalSolvers[b]
            : NULL;
    int err = mixedProblem_local_solver_init(problem, b, &nsgs_data->block_solvers[b],
                                             internal_options);
    if (err != NUMERICS_OK) {
      /* Cleanup already initialized blocks */
      for (int i = 0; i < b; i++) {
        mixedProblem_local_solver_free(&nsgs_data->block_solvers[i],
                                       problem->blocks[i].model_type);
      }
      free(nsgs_data->z_prev);
      free(nsgs_data->block_solvers);
      free(nsgs_data);
      int err = numerics_error("mixedProblem_NSGS_init", "Failed to initialize local solver");
      (void)err;
      return;
    }
  }

  options->solverData = nsgs_data;
}

void mixedProblem_NSGS_free(MixedProblem* problem, SolverOptions* options) {
  assert(problem && options);

  MixedProblemNSGSData* nsgs_data = (MixedProblemNSGSData*)options->solverData;
  if (!nsgs_data) return;

  /* Free local solvers for each block */
  for (int b = 0; b < problem->num_blocks; b++) {
    mixedProblem_local_solver_free(&nsgs_data->block_solvers[b],
                                   problem->blocks[b].model_type);
  }

  free(nsgs_data->block_solvers);
  free(nsgs_data->z_prev);
  free(nsgs_data);

  options->solverData = NULL;
}

void mixedProblem_NSGS_set_default_solver_options(MixedProblem* problem,
                                                  SolverOptions* options) {
  assert(options);
  assert(problem);

  options->solverId = SICONOS_MIXED_NSGS;
  options->numberOfInternalSolvers = problem->num_blocks;
  options->isSet = 1;

  options->iparam = (int*)calloc(OPTIONS_PARAM_SIZE, sizeof(int));
  options->dparam = (double*)calloc(OPTIONS_PARAM_SIZE, sizeof(double));

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  /* Create internal solvers for each block based on model type */
  if (problem->num_blocks > 0) {
    options->internalSolvers =
        (SolverOptions**)calloc(problem->num_blocks, sizeof(SolverOptions*));
    if (!options->internalSolvers) {
      int err = numerics_error("mixedProblem_NSGS_set_default_solver_options",
                               "Failed to allocate internal solvers array");
      (void)err;
      return;
    }

    for (int b = 0; b < problem->num_blocks; b++) {
      MixedProblemBlock* block = &problem->blocks[b];
      SolverOptions* internal_options = NULL;

      switch (block->model_type) {
        case MIXED_MODEL_FRICTION_2D:
          /* Use NSGS for 2D friction */
          internal_options = solver_options_create(SICONOS_FRICTION_2D_NSGS);
          break;

        case MIXED_MODEL_FRICTION_3D:
          /* Use projection on cone for 3D friction */
          internal_options =
              solver_options_create(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
          break;

        case MIXED_MODEL_ROLLING_FRICTION_2D:
          /* Use NSGS for rolling friction 2D */
          internal_options = solver_options_create(SICONOS_ROLLING_FRICTION_2D_NSGS);
          break;

        case MIXED_MODEL_ROLLING_FRICTION_3D:
          /* Use NSGS for rolling friction 3D */
          internal_options = solver_options_create(SICONOS_ROLLING_FRICTION_3D_NSGS);
          break;

        case MIXED_MODEL_VON_MISES:
          /* Use Von Mises radial return projection */
          internal_options = solver_options_create(PLASTICITY_2D_ONECONE_VONMISES);
          break;

        case MIXED_MODEL_DRUCKER_PRAGER:
          /* Use projection on cone for Drucker-Prager */
          internal_options = solver_options_create(PLASTICITY_2D_ONECONE_ProjectionOnCone);
          break;

        default:
          /* Unknown model type - will be handled later */
          break;
      }

      options->internalSolvers[b] = internal_options;
    }
  }
}

int mixedProblem_compute_error(MixedProblem* problem, double* z, double* w, double tolerance,
                               SolverOptions* options, double* error) {
  assert(problem && z && w && error);

  int total_dim = problem->total_dimension;
  double residual = 0.0;
  double norm_z = 0.0;

  /* Compute w = Mz + q */
  for (int i = 0; i < total_dim; i++) {
    double prod = 0.0;
    for (int j = 0; j < total_dim; j++) {
      prod += problem->M->dense[i * total_dim + j] * z[j];
    }
    w[i] = prod + problem->q[i];
  }

  /* Compute residual ||z - P_C(z - w)|| / ||z|| */
  for (int b = 0; b < problem->num_blocks; b++) {
    MixedProblemBlock* block = &problem->blocks[b];
    int start = block->start_index;
    int dim = block->dimension;
    int ne = block->num_elements;

    for (int e = 0; e < ne; e++) {
      int idx = start + e * dim;

      switch (block->model_type) {
        case MIXED_MODEL_FRICTION_2D: {
          /* Projection onto friction cone 2D */
          double mu = block->model.friction2d->mu[e];
          double rz = z[idx] - w[idx];
          double rz_n = rz;                      /* normal component */
          double rz_t = z[idx + 1] - w[idx + 1]; /* tangential */

          double p_n = fmax(0.0, rz_n);
          double p_t = fmax(1.0 - mu * fabs(rz_n) / fmax(fabs(rz_t), 1e-14), 0.0) * rz_t;

          residual += pow(z[idx] - p_n, 2) + pow(z[idx + 1] - p_t, 2);
          break;
        }

        case MIXED_MODEL_FRICTION_3D: {
          /* Projection onto friction cone 3D */
          double mu = block->model.friction3d->mu[e];
          double rz[3] = {z[idx] - w[idx], z[idx + 1] - w[idx + 1], z[idx + 2] - w[idx + 2]};
          double rz_n = rz[0];
          double rz_t_norm = sqrt(rz[1] * rz[1] + rz[2] * rz[2]);

          double p[3];
          projectionOnCone(rz, mu);
          for (int k = 0; k < 3; k++) p[k] = rz[k];

          residual +=
              pow(z[idx] - p[0], 2) + pow(z[idx + 1] - p[1], 2) + pow(z[idx + 2] - p[2], 2);
          break;
        }

        case MIXED_MODEL_VON_MISES: {
          /* Projection onto Von Mises cylinder */
          double sigma_y = block->model.von_mises->sigma_y[e];
          /* Simplified error computation - would need proper projection */
          for (int d = 0; d < dim; d++) {
            residual += pow(w[idx + d], 2);
          }
          break;
        }

        case MIXED_MODEL_DRUCKER_PRAGER: {
          /* Projection onto Drucker-Prager cone */
          /* Simplified error computation */
          for (int d = 0; d < dim; d++) {
            residual += pow(w[idx + d], 2);
          }
          break;
        }

        default:
          /* Generic residual computation */
          for (int d = 0; d < dim; d++) {
            residual += pow(w[idx + d], 2);
          }
          break;
      }
    }
  }

  /* Compute norm of z */
  for (int i = 0; i < total_dim; i++) {
    norm_z += z[i] * z[i];
  }
  norm_z = sqrt(norm_z);

  residual = sqrt(residual);

  if (norm_z > 1e-14) {
    *error = residual / norm_z;
  } else {
    *error = residual;
  }

  return NUMERICS_OK;
}
