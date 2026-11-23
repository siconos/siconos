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
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"

/* Include friction contact solvers */
#include "FrictionContactProblem.h"
#include "fc2d_Solvers.h"
#include "fc2d_compute_error.h"
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
#include "plasticity_2d_druckerprager_projection.h"
#include "plasticity_2d_solvers.h"
#include "plasticity_2d_vonmises.h"
#include "rolling_friction_2d_compute_error.h"
#include "rolling_friction_3d_compute_error.h"

/* Local solver data for each block */
typedef struct MixedProblemLocalSolverData {
  void* local_problem;
  SolverOptions* local_options;
  int owns_local_options; /* local_options created here (rolling blocks), freed with the data
                           */
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
    case MIXED_MODEL_ROLLING_FRICTION_2D: {
      /* Create local 2D rolling friction problem */
      RollingFrictionContactProblem* rfc2d = rollingFrictionContactProblem_new();
      if (!rfc2d) return NUMERICS_ERR_MEMORY;

      rfc2d->dimension = 3; /* (r_n, r_t, r_rolling) */
      rfc2d->numberOfContacts = ne;
      rfc2d->mu = block->model.rolling_friction2d->mu;
      rfc2d->mu_r = block->model.rolling_friction2d->mur;

      rfc2d->M = NM_create_with_data(NM_DENSE, local_size, local_size,
                                     calloc((size_t)local_size * local_size, sizeof(double)));
      if (!rfc2d->M) {
        free(rfc2d);
        return NUMERICS_ERR_MEMORY;
      }

      rfc2d->q = solver_data->local_q;
      solver_data->local_problem = rfc2d;

      if (!solver_data->local_options) {
        solver_data->local_options = solver_options_create(SICONOS_ROLLING_FRICTION_2D_NSGS);
        if (!solver_data->local_options) return NUMERICS_ERR_MEMORY;
        solver_data->owns_local_options = 1;
      }

      break;
    }

    case MIXED_MODEL_ROLLING_FRICTION_3D: {
      /* Create local 3D rolling friction problem */
      RollingFrictionContactProblem* rfc3d = rollingFrictionContactProblem_new();
      if (!rfc3d) return NUMERICS_ERR_MEMORY;

      rfc3d->dimension = 5; /* (r_n, r_t1, r_t2, r_roll1, r_roll2) */
      rfc3d->numberOfContacts = ne;
      rfc3d->mu = block->model.rolling_friction3d->mu;
      rfc3d->mu_r = block->model.rolling_friction3d->mur;

      rfc3d->M = NM_create_with_data(NM_DENSE, local_size, local_size,
                                     calloc((size_t)local_size * local_size, sizeof(double)));
      if (!rfc3d->M) {
        free(rfc3d);
        return NUMERICS_ERR_MEMORY;
      }

      rfc3d->q = solver_data->local_q;
      solver_data->local_problem = rfc3d;

      if (!solver_data->local_options) {
        solver_data->local_options = solver_options_create(SICONOS_ROLLING_FRICTION_3D_NSGS);
        if (!solver_data->local_options) return NUMERICS_ERR_MEMORY;
        solver_data->owns_local_options = 1;
      }

      break;
    }

    case MIXED_MODEL_FRICTION_2D:
    case MIXED_MODEL_FRICTION_3D:
    case MIXED_MODEL_VON_MISES:
    case MIXED_MODEL_DRUCKER_PRAGER:
      /* element-wise Gauss-Seidel on the global problem: no local problem */
      solver_data->local_problem = NULL;
      break;

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

  /* local_options belong to the internalSolvers of the main SolverOptions, unless they were
   * created by mixedProblem_local_solver_init */
  if (solver_data->owns_local_options) solver_options_delete(solver_data->local_options);
  solver_data->local_options = NULL;

  if (solver_data->local_problem) {
    /* Before freeing the problem, set shared pointers to NULL to avoid
     * double-free. These pointers point to memory owned by MixedProblem
     * or already freed above (local_q). */
    switch (model_type) {
      case MIXED_MODEL_ROLLING_FRICTION_2D: {
        RollingFrictionContactProblem* rfc2d =
            (RollingFrictionContactProblem*)solver_data->local_problem;
        rfc2d->mu = NULL;   /* Don't free */
        rfc2d->mu_r = NULL; /* Don't free */
        rfc2d->q = NULL;    /* Don't free - already freed above */
        rollingFrictionContactProblem_free(rfc2d);
        break;
      }
      case MIXED_MODEL_ROLLING_FRICTION_3D: {
        RollingFrictionContactProblem* rfc3d =
            (RollingFrictionContactProblem*)solver_data->local_problem;
        rfc3d->mu = NULL;   /* Don't free */
        rfc3d->mu_r = NULL; /* Don't free */
        rfc3d->q = NULL;    /* Don't free - already freed above */
        rollingFrictionContactProblem_free(rfc3d);
        break;
      }
      default:
        break;
    }
  }
}

/* Dense local matrix (diagonal block of M) stored in the local problem of a block. */
static NumericsMatrix* mixedProblem_local_M(MixedModelType model_type,
                                            MixedProblemLocalSolverData* solver_data) {
  switch (model_type) {
    case MIXED_MODEL_FRICTION_2D:
    case MIXED_MODEL_FRICTION_3D:
      return ((FrictionContactProblem*)solver_data->local_problem)->M;
    case MIXED_MODEL_ROLLING_FRICTION_2D:
    case MIXED_MODEL_ROLLING_FRICTION_3D:
      return ((RollingFrictionContactProblem*)solver_data->local_problem)->M;
    case MIXED_MODEL_VON_MISES:
    case MIXED_MODEL_DRUCKER_PRAGER:
      return ((PlasticityProblem*)solver_data->local_problem)->M;
    default:
      return NULL;
  }
}

/* Copy the diagonal block of the global M into the local M (dense, column-major).
 * M does not change during the iterations: done once, before the NSGS loop. */
static void mixedProblem_fill_local_M(MixedProblem* problem, int block_index,
                                      MixedProblemLocalSolverData* solver_data) {
  MixedProblemBlock* block = &problem->blocks[block_index];
  int start = block->start_index;
  int local_size = block->dimension * block->num_elements;
  NumericsMatrix* local_M = mixedProblem_local_M(block->model_type, solver_data);
  if (!local_M) return;
  /* NM_get_value works for all the storages of the global M */
  for (int j = 0; j < local_size; j++)
    for (int i = 0; i < local_size; i++)
      local_M->dense[i + j * local_size] = NM_get_value(problem->M, start + i, start + j);
}

/* Copy the block of the global z into the local z. */
static void mixedProblem_extract_local_block(MixedProblem* problem, int block_index, double* z,
                                             MixedProblemLocalSolverData* solver_data) {
  MixedProblemBlock* block = &problem->blocks[block_index];
  int start = block->start_index;
  int local_size = block->dimension * block->num_elements;
  for (int i = 0; i < local_size; i++) solver_data->local_z[i] = z[start + i];
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

/* Solve the one-element problem (dimension dim) of element e of a block:
 *   w_e = Me z_e + qe,  Me (dim x dim, dense, column-major) the diagonal block of element e.
 * z_e is updated in place. */
static int mixedProblem_solve_element(MixedProblemBlock* block, int e, double* Me, double* qe,
                                      double* ze, SolverOptions* options) {
  int dim = block->dimension;
  NumericsMatrix Me_nm;
  NM_null(&Me_nm);
  Me_nm.storageType = NM_DENSE;
  Me_nm.size0 = dim;
  Me_nm.size1 = dim;
  Me_nm.dense = Me;

  switch (block->model_type) {
    case MIXED_MODEL_FRICTION_2D: {
      /* exact solution of the 2D one-contact problem (enumeration of the cases) */
      double D = Me[0] * Me[3] - Me[2] * Me[1];
      fc2d_nsgs_local_solve(Me, D, qe, block->model.friction2d->mu[e], ze);
      return 0;
    }
    case MIXED_MODEL_FRICTION_3D: {
      FrictionContactProblem contact = {3, 1, &Me_nm, qe, &block->model.friction3d->mu[e]};
      return fc3d_projectionOnCone_solve(&contact, ze, options);
    }
    case MIXED_MODEL_DRUCKER_PRAGER: {
      Plasticity_DruckerPrager_model dp = {&block->model.drucker_prager->eta[e],
                                           &block->model.drucker_prager->theta[e], NULL};
      /* cohesion c: the element problem is solved for sigma~ = sigma + a e_N, a = c / eta,
       * with q~ = qe - M_ee a e_N (cone without cohesion) */
      const double* c = block->model.drucker_prager->cohesion;
      double a = (c && c[e] != 0.0) ? c[e] / block->model.drucker_prager->eta[e] : 0.0;
      double qt[3] = {qe[0] - Me[0] * a, qe[1] - Me[1] * a, qe[2] - Me[2] * a};
      PlasticityProblem cone;
      cone.dimension = 3;
      cone.numberOfCones = 1;
      cone.M = &Me_nm;
      cone.q = qt;
      cone.model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
      cone.model.drucker_prager = &dp;
      ze[0] += a;
      int info = plasticity_2d_projectionOnCone_solve(&cone, ze, options);
      ze[0] -= a;
      return info;
    }
    case MIXED_MODEL_VON_MISES:
      /* exact solution of the one-element Von Mises (J2) problem */
      return (plasticity_2d_vonMises_return_mapping(Me, qe, block->model.von_mises->sigma_y[e],
                                                    ze) == 1);
    default:
      return NUMERICS_ERR_INVALID_ARGUMENT;
  }
}

/* Diagonal block (dim x dim) of element `elem` (global element number, i.e. block row number
 * in a sparse-block M) starting at row `row` of M. The specialised
 * NM_extract_diag_block{2,3,5} assume row = dim * elem, which holds when all the previous
 * elements have the same dimension; otherwise the generic NM_extract_diag_block is used. For a
 * sparse-block M the returned pointer is the block of the SBM (no copy), otherwise `buffer` is
 * filled. */
static double* mixedProblem_extract_diag_block(NumericsMatrix* M, int elem, size_t row,
                                               int dim, double* buffer) {
  double* B = buffer;
  int same_size = (row == (size_t)dim * (size_t)elem);
  if (same_size && dim == 2)
    NM_extract_diag_block2(M, elem, &B);
  else if (same_size && dim == 3)
    NM_extract_diag_block3(M, elem, &B);
  else if (same_size && dim == 5)
    NM_extract_diag_block5(M, elem, &B);
  else
    NM_extract_diag_block(M, elem, row, dim, &B);
  return B;
}

/* qe = q(row:row+dim) + (row of M) z without the diagonal block of element `elem`. */
static void mixedProblem_element_q(MixedProblem* problem, int elem, size_t row, int dim,
                                   double* z, double* qe) {
  size_t n = (size_t)problem->total_dimension;
  for (int i = 0; i < dim; i++) qe[i] = problem->q[row + i];
  switch (dim) {
    case 1:
      NM_row_prod_no_diag1x1(n, elem, row, problem->M, z, qe, false);
      break;
    case 2:
      NM_row_prod_no_diag2(n, elem, row, problem->M, z, qe, false);
      break;
    case 3:
      NM_row_prod_no_diag3(n, elem, row, problem->M, z, qe, false);
      break;
    default:
      NM_row_prod_no_diag(n, (size_t)dim, elem, row, problem->M, z, qe, NULL, false);
      break;
  }
}

/* One Gauss-Seidel sweep over the elements (contacts, cones) of a block, directly on the
 * global problem w = M z + q (any storage of M). For element e (rows r..r+dim-1):
 *   M_ee z_e + q_e + sum_{f != e} M_ef z_f
 * with the z_f already updated (previous elements of this block and previous blocks), solved
 * in place by mixedProblem_solve_element:
 *  - Friction 2D: fc2d_nsgs_local_solve (exact);
 *  - Friction 3D: fc3d_projectionOnCone_solve (one projection step);
 *  - Drucker-Prager (Mohr-Coulomb): plasticity_2d_projectionOnCone_solve (one projection step
 *    on K_eta with the dilatancy term theta);
 *  - Von Mises: plasticity_2d_vonMises_return_mapping (exact).
 * first_elem is the global number of the first element of the block (block row number of a
 * sparse-block M, which must have one diagonal block per element).
 */
static int mixedProblem_block_sweep(MixedProblem* problem, int block_index, int first_elem,
                                    double* z, SolverOptions* options) {
  enum { DIM_MAX = 5 };
  MixedProblemBlock* block = &problem->blocks[block_index];
  int dim = block->dimension;
  assert(dim <= DIM_MAX);
  double Me_buffer[DIM_MAX * DIM_MAX], qe[DIM_MAX];

  int info = 0;
  for (int e = 0; e < block->num_elements; e++) {
    size_t row = (size_t)(block->start_index + dim * e);
    int elem = first_elem + e;
    double* Me = mixedProblem_extract_diag_block(problem->M, elem, row, dim, Me_buffer);
    mixedProblem_element_q(problem, elem, row, dim, z, qe);
    info |= mixedProblem_solve_element(block, e, Me, qe, z + row, options);
  }
  return info;
}

/* Rolling-friction blocks: the block problem (local M = diagonal block of M, local q with the
 * coupling to the other blocks) is solved by the rolling-friction NSGS solver.
 * With a sparse-block M, NM_row_prod_no_diag* remove a single SBM block: the rolling block
 * must then have one contact. */
static int mixedProblem_solve_rolling_block(MixedProblem* problem, int block_index,
                                            int first_elem,
                                            MixedProblemLocalSolverData* solver_data,
                                            double* z_global) {
  MixedProblemBlock* block = &problem->blocks[block_index];
  int info = 0;
  int start = block->start_index;
  int local_size = block->dimension * block->num_elements;

  /* coupling = (rows of the block of M) z without the diagonal block (local_w as storage) */
  double* coupling = solver_data->local_w;

  switch (local_size) {
    case 1:
      NM_row_prod_no_diag1x1((size_t)problem->total_dimension, first_elem, (size_t)start,
                             problem->M, z_global, coupling, true);
      break;
    case 2:
      NM_row_prod_no_diag2((size_t)problem->total_dimension, first_elem, (size_t)start,
                           problem->M, z_global, coupling, true);
      break;
    case 3:
      NM_row_prod_no_diag3((size_t)problem->total_dimension, first_elem, (size_t)start,
                           problem->M, z_global, coupling, true);
      break;
    default:
      NM_row_prod_no_diag((size_t)problem->total_dimension, (size_t)local_size, first_elem,
                          (size_t)start, problem->M, z_global, coupling, NULL, true);
      break;
  }

  /* local_q = global_q + coupling: w_block = M_bb z_b + (q_b + sum_{c != b} M_bc z_c) */
  for (int i = 0; i < local_size; i++) {
    solver_data->local_q[i] = problem->q[start + i] + coupling[i];
  }

  switch (block->model_type) {
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

  int max_iter = options->iparam[SICONOS_MIXED_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_MIXED_DPARAM_TOL];
  int error_evaluation = options->iparam[SICONOS_MIXED_IPARAM_ERROR_EVALUATION];
  int frequency = options->iparam[SICONOS_MIXED_IPARAM_ERROR_EVALUATION_FREQUENCY];
  if (frequency <= 0) frequency = 1;
  double error = INFINITY;
  int iter = 0;
  *info = 1;

  /* Initial guess z = 0 */
  for (int i = 0; i < total_dim; i++) z[i] = 0.0;

  /* Diagonal blocks of the rolling-friction blocks (solved block-wise) */
  for (int b = 0; b < num_blocks; b++) {
    MixedModelType model = problem->blocks[b].model_type;
    if (model == MIXED_MODEL_ROLLING_FRICTION_2D || model == MIXED_MODEL_ROLLING_FRICTION_3D)
      mixedProblem_fill_local_M(problem, b, &nsgs_data->block_solvers[b]);
  }

  while (iter < max_iter) {
    memcpy(nsgs_data->z_prev, z, total_dim * sizeof(double));

    /* Gauss-Seidel sweep over the blocks */
    int first_elem = 0; /* global number of the first element of the block */
    for (int b = 0; b < num_blocks; b++) {
      MixedProblemBlock* block = &problem->blocks[b];
      MixedProblemLocalSolverData* solver_data = &nsgs_data->block_solvers[b];
      int local_info;
      if (block->model_type == MIXED_MODEL_ROLLING_FRICTION_2D ||
          block->model_type == MIXED_MODEL_ROLLING_FRICTION_3D) {
        mixedProblem_extract_local_block(problem, b, z, solver_data);
        local_info = mixedProblem_solve_rolling_block(problem, b, first_elem, solver_data, z);
        mixedProblem_update_global_block(problem, b, z, w, solver_data);
      } else {
        local_info = mixedProblem_block_sweep(problem, b, first_elem, z, options);
      }
      if (local_info != 0 && verbose > 0)
        printf("MixedProblem_NSGS: local solver failed for block %d\n", b);
      first_elem += block->num_elements;
    }
    iter++;

    /* Light error: relative change of z */
    double diff2 = 0.0, norm_z2 = 0.0;
    for (int i = 0; i < total_dim; i++) {
      double d = z[i] - nsgs_data->z_prev[i];
      diff2 += d * d;
      norm_z2 += z[i] * z[i];
    }
    double light_error = sqrt(diff2);
    if (norm_z2 > DBL_EPSILON * DBL_EPSILON) light_error /= sqrt(norm_z2);

    int converged = 0;
    switch (error_evaluation) {
      case SICONOS_MIXED_ERROR_EVALUATION_LIGHT:
        error = light_error;
        converged = (error < tolerance);
        break;
      case SICONOS_MIXED_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL:
        error = light_error;
        if (light_error < tolerance)
          converged = !mixedProblem_compute_error(problem, z, w, tolerance, options, &error);
        break;
      case SICONOS_MIXED_ERROR_EVALUATION_FULL:
      default:
        if (iter % frequency == 0)
          converged = !mixedProblem_compute_error(problem, z, w, tolerance, options, &error);
        break;
    }
    if (verbose > 0)
      printf("MixedProblem_NSGS: iter %d, light error = %e, error = %e\n", iter, light_error,
             error);
    if (converged) {
      *info = 0;
      break;
    }
  }

  /* w = M z + q (and the full error, unless only the light error is asked for) */
  if (error_evaluation == SICONOS_MIXED_ERROR_EVALUATION_LIGHT) {
    cblas_dcopy(total_dim, problem->q, 1, w, 1);
    NM_gemv(1.0, problem->M, z, 1.0, w);
  } else {
    mixedProblem_compute_error(problem, z, w, tolerance, options, &error);
  }

  options->dparam[SICONOS_MIXED_DPARAM_RESIDU] = error;
  options->iparam[SICONOS_MIXED_IPARAM_ITER_DONE] = iter;
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
        (options->internalSolvers && (size_t)b < options->numberOfInternalSolvers)
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
  if (!problem || !options) return;

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

  /* Zero-initialize all fields first to ensure safe cleanup */
  memset(options, 0, sizeof(SolverOptions));

  options->solverId = SICONOS_MIXED_NSGS;
  options->numberOfInternalSolvers = 0;
  options->internalSolvers = NULL;
  options->isSet = true;

  options->iparam = (int*)calloc(OPTIONS_PARAM_SIZE, sizeof(int));
  options->dparam = (double*)calloc(OPTIONS_PARAM_SIZE, sizeof(double));

  options->iparam[SICONOS_MIXED_IPARAM_MAX_ITER] = 1000;
  options->iparam[SICONOS_MIXED_IPARAM_ITER_DONE] = 0;
  /* stop on the relative change of z, confirmed by the full (natural-map) error */
  options->iparam[SICONOS_MIXED_IPARAM_ERROR_EVALUATION] =
      SICONOS_MIXED_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_MIXED_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1;
  options->iparam[SICONOS_MIXED_IPARAM_SHUFFLE] = 0;

  options->dparam[SICONOS_MIXED_DPARAM_TOL] = 1e-6;
  options->dparam[SICONOS_MIXED_DPARAM_RESIDU] = 0.0;

  /* Note: We use direct local solvers that don't need internal solver options.
   * For FC2D: fc2d_nsgs_local_solve (direct, no options needed)
   * For FC3D: fc3d_projectionOnCone_solve (direct, no options needed)
   * For plasticity: direct projection functions (no options needed)
   */
}

int mixedProblem_compute_error(MixedProblem* problem, double* z, double* w, double tolerance,
                               SolverOptions* options, double* error) {
  assert(problem && z && w && error);
  (void)options;
  int n = problem->total_dimension;

  /* w = M z + q */
  cblas_dcopy(n, problem->q, 1, w, 1);
  NM_gemv(1.0, problem->M, z, 1.0, w);

  /* Sum over the elements of the squared natural-map residuals, with the same convention as
   * the single-model solvers (fc2d/fc3d/rolling/plasticity compute_error):
   *   z_e - P_K(z_e - (w_e + De Saxce/dilatancy term)) */
  double err = 0.0, worktmp[5];
  for (int b = 0; b < problem->num_blocks; b++) {
    MixedProblemBlock* block = &problem->blocks[b];
    for (int e = 0; e < block->num_elements; e++) {
      int idx = block->start_index + e * block->dimension;
      double* ze = z + idx;
      double* we = w + idx;
      switch (block->model_type) {
        case MIXED_MODEL_FRICTION_2D:
          fc2d_unitary_compute_and_add_error(ze, we, block->model.friction2d->mu[e], &err,
                                             worktmp);
          break;
        case MIXED_MODEL_FRICTION_3D:
          fc3d_unitary_compute_and_add_error(ze, we, block->model.friction3d->mu[e], &err,
                                             worktmp);
          break;
        case MIXED_MODEL_ROLLING_FRICTION_2D:
          rolling_friction_2d_unitary_compute_and_add_error(
              ze, we, block->model.rolling_friction2d->mu[e],
              block->model.rolling_friction2d->mur[e], &err, worktmp);
          break;
        case MIXED_MODEL_ROLLING_FRICTION_3D:
          rolling_friction_3d_unitary_compute_and_add_error(
              ze, we, block->model.rolling_friction3d->mu[e],
              block->model.rolling_friction3d->mur[e], &err, worktmp);
          break;
        case MIXED_MODEL_VON_MISES:
          /* z - P_C(z - w), C the Von Mises cylinder */
          plasticity_2d_vonmises_unitary_compute_and_add_error(
              ze, we, block->model.von_mises->sigma_y[e], &err, worktmp);
          break;
        case MIXED_MODEL_DRUCKER_PRAGER: {
          /* cohesion c: error for sigma + c/eta e_N (cone without cohesion) */
          const double* c = block->model.drucker_prager->cohesion;
          double eta = block->model.drucker_prager->eta[e];
          double r[3] = {ze[0], ze[1], ze[2]};
          if (c && c[e] != 0.0) r[0] += c[e] / eta;
          plasticity_2d_unitary_compute_and_add_error(
              r, we, eta, block->model.drucker_prager->theta[e], &err, worktmp);
          break;
        }
        default:
          return numerics_error("mixedProblem_compute_error", "unknown model type %d",
                                block->model_type);
      }
    }
  }
  err = sqrt(err);

  /* relative error (|q| scaling, as in fc3d_compute_error) */
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  if (norm_q > DBL_EPSILON) err /= norm_q;
  *error = err;
  return (err > tolerance);
}
