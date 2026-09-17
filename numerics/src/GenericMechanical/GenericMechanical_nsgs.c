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

/*! \file GenericMechanical_nsgs.c
 *  \brief Non-Smooth Gauss-Seidel solver for GenericMechanicalProblem.
 */

#include <assert.h>  // for assert
#include <math.h>    // for fabs, sqrt, isnan
#ifndef __cplusplus
#include <stdbool.h>  // for true
#endif
#include <stdio.h>   // for printf, size_t
#include <stdlib.h>  // for free, malloc
#include <string.h>  // for NULL, memcpy

#include "FrictionContactProblem.h"     // IWYU pragma: keep
#include "GenericMechanicalProblem.h"   // for GMP_LocalProblem, GenericM...
#include "GenericMechanical_Solvers.h"  // for gmp_compute_error, gmp_gauss_seidel
#include "GenericMechanical_cst.h"      // for SICONOS_GENERIC_MECHANICAL...
#include "LCP_Solvers.h"                // for lcp_compute_error_only
#include "NonSmoothDrivers.h"           // for fc3d_driver, linearComplem...
#include "NumericsFwd.h"                // for SolverOptions
#include "NumericsMatrix.h"             // for NM_row_prod_no_diag, NM_ex...
#include "Relay_Solvers.h"              // for relay_compute_error
#include "SiconosBlas.h"                // for cblas_dnrm2, cblas_dgemv
#include "SolverOptions.h"              // for SolverOptions, solver_opti...
#include "fc2d_compute_error.h"         // for fc2d_unitary_compute_and_a...
#include "fc3d_compute_error.h"         // for fc3d_unitary_compute_and_a...
#include "numerics_verbose.h"

#include "numerics_errors.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"  // for DEBUG_PRINTF, DEBUG_EXPR

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

/* ===========================================================================
 * Local problem operations (vtable)
 * ===========================================================================
 * Each supported problem type registers four callbacks:
 *  - solve_local:     solve one block given its diagonal block and q_local
 *  - compute_error:   evaluate the local residual
 *  - detach_diag_block: clean up the view on the global diagonal block
 * ===========================================================================
 */

/** Operations for one local problem type inside a GenericMechanicalProblem. */
typedef struct {
  const char* name;
  int solver_options_index; /* index in options->internalSolvers used by this type */
  int (*solve_local)(GMP_LocalProblem* local, double* diag_block, double* sol, double* w,
                     SolverOptions* options);
  int (*compute_error)(GMP_LocalProblem* local, const double* reaction, const double* velocity,
                       double tol, double* err);
  void (*detach_diag_block)(GMP_LocalProblem* local);
} GMP_ProblemOps;

/* EQUALITY: linear system -M R = q_local */
static int gmp_equality_solve_local(GMP_LocalProblem* local, double* diag_block, double* sol,
                                    double* w, SolverOptions* options) {
  (void)w;
  (void)options;
  NumericsMatrix M;
  NM_fill(&M, NM_DENSE, local->size, local->size, diag_block);

  for (size_t i = 0; i < local->size; ++i) sol[i] = -local->q_local[i];

  int info = NM_LU_solve(NM_preserve(&M), sol, 1);
  NM_unpreserve(&M);
  M.matrix0 = NULL;
  NM_clear(&M);
  return info;
}

static int gmp_equality_compute_error(GMP_LocalProblem* local, const double* reaction,
                                      const double* velocity, double tol, double* err) {
  (void)local;
  (void)reaction;
  (void)tol;
  double localError = 0.;
  for (size_t i = 0; i < local->size; i++) {
    if (fabs(velocity[i]) > localError) localError = fabs(velocity[i]);
  }
  if (localError > *err) *err = localError;
  return 0;
}

static void gmp_equality_detach_diag_block(GMP_LocalProblem* local) { (void)local; }

/* LCP */
static int gmp_lcp_solve_local(GMP_LocalProblem* local, double* diag_block, double* sol,
                               double* w, SolverOptions* options) {
  LinearComplementarityProblem* lcpProblem = (LinearComplementarityProblem*)local->problem;
  lcpProblem->M->matrix0 = diag_block;
  return linearComplementarity_driver(lcpProblem, sol, w, options);
}

static int gmp_lcp_compute_error(GMP_LocalProblem* local, const double* reaction,
                                 const double* velocity, double tol, double* err) {
  (void)tol;
  double localError = 0.;
  lcp_compute_error_only(local->size, reaction, velocity, &localError);
  localError = localError / (1 + cblas_dnrm2(local->size, local->q_local, 1));
  if (localError > *err) *err = localError;
  return 0;
}

static void gmp_lcp_detach_diag_block(GMP_LocalProblem* local) {
  LinearComplementarityProblem* lcpProblem = (LinearComplementarityProblem*)local->problem;
  lcpProblem->M->matrix0 = NULL;
}

/* RELAY */
static int gmp_relay_solve_local(GMP_LocalProblem* local, double* diag_block, double* sol,
                                 double* w, SolverOptions* options) {
  RelayProblem* relayProblem = (RelayProblem*)local->problem;
  relayProblem->M->matrix0 = diag_block;
  return relay_driver(relayProblem, sol, w, options);
}

static int gmp_relay_compute_error(GMP_LocalProblem* local, const double* reaction,
                                   const double* velocity, double tol, double* err) {
  double localError = 0.;
  relay_compute_error((RelayProblem*)local->problem, reaction, velocity, tol, &localError);
  localError = localError / (1 + cblas_dnrm2(local->size, local->q_local, 1));
  if (localError > *err) *err = localError;
  return 0;
}

static void gmp_relay_detach_diag_block(GMP_LocalProblem* local) {
  RelayProblem* relayProblem = (RelayProblem*)local->problem;
  relayProblem->M->matrix0 = NULL;
}

/* FC3D */
static int gmp_fc3d_solve_local(GMP_LocalProblem* local, double* diag_block, double* sol,
                                double* w, SolverOptions* options) {
  FrictionContactProblem* fcProblem = (FrictionContactProblem*)local->problem;
  assert(fcProblem);
  assert(fcProblem->M);
  assert(fcProblem->q);
  fcProblem->M->matrix0 = diag_block;
  return fc3d_driver(fcProblem, sol, w, options);
}

static int gmp_fc3d_compute_error(GMP_LocalProblem* local, const double* reaction,
                                  const double* velocity, double tol, double* err) {
  (void)tol;
  FrictionContactProblem* fcProblem = (FrictionContactProblem*)local->problem;
  double localError = 0.;
  double worktmp[3];
  fc3d_unitary_compute_and_add_error(reaction, velocity, fcProblem->mu[0], &localError, worktmp);
  localError = sqrt(localError) / (1 + cblas_dnrm2(local->size, local->q_local, 1));
  if (localError > *err) *err = localError;
  return 0;
}

static void gmp_fc3d_detach_diag_block(GMP_LocalProblem* local) {
  FrictionContactProblem* fcProblem = (FrictionContactProblem*)local->problem;
  fcProblem->M->matrix0 = NULL;
}

/* FC2D */
static int gmp_fc2d_solve_local(GMP_LocalProblem* local, double* diag_block, double* sol,
                                double* w, SolverOptions* options) {
  FrictionContactProblem* fcProblem = (FrictionContactProblem*)local->problem;
  assert(fcProblem);
  assert(fcProblem->M);
  assert(fcProblem->q);
  fcProblem->M->matrix0 = diag_block;
  return fc2d_driver(fcProblem, sol, w, options);
}

static int gmp_fc2d_compute_error(GMP_LocalProblem* local, const double* reaction,
                                  const double* velocity, double tol, double* err) {
  (void)tol;
  FrictionContactProblem* fcProblem = (FrictionContactProblem*)local->problem;
  double localError = 0.;
  double worktmp[2];
  fc2d_unitary_compute_and_add_error(reaction, velocity, fcProblem->mu[0], &localError, worktmp);
  localError = sqrt(localError) / (1 + cblas_dnrm2(local->size, local->q_local, 1));
  if (localError > *err) *err = localError;
  return 0;
}

static void gmp_fc2d_detach_diag_block(GMP_LocalProblem* local) {
  FrictionContactProblem* fcProblem = (FrictionContactProblem*)local->problem;
  fcProblem->M->matrix0 = NULL;
}

/** Lookup table of operations, indexed by SICONOS_NUMERICS_PROBLEM_TYPE. */
static const GMP_ProblemOps* gmp_get_ops(int type) {
  static const GMP_ProblemOps equality_ops = {"EQUALITY", -1, gmp_equality_solve_local,
                                              gmp_equality_compute_error,
                                              gmp_equality_detach_diag_block};
  static const GMP_ProblemOps lcp_ops = {"LCP", 0, gmp_lcp_solve_local,
                                         gmp_lcp_compute_error, gmp_lcp_detach_diag_block};
  static const GMP_ProblemOps relay_ops = {"RELAY", 2, gmp_relay_solve_local,
                                           gmp_relay_compute_error, gmp_relay_detach_diag_block};
  static const GMP_ProblemOps fc2d_ops = {"FC2D", 3, gmp_fc2d_solve_local,
                                          gmp_fc2d_compute_error, gmp_fc2d_detach_diag_block};
  static const GMP_ProblemOps fc3d_ops = {"FC3D", 1, gmp_fc3d_solve_local,
                                          gmp_fc3d_compute_error, gmp_fc3d_detach_diag_block};

  switch (type) {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
      return &equality_ops;
    case SICONOS_NUMERICS_PROBLEM_LCP:
      return &lcp_ops;
    case SICONOS_NUMERICS_PROBLEM_RELAY:
      return &relay_ops;
    case SICONOS_NUMERICS_PROBLEM_FC2D:
      return &fc2d_ops;
    case SICONOS_NUMERICS_PROBLEM_FC3D:
      return &fc3d_ops;
    default:
      return NULL;
  }
}

/** Build the local right-hand side q for one block of a GenericMechanicalProblem.
 *
 *  For the global problem
 *      M reaction + q = velocity
 *  the local problem associated with block \a block_row is built from the
 *  off-diagonal contribution of the current iterate \a reaction:
 *
 *      q_local = q_global(block_row) + M_offdiag(block_row,:) * reaction
 *
 *  In practice, this function copies the corresponding block of the global
 *  vector problem->q into \a q_local and then adds the row-block product of M
 *  excluding the diagonal block (via NM_row_prod_no_diag).
 *
 *  \param[in]  problem         the global GenericMechanicalProblem
 *  \param[in]  block_row    block-row index in M (used for SBM storage)
 *  \param[in]  row_start    first row of the block in dense storage (unused for SBM)
 *  \param[in]  block_size   size of the local block
 *  \param[in]  reaction     current global reaction iterate
 *  \param[out] q_local      local right-hand side, must be allocated with size >= block_size
 */
static void gmp_build_local_q(const GenericMechanicalProblem* problem, int block_row,
                              size_t row_start, size_t block_size, const double* reaction,
                              double* q_local) {
  assert(problem);
  assert(problem->q);
  assert(reaction);
  assert(q_local);

  memcpy(q_local, &(problem->q[row_start]), block_size * sizeof(double));
  /* Add the off-diagonal block row product.
   * NM_row_prod_no_diag is not const-correct for the x argument but does not
   * modify reaction when called with init=0 and xsave=NULL. */
  NM_row_prod_no_diag(problem->globalSize, block_size, block_row, row_start, problem->M,
                      (double*)reaction, q_local, NULL, 0);
}

/* ===========================================================================
 * Gauss-Seidel workspace
 * ===========================================================================
 */

/** Workspace used by the Non-Smooth Gauss-Seidel iteration.
 *
 *  If \a owns_memory is false, the buffers are views into options->dWork and
 *  must not be freed by the workspace destructor.
 */
typedef struct {
  double* prev_reaction;
  double* buff_velocity;
  double* diag_dense_buffer;
  int owns_memory;
} GMP_GS_Workspace;

/** Create a workspace for gmp_gauss_seidel.
 *
 *  If \a options->dWork is non-NULL it is reused (prev_reaction/buff_velocity);
 *  otherwise memory is allocated and \a owns_memory is set to true.
 *  The dense diagonal buffer is allocated only when needed.
 */
static GMP_GS_Workspace gmp_gs_workspace_create(GenericMechanicalProblem* problem,
                                                SolverOptions* options,
                                                NM_types storageType) {
  GMP_GS_Workspace ws = {0};
  if (options->dWork) {
    ws.prev_reaction = options->dWork;
    ws.owns_memory = 0;
  } else {
    ws.prev_reaction = (double*)malloc(gmp_get_nb_dwork(problem, options) * sizeof(double));
    ws.owns_memory = 1;
  }
  ws.buff_velocity = ws.prev_reaction + problem->globalSize;

  if (storageType == NM_DENSE) {
    ws.diag_dense_buffer = (double*)malloc(problem->maxLocalSize * problem->maxLocalSize *
                                           sizeof(double));
  }
  return ws;
}

/** Destroy a Gauss-Seidel workspace.
 *
 *  Only frees buffers that were allocated by gmp_gs_workspace_create.
 */
static void gmp_gs_workspace_destroy(GMP_GS_Workspace* ws) {
  if (!ws) return;
  if (ws->diag_dense_buffer) {
    free(ws->diag_dense_buffer);
    ws->diag_dense_buffer = NULL;
  }
  if (ws->owns_memory && ws->prev_reaction) {
    free(ws->prev_reaction);
    ws->prev_reaction = NULL;
  }
  ws->buff_velocity = NULL;
}

int gmp_compute_error(const GenericMechanicalProblem* problem, double* reaction, double* velocity,
                      double tol, SolverOptions* options, double* err) {
  (void)options;
  GMP_LocalProblem* local = problem->firstLocal;
  NM_types storageType = problem->M->storageType;
  NumericsMatrix* numMat = problem->M;
  size_t block_row = 0;
  size_t local_size = 0;
  *err = 0.0;
  double* bufForLocalProblemDense =
      (storageType == NM_DENSE)
          ? (double*)malloc(problem->maxLocalSize * problem->maxLocalSize * sizeof(double))
          : 0;

  DEBUG_PRINT("GenericMechanical compute_error BEGIN\n");
  /* Update each local problem->q and compute V = M*R + Q of the GMP. */
  size_t global_offset = 0;
  while (local) {
    local_size = local->size;

    gmp_build_local_q(problem, block_row, global_offset, local_size, reaction, local->q_local);
    DEBUG_PRINT_MAT_STR("q", local->q_local, (unsigned)local_size, 1);
    DEBUG_PRINT_MAT_STR("reaction", reaction, (unsigned)problem->globalSize, 1);
    DEBUG_PRINT_MAT_STR("qnodiag", local->q_local, (unsigned)local_size, 1);
    /* Computation of the velocity of the GMP. */
    memcpy(velocity + global_offset, local->q_local, local_size * sizeof(double));
    /* Add the missing diagonal product to the velocity. */

    double* diagBlock = 0;
    if (storageType == NM_DENSE) {
      NM_extract_diag_block(numMat, block_row, global_offset, local_size,
                            &bufForLocalProblemDense);
      diagBlock = bufForLocalProblemDense;
    } else {
      NM_extract_diag_block(numMat, block_row, global_offset, local_size, &diagBlock);
    }
    DEBUG_PRINT_MAT_STR("diagBlock", diagBlock, (unsigned)local_size, (unsigned)local_size);
    DEBUG_PRINT_MAT_STR("Rlocal", reaction + global_offset, (unsigned)local_size, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, local_size, local_size, 1.0, diagBlock, local_size,
                reaction + global_offset, 1, 1.0, velocity + global_offset, 1);
    DEBUG_PRINT_MAT_STR("velocity", velocity + global_offset, (unsigned)local_size, 1);
    /* Next block. */
    global_offset += local->size;
    local = local->next;
    block_row++;
  }

  /* For each sub-problem, call the corresponding function computing the error. */
  global_offset = 0;
  block_row = 0;
  local = problem->firstLocal;
  while (local) {
    local_size = local->size;
    double* Vl = velocity + global_offset;
    double* Rl = reaction + global_offset;
    for (size_t ii = 0; ii < local_size; ii++)
      if (isnan(Vl[ii]) || isnan(Rl[ii])) {
        *err = 10;
        return 1;
      }

    const GMP_ProblemOps* ops = gmp_get_ops(local->type);
    if (ops) {
      ops->compute_error(local, Rl, Vl, tol, err);
      DEBUG_PRINTF("GenericMechanical_driver, localerror of %s: %e\n", ops->name, *err);
    } else {
      numerics_printf("Numerics : gmp_compute_error unknown problem type %d.\n",
                      local->type);
    }

    /* Next block. */
    global_offset += local->size;
    local = local->next;
    block_row++;
  }

  if (storageType == NM_DENSE) free(bufForLocalProblemDense);
  bufForLocalProblemDense = NULL;
  if (*err > tol)
    return 1;
  else
    return 0;
}

static void gmp_gauss_seidel_internal(GenericMechanicalProblem* problem, double* reaction,
                                      double* velocity, int* info, SolverOptions* options,
                                      GMP_GS_Workspace* ws) {
  DEBUG_BEGIN("gmp_gauss_seidel_internal(...)\n");

  GMP_LocalProblem* local = 0;
  NM_types storageType = problem->M->storageType;
  NumericsMatrix* numMat = problem->M;
  int iterMax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  int it = 0;
  size_t block_row = 0;
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  double* err = &(options->dparam[SICONOS_DPARAM_RESIDU]);
  double* errLS = &(options->dparam[SICONOS_DPARAM_GMP_ERROR_LS]);
  int tolViolate = 1;
  int tolViolateLS = 1;
  double* sol = 0;
  double* w = 0;
  int resLocalSolver = 0;
  int local_solver_error_occurred = 0;
  int withLS = options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_LINESEARCH];
  double* pCoefLS = &(options->dparam[SICONOS_DPARAM_GMP_COEFF_LS]);

  while (it < iterMax && tolViolate) {
    memcpy(ws->prev_reaction, reaction, problem->globalSize * sizeof(double));
    block_row = 0;
    local = problem->firstLocal;
    size_t global_offset = 0;
    size_t local_size = 0;

    DEBUG_PRINTF("GS it %d, initial value:\n", it);
    DEBUG_EXPR(for (size_t ii = 0; ii < problem->globalSize; ii++) numerics_printf(
                   "R[%i]=%e | V[%i]=%e ", ii, reaction[ii], ii, velocity[ii]););

    while (local) {
      numerics_printf_verbose(1, "Gauss-Seidel iteration %d. Problem (row) number %d ", it,
                              block_row);
      local_size = local->size;
      local->error = 0;

      /* Extract the diagonal block for the local solver. */
      double* diagBlock = 0;
      if (storageType == NM_DENSE) {
        NM_extract_diag_block(numMat, block_row, global_offset, local_size,
                              &ws->diag_dense_buffer);
        diagBlock = ws->diag_dense_buffer;
      } else {
        NM_extract_diag_block(numMat, block_row, global_offset, local_size, &diagBlock);
      }

      sol = reaction + global_offset;
      w = velocity + global_offset;

      /* Build local q and solve the local problem. */
      gmp_build_local_q(problem, block_row, global_offset, local_size, reaction, local->q_local);

      const GMP_ProblemOps* ops = gmp_get_ops(local->type);
      if (ops) {
        numerics_printf_verbose(1, "solve SICONOS_NUMERICS_PROBLEM_%s", ops->name);
        SolverOptions* local_opts =
            (ops->solver_options_index >= 0)
                ? options->internalSolvers[ops->solver_options_index]
                : NULL;
        resLocalSolver = ops->solve_local(local, diagBlock, sol, w, local_opts);
      } else {
        numerics_printf("genericMechanical_GS Numerics : gmp_gauss_seidel unknown problem type "
                        "%d.\n",
                        local->type);
        resLocalSolver = 1;
      }

      if (resLocalSolver) {
        local->error = 1;
        local_solver_error_occurred = 1;
      }

      global_offset += local->size;
      local = local->next;
      block_row++;
    }
    /* Compute global error. */

    if (withLS) {
      tolViolate = gmp_compute_error(problem, reaction, ws->buff_velocity, tol, options, err);
      for (size_t i = 0; i < problem->globalSize; i++)
        ws->prev_reaction[i] = reaction[i] + (*pCoefLS) * (reaction[i] - ws->prev_reaction[i]);
      tolViolateLS = gmp_compute_error(problem, ws->prev_reaction, velocity, tol, options, errLS);

      DEBUG_PRINTF("GMP :noscale error=%e error LS=%e\n", *err, *errLS);
      DEBUG_PRINTF("GMP :scale coeff=%e\n", *pCoefLS);

      if (*errLS < *err) {
        if ((*pCoefLS) < 10.0) (*pCoefLS) = 1.0 + (*pCoefLS);
        memcpy(reaction, ws->prev_reaction, problem->globalSize * sizeof(double));
        tolViolate = tolViolateLS;
        *err = *errLS;
      } else {
        *pCoefLS = 1.0;
        memcpy(velocity, ws->buff_velocity, problem->globalSize * sizeof(double));
      }
    } else {
      tolViolate = gmp_compute_error(problem, reaction, velocity, tol, options, err);
    }
    numerics_printf_verbose(
        1, "--------------- GMP - GS - Iteration %i Residual = %14.7e <= %7.3e\n", it, *err,
        options->dparam[SICONOS_DPARAM_TOL]);

    /* Next GS iteration. */
    it++;
  }

  /* Detach diagonal blocks that were only views on the global matrix. */
  local = problem->firstLocal;
  while (local) {
    const GMP_ProblemOps* ops = gmp_get_ops(local->type);
    if (ops) ops->detach_diag_block(local);
    local = local->next;
  }

  if (tolViolate) {
    if (verbose > 0)
      numerics_printf("gmp_gauss_seidel failed with Iteration %i Residual = %14.7e <= %7.3e\n",
                      it, *err, options->dparam[SICONOS_DPARAM_TOL]);
  }

  if (local_solver_error_occurred) {
    block_row = 0;
    local = problem->firstLocal;
    while (local) {
      if (local->error && verbose)
        numerics_printf(
            "genericMechanical_GS Numerics : Local solver FAILED row %d of type %s\n",
            block_row, ns_problem_id_to_name(local->type));
      local = local->next;
      block_row++;
    }
  }

  *info = tolViolate;
  if (local_solver_error_occurred && !*info) *info = 1;

  options->iparam[SICONOS_IPARAM_ITER_DONE] = it;

  DEBUG_END("gmp_gauss_seidel_internal(...)\n");
}

/** Print a summary of the local problems and their selected solvers.
 *
 *  Output is produced only when verbose > 0.
 */
static void gmp_print_local_problem_summary(GenericMechanicalProblem* problem,
                                            SolverOptions* options) {
  if (verbose <= 0) return;

  int count = 0;
  size_t total_size = 0;
  GMP_LocalProblem* local = problem->firstLocal;
  while (local) {
    count++;
    total_size += local->size;
    local = local->next;
  }

  numerics_printf("GenericMechanical NSGS: %d local problem(s), total size %zu", count,
                  total_size);

  local = problem->firstLocal;
  int block = 0;
  while (local) {
    const GMP_ProblemOps* ops = gmp_get_ops(local->type);
    const char* type_name = ns_problem_id_to_name(local->type);
    if (ops && ops->solver_options_index >= 0) {
      SolverOptions* local_opts = options->internalSolvers[ops->solver_options_index];
      numerics_printf("  block %d: %s size %zu (%s)", block, type_name, local->size,
                      solver_options_id_to_name(local_opts->solverId));
    } else {
      numerics_printf("  block %d: %s size %zu (direct solve)", block, type_name, local->size);
    }
    local = local->next;
    block++;
  }
}

void gmp_gauss_seidel(GenericMechanicalProblem* problem, double* reaction, double* velocity,
                      int* info, SolverOptions* options) {
  DEBUG_BEGIN("gmp_gauss_seidel(...)\n");
  DEBUG_EXPR_WE({
    FILE* toto1 = fopen("GMP_CURRENT.txt", "w");
    if (toto1) {
      genericMechanicalProblem_printInFile(problem, toto1);
      fclose(toto1);
    }
  });

  gmp_print_local_problem_summary(problem, options);

  NM_types storageType = problem->M->storageType;
  GMP_GS_Workspace ws = gmp_gs_workspace_create(problem, options, storageType);
  gmp_gauss_seidel_internal(problem, reaction, velocity, info, options, &ws);
  gmp_gs_workspace_destroy(&ws);

  DEBUG_END("gmp_gauss_seidel(...)\n");
}

/** Allocate working memory if not already present.
 *
 *  This function allocates options->dWork only if it is currently NULL.
 *  The layout expected by the GS workspace is:
 *    [0 .. globalSize-1]              prev_reaction
 *    [globalSize .. 2*globalSize-1]   buff_velocity
 */
int gmp_working_memory_alloc(GenericMechanicalProblem* problem, SolverOptions* options) {
  if (options->dWork) return 0;

  options->dWork = (double*)malloc(gmp_get_nb_dwork(problem, options) * sizeof(double));
  return 1;
}

/** Return the size (number of doubles) needed for options->dWork. */
int gmp_get_nb_dwork(GenericMechanicalProblem* problem, SolverOptions* options) {
  (void)options;
  return 2 * problem->globalSize;
}
