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
#include "fc3d_compute_rho.h"

#include <assert.h>  // for assert
#include <math.h>    // for sqrt
#include <stddef.h>  // for size_t

#include "FrictionContactProblem.h"    // for FrictionContactProblem
#include "NumericsMatrix.h"            // for NM_extract_diagonal_blocks
#include "NumericsVector.h"            // for NM_dense_display
#include "SparseBlockMatrix.h"         // for SBM_free, SBM_FREE_ALL
#include "fc3d_local_problem_tools.h"  // for fc3d_local_problem_allocate, ...
#include "numerics_verbose.h"          // for numerics_error, numerics_printf_verbose
#include "op3x3.h"                     // for eig_3x3
#include "siconos_debug.h"             // for DEBUG_PRINTF, DEBUG_EXPR

void compute_rho_split_spectral_norm_cond(FrictionContactProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->dense;
  assert(MLocal[0 + 0 * 3] > 0);

  DEBUG_EXPR(NM_dense_display(MLocal, 3, 3, 3););
  double sw = MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3];

  double dw = sw * sw - 4.0 * (sw - MLocal[2 + 1 * 3] + MLocal[1 + 2 * 3]);
  DEBUG_PRINTF("dw = %e\n", dw);
  if (dw > 0.0)
    dw = sqrt(dw);
  else
    dw = 0.0;

  rho[0] = 1.0 / MLocal[0 + 0 * 3];
  rho[1] = 2.0 * (sw - dw) / ((sw + dw) * (sw + dw));
  rho[2] = rho[1];

  assert(rho[0] > 0);
  assert(rho[1] > 0);
  assert(rho[2] > 0);

  DEBUG_PRINTF("sw=%le\t  ", sw);
  DEBUG_PRINTF("dw=%le\n ", dw);
  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\n", rho[2]);
}

void compute_rho_split_spectral_norm(FrictionContactProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->dense;
  assert(MLocal[0 + 0 * 3] > 0);

  DEBUG_EXPR(NM_dense_display(MLocal, 3, 3, 3););
  double sw = MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3];

  double dw = sw * sw - 4.0 * (sw - MLocal[2 + 1 * 3] + MLocal[1 + 2 * 3]);
  DEBUG_PRINTF("dw = %e\n", dw);
  if (dw > 0.0)
    dw = sqrt(dw);
  else
    dw = 0.0;

  rho[0] = 1.0 / MLocal[0 + 0 * 3];

  rho[1] = 2.0 / (sw + dw);
  rho[2] = rho[1];

  assert(rho[0] > 0);
  assert(rho[1] > 0);
  assert(rho[2] > 0);

  DEBUG_PRINTF("sw=%le\t  ", sw);
  DEBUG_PRINTF("dw=%le\n ", dw);
  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\n", rho[2]);
}

void compute_rho_spectral_norm(FrictionContactProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->dense;
  double worktmp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eig[3] = {0.0, 0.0, 0.0};
  if (eig_3x3(MLocal, worktmp, eig)) numerics_printf("compute_rho_spectral_norm : failed");
  DEBUG_PRINTF("eig[0] = %4.2e, eig[1] = %4.2e, eig[2] = %4.2e", eig[0], eig[1], eig[2]);
  DEBUG_PRINTF("1/eig[0] = %4.2e, 1/eig[1] = %4.2e, 1/eig[2] = %4.2e", 1.0 / eig[0],
               1.0 / eig[1], 1.0 / eig[2]);
  rho[0] = 1.0 / eig[0];
  rho[1] = rho[0];
  rho[2] = rho[0];

  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\n", rho[2]);
}

int fc3d_compute_rho(FrictionContactProblem* problem, int rho_strategy, double rho_value,
                     double* rho) {
  assert(problem);
  assert(problem->dimension == 3);

  int nc = problem->numberOfContacts;
  if (nc == 0) return 0;

  FrictionContactProblem* localproblem = fc3d_local_problem_allocate(problem->M->storageType);

  SparseBlockStructuredMatrix* matrix1 = NULL;
  if (problem->M->storageType == NM_SPARSE) {
    if (problem->M->sparse_block) {
      numerics_printf("fc3d_compute_rho: Warning matrix 1 different from NULL");
    }
    matrix1 = problem->M->sparse_block;
    problem->M->sparse_block =
        NM_extract_diagonal_blocks(problem->M, (size_t)problem->dimension);
  }

  for (int contact = 0; contact < nc; ++contact) {
    double* rho_i = &rho[3 * contact];

    switch (rho_strategy) {
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND: {
        fc3d_local_problem_fill_M(problem, localproblem, contact);
        compute_rho_split_spectral_norm_cond(localproblem, rho_i);
        break;
      }
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM: {
        fc3d_local_problem_fill_M(problem, localproblem, contact);
        compute_rho_split_spectral_norm(localproblem, rho_i);
        break;
      }
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM: {
        fc3d_local_problem_fill_M(problem, localproblem, contact);
        compute_rho_spectral_norm(localproblem, rho_i);
        break;
      }
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT: {
        rho_i[0] = rho_value;
        rho_i[1] = rho_value;
        rho_i[2] = rho_value;
        break;
      }
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE: {
        return numerics_error("fc3d_compute_rho",
                              "Adaptive strategy for computing rho not yet implemented");
      }
      default: {
        return numerics_error("fc3d_compute_rho", "unknown strategy for computing rho");
      }
    }

    numerics_printf_verbose(2,
                            "fc3d_compute_rho"
                            "contact = %i, rho[0] = %4.2e, rho[1] = %4.2e, rho[2] = %4.2e",
                            contact, rho_i[0], rho_i[1], rho_i[2]);
  }

  if (problem->M->storageType == NM_SPARSE) {
    SBM_free(problem->M->sparse_block, SBM_FREE_ALL);
    problem->M->sparse_block = matrix1;
    localproblem->M->dense = NULL;
  }
  if (problem->M->storageType == NM_SPARSE_BLOCK) {
    /* we release the pointer to avoid deallocation of the diagonal blocks of the original
     * matrix of the problem*/
    localproblem->M->dense = NULL;
  }
  fc3d_local_problem_free(localproblem, problem);

  return 0;
}
