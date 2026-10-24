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
#include "plasticity_2d_local_problem_tools.h"

#include "NM_types.h"
#ifndef __cplusplus
#include <stdbool.h>  // for false
#endif
#include <stdlib.h>  // for malloc, NULL

#include "NumericsMatrix.h"         // for NM_create_with_data, NumericsMatrix
#include "PlasticityProblem.h"      // for PlasticityProblem,
#include "SparseBlockMatrix.h"      // IWYU pragma: keep
#include "local_problem_toolkit.h"  // for local_problem_toolkit_new
#include "numerics_errors.h"
#include "numerics_verbose.h"

struct LocalPLASTICITY_2DProblemFunctionToolkit* localPLASTICITY_2DProblemFunctionToolkit_new(
    void) {
  return (struct LocalPLASTICITY_2DProblemFunctionToolkit*)local_problem_toolkit_new(
      sizeof(struct LocalPLASTICITY_2DProblemFunctionToolkit));
}

void localPLASTICITY_2DProblemFunctionToolkit_display(
    struct LocalPLASTICITY_2DProblemFunctionToolkit* lpft) {
  printf("local_solver %p\n ", lpft->local_solver);
  printf("update_local_problem %p\n ", lpft->update_local_problem);
  printf("post_processed_local_result %p\n ", lpft->post_processed_local_result);
  printf("free_local_solver %p\n ", lpft->free_local_solver);
};
void plasticity_2d_local_problem_compute_q(PlasticityProblem* problem,
                                           PlasticityProblem* localproblem, double* reaction,
                                           int contact) {
  double* qLocal = localproblem->q;

  int n = 3 * problem->numberOfCones;

  int in = 3 * contact, it = in + 1, is = it + 1;

  /* qLocal computation*/
  qLocal[0] = problem->q[in];
  qLocal[1] = problem->q[it];
  qLocal[2] = problem->q[is];

  CHECK_RETURN(NM_row_prod_no_diag3(n, contact, 3 * contact, problem->M, reaction, qLocal,
                                    false) == NUMERICS_OK);
}

void plasticity_2d_local_problem_fill_M(PlasticityProblem* problem,
                                        PlasticityProblem* localproblem, int contact) {
  if (problem->M->storageType == NM_SPARSE) {
    localproblem->M->dense = problem->M->sparse_block->block[contact];
  } else
    CHECK_RETURN(NM_extract_diag_block3(problem->M, contact, &localproblem->M->dense) ==
                 NUMERICS_OK);
}

PlasticityProblem* plasticity_2d_local_problem_allocate(PlasticityProblem* problem) {
  /* Connect local solver and local problem*/
  PlasticityProblem* localproblem = (PlasticityProblem*)malloc(sizeof(PlasticityProblem));
  localproblem->numberOfCones = 1;
  localproblem->dimension = 3;
  localproblem->q = (double*)calloc(3, sizeof(double));

  /* Copy model type from parent problem */
  localproblem->model_type = problem->model_type;

  /* Allocate model-specific data based on model type */
  switch (problem->model_type) {
    case PLASTICITY_MODEL_DRUCKER_PRAGER: {
      localproblem->model.drucker_prager =
          (Plasticity_DruckerPrager_model*)malloc(sizeof(Plasticity_DruckerPrager_model));
      localproblem->model.drucker_prager->theta = (double*)malloc(sizeof(double));
      localproblem->model.drucker_prager->eta = (double*)malloc(sizeof(double));
      break;
    }
    case PLASTICITY_MODEL_VON_MISES: {
      localproblem->model.von_mises =
          (Plasticity_VonMises_model*)malloc(sizeof(Plasticity_VonMises_model));
      localproblem->model.von_mises->sigma_y = (double*)malloc(sizeof(double));
      break;
    }
    default:
      localproblem->model.generic = NULL;
      break;
  }

  if (problem->M->storageType == NM_DENSE) {
    localproblem->M = NM_create_with_data(NM_DENSE, 3, 3, calloc(9, sizeof(double)));
  } else /* NM_SPARSE_BLOCK or NM_SPARSE*/
  {
    localproblem->M = NM_create(NM_DENSE, 3, 3);
    // We need to create "localproblem->M" but not matrix0, since it will be set with a link to
    // problem->M blocks
  }
  return localproblem;
}

void plasticity_2d_local_problem_free(PlasticityProblem* localproblem,
                                      PlasticityProblem* problem) {
  if (problem->M->storageType == NM_SPARSE || problem->M->storageType == NM_SPARSE_BLOCK) {
    // localproblem->M has been set using plasticity_2d_local_problem_fill_M, with
    // a pointer link when problem->M is NM_SPARSE.
    // So we need to release the pointer to avoid double free deallocation of the diagonal
    // blocks of the original matrix of the problem.
    localproblem->M->dense = NULL;
  }
  plasticity2DProblem_free(localproblem);
}
