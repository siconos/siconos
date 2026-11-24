/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdio.h>   // for NULL
#include <stdlib.h>  // for malloc

#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for LinearComplementarityProblem
#include "NumericsMatrix.h"                // for NumericsMatrix, RawNumeric...
#include "fc2d_Solvers.h"                  // for fc2d_tolcp
#include "numerics_errors.h"
#include "numerics_verbose.h"  // for numerics_error

int fc2d_tolcp(FrictionContactProblem* problem, LinearComplementarityProblem* lcp_problem) {
  if (problem->dimension != 2) {
    return numerics_error(
        "fc2d_tolcp",
        "Dimension of the problem : problem-> dimension is not compatible or is not set");
  }
  int nc = problem->numberOfContacts;
  lcp_problem->size = 3 * nc;
  lcp_problem->M = NM_new();
  lcp_problem->M->size0 = 3 * nc;
  lcp_problem->M->size1 = 3 * nc;

  lcp_problem->M->storageType = 0;
  lcp_problem->M->sparse_block = NULL;
  lcp_problem->M->sparse = NULL;
  lcp_problem->M->internalData = NULL;
  lcp_problem->M->dense =
      (double*)malloc(lcp_problem->size * lcp_problem->size * sizeof(double));
  ;
  lcp_problem->q = (double*)malloc(lcp_problem->size * sizeof(double));
  int n = 2 * nc;
  int i, j;
  for (i = 0; i < nc; i++) {
    for (j = 0; j < nc; j++) {
      /* first Column */
      lcp_problem->M->dense[3 * i + 3 * j * lcp_problem->size] =
          problem->M->dense[2 * i + 2 * j * n] -
          problem->mu[i] * problem->M->dense[2 * i + (2 * j + 1) * n];  // compute W_NN-mu W_NT

      lcp_problem->M->dense[3 * i + 1 + 3 * j * lcp_problem->size] =
          problem->M->dense[(2 * i + 1) + 2 * j * n] -
          problem->mu[i] *
              problem->M->dense[(2 * i + 1) + (2 * j + 1) * n];  // compute W_TN-mu W_TT

      if (i == j) {
        lcp_problem->M->dense[3 * i + 2 + 3 * j * lcp_problem->size] = 2.0 * problem->mu[i];
      } else {
        lcp_problem->M->dense[3 * i + 2 + 3 * j * lcp_problem->size] = 0.0;
      }

      /* second Column */
      lcp_problem->M->dense[3 * i + (3 * j + 1) * lcp_problem->size] =
          problem->M->dense[2 * i + (2 * j + 1) * n];  // set W_NT

      lcp_problem->M->dense[3 * i + 1 + (3 * j + 1) * lcp_problem->size] =
          problem->M->dense[(2 * i + 1) + (2 * j + 1) * n];  // set WTT

      if (i == j) {
        lcp_problem->M->dense[3 * i + 2 + (3 * j + 1) * lcp_problem->size] = -1.0;
      } else {
        lcp_problem->M->dense[3 * i + 2 + (3 * j + 1) * lcp_problem->size] = 0.0;
      }

      /* Third Column */
      lcp_problem->M->dense[3 * i + (3 * j + 2) * lcp_problem->size] = 0.0;

      if (i == j) {
        lcp_problem->M->dense[3 * i + 1 + (3 * j + 2) * lcp_problem->size] = 1.0;
      } else {
        lcp_problem->M->dense[3 * i + 1 + (3 * j + 2) * lcp_problem->size] = 0.0;
      }

      lcp_problem->M->dense[3 * i + 2 + (3 * j + 2) * lcp_problem->size] = 0.0;
    }
  }

  for (i = 0; i < nc; i++) {
    lcp_problem->q[3 * i] = problem->q[2 * i];
    lcp_problem->q[3 * i + 1] = problem->q[2 * i + 1];
    lcp_problem->q[3 * i + 2] = 0.0;
    ;
  }

  return 0;
}
