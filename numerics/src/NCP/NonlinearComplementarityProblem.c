/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "NonlinearComplementarityProblem.h"

#include <stdlib.h>  // for free, NULL, malloc

#include "NumericsMatrix.h"  // for NM_clear
#include "numerics_errors.h"

void NCP_free(NonlinearComplementarityProblem* ncp) {
  if (ncp->nabla_F) {
    NM_clear(ncp->nabla_F);
    free(ncp->nabla_F);
    ncp->nabla_F = NULL;
  }

  free(ncp);
}

NonlinearComplementarityProblem* newNCP(void) {
  NonlinearComplementarityProblem* ncp =
      (NonlinearComplementarityProblem*)malloc(sizeof(NonlinearComplementarityProblem));

  ncp->n = 0;
  ncp->compute_F = NULL;
  ncp->compute_nabla_F = NULL;
  ncp->nabla_F = NULL;
  ncp->env = NULL;

  return ncp;
}
