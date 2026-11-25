// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdio.h>  // for printf, NULL

#include "NumericsFwd.h"    // for SolverOptions
#include "SolverOptions.h"  // for solver_options_delete, solver_options_create
#include "lcp_cst.h"        // for SICONOS_LCP_NSGS_SBM, SICONOS_LCP_AVI_CAO...

int main(void) {
  printf("\n Start of test on Default SolverOptions\n");
  int info = 0;

  int solvers[] = {SICONOS_LCP_NSGS_SBM,
                   SICONOS_LCP_LEMKE,
                   SICONOS_LCP_PGS,
                   SICONOS_LCP_CPG,
                   SICONOS_LCP_LATIN,
                   SICONOS_LCP_LATIN_W,
                   SICONOS_LCP_QP,
                   SICONOS_LCP_NSQP,
                   SICONOS_LCP_NEWTONMIN,
                   SICONOS_LCP_NEWTON_FB_FBLSA,
                   SICONOS_LCP_PSOR,
                   SICONOS_LCP_RPGS,
                   SICONOS_LCP_PATH,
                   SICONOS_LCP_ENUM,
                   SICONOS_LCP_AVI_CAOFERRIS,
                   SICONOS_LCP_PIVOT,
                   SICONOS_LCP_BARD,
                   SICONOS_LCP_MURTY,
                   SICONOS_LCP_NEWTON_MIN_FBLSA,
                   SICONOS_LCP_PATHSEARCH,
                   SICONOS_LCP_PIVOT_LUMOD,
                   SICONOS_LCP_GAMS,
                   SICONOS_LCP_CONVEXQP_PG};
  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));
  SolverOptions* options = NULL;
  for (int s = 0; s < n_solvers; ++s) {
    options = solver_options_create(solvers[s]);
    solver_options_print(options);
    solver_options_delete(options);
  }

  printf("\n End of test on Default SolverOptions\n");
  return info;
}
