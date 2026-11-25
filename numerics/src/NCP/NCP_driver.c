/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <assert.h>  // for assert
#include <stdio.h>   // for fprintf, stderr
#include <stdlib.h>  // for exit, EXIT_FAILURE

#include "NCP_Solvers.h"                      // for ncp_compute_error, ncp_...
#include "NCP_cst.h"                          // for SICONOS_NCP_NEWTON_FB_F...
#include "NonSmoothDrivers.h"                 // for ncp_driver
#include "NonlinearComplementarityProblem.h"  // for NonlinearComplementarit...
#include "NumericsFwd.h"                      // for SolverOptions, Nonlinea...
#include "SolverOptions.h"                    // for SolverOptions
#include "numerics_errors.h"
#include "sn_error_handling.h"  // for sn_fatal_error_msg, SN_...
int ncp_driver(NonlinearComplementarityProblem* problem, double* z, double* F,
               SolverOptions* options) {
  /* Input validation using standardized macros */
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(F);
  CHECK_OPTIONS(options);

  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  int info_jmp = SN_SETJMP_INTERNAL_START;
  if (info_jmp == SN_NO_ERROR) {
    switch (options->solverId) {
      case SICONOS_NCP_NEWTON_FB_FBLSA:  // Fischer-Burmeister + Newton w/ LS
        ncp_newton_FBLSA(problem, z, F, &info, options);
        break;
      case SICONOS_NCP_NEWTON_MIN_FBLSA:  // min (+ FB as backup) + Newton w/ LS
        ncp_newton_minFBLSA(problem, z, F, &info, options);
        break;
      case SICONOS_NCP_PATHSEARCH:  // pathsearch method
        ncp_pathsearch(problem, z, F, &info, options);
        break;
      case SICONOS_NCP_PATH:  // PATH method
        ncp_path(problem, z, F, &info, options);
        break;
      default:
        CHECK_ARG(0, "ncp_driver error: unknown solver id: %d\n");
    }

    /* check the conditions 0 <= z _|_ F(z) >= 0 */
    if (options->filterOn > 0) {
      int info_ =
          ncp_compute_error(problem->n, z, F, options->dparam[0], &(options->dparam[1]));
      if (info <= 0) /* info was not set or the solver was happy */
        info = info_;
    }
    SN_SETJMP_INTERNAL_STOP
  } else {
    fprintf(stderr, "Fatal error in ncp_driver: %s", sn_fatal_error_msg());
    info = info_jmp;
  }

  return info;
}
