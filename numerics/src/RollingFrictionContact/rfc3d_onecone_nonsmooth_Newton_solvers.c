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
//#include <H5Apublic.h>
#include "rfc3d_onecone_nonsmooth_Newton_solvers.h"  // for computeNonsmoo...

#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs, isinf
#include <stdio.h>   // for NULL, printf
#include <stdlib.h>  // for calloc, realloc

#include "FrictionContact_options.h"
#include "NSSTools.h"                          // for max
#include "NonSmoothSolvers/NonSmoothNewton.h"  // for nonSmoothDirec...
#include "NumericsFwd.h"                       // for SolverOptions
#include "NumericsMatrix.h"                    // for NumericsMatrix
#include "RollingFrictionContactProblem.h"
#include "SiconosBlas.h"    // for cblas_ddot
#include "SolverOptions.h"  // for SolverOptions
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "op3x3.h"  // for cpy3, mvp3x3
#include "op5x5.h"  // for cpy3, mvp3x3
#include "rfc3d_short_names.h"
#include "rolling_fc2d_projection.h"
#include "rolling_fc3d_compute_error.h"
#include "rolling_fc3d_local_problem_tools.h"
#include "rolling_fc3d_projection.h"
#include "rolling_naturalmap_functions.h"

/* #define DEBUG_CHECK */
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"  // for DEBUG_PRINTF

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif
static computeNonsmoothFunction Function = NULL;
static NewtonFunctionPtr F = NULL;
static NewtonFunctionPtr jacobianF = NULL;

/* size of a block */
static int Fsize;

static double rfc3d_compute_local_error(RollingFrictionContactProblem* localproblem,
                                        double* local_reaction) {
  double* local_q = localproblem->q;
  double norm_q = cblas_dnrm2(5, localproblem->q, 1);
  double* local_M = localproblem->M->matrix0;
  double worktmp[5];
  double local_velocity[5] = {0., 0., 0., 0., 0.};
  cpy5(local_q, local_velocity);
  mvp5x5(local_M, local_reaction, local_velocity);
  DEBUG_EXPR(NV_display(local_q, 5););
  DEBUG_EXPR(NV_display(local_velocity, 5););
  DEBUG_EXPR(NV_display(local_reaction, 5););

  double current_error = 0.0;

  rolling_fc3d_unitary_compute_and_add_error(local_reaction, local_velocity,
                                             localproblem->mu[0], localproblem->mu_r[0],
                                             &current_error, worktmp);
  current_error = sqrt(current_error);
  DEBUG_PRINTF("absolute local error = %e", current_error);
  if (fabs(norm_q) > DBL_EPSILON) current_error /= norm_q;
  DEBUG_PRINTF("relative local error = %e", current_error);
  return current_error;
}

static void compute_rho_spectral_norm(RollingFrictionContactProblem* localproblem,
                                      double* rho) {
  double* MLocal = localproblem->M->matrix0;
  double eigvals[5];
  int info = eigvals_sym5x5_for_loop(MLocal, eigvals);

  DEBUG_PRINTF("compute_rho_spectral_norm : eigvals\n");
  // display5(eigvals);
  rho[0] = 1.0 / eigvals[4];
}

static void rfc3d_onecone_nonsmooth_Newton_initialize(
    RollingFrictionContactProblem* problem, RollingFrictionContactProblem* localproblem,
    SolverOptions* options) {
  /** In initialize, these operators are "connected" to their corresponding static variables,
   * that will be used to build local problem for each considered cone.
   * Local problem is built during call to update (which depends on the storage type for M).
   */

  DEBUG_PRINTF(
      "rfc3d_onecone_nonsmooth_Newton_initialize starts with "
      "options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = "
      "%i\n",
      options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]);

  if (options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] ==
      SICONOS_FRICTION_3D_NSN_FORMULATION_NATURALMAP) {
    numerics_printf(
        " SICONOS_FRICTION_3D_NSN_FORMULATION =   "
        "SICONOS_FRICTION_3D_NSN_FORMULATION_NATURALMAP");
    Function = &(rolling_friction_3D_computeNaturalMap);
  } else if (options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] ==
             SICONOS_FRICTION_3D_NSN_FORMULATION_NULL) {
    Function = NULL;
  }

  /* Compute and store default value of rho value */
  size_t nc = problem->numberOfContacts;

  if (options->solverId == SICONOS_ONECONE_NSN ||
      options->solverId == SICONOS_ONECONE_NSN_GP) {
    if (!options->dWork || options->dWorkSize < nc) {
      options->dWork = (double*)realloc(options->dWork, nc * sizeof(double));
      options->dWorkSize = nc;
    }
  } else if (options->solverId == SICONOS_ONECONE_NSN_GP_HYBRID) {
    if (!options->dWork || options->dWorkSize < 2 * nc) {
      options->dWork = (double*)realloc(options->dWork, 2 * nc * sizeof(double));
      options->dWorkSize = 2 * nc;
    }
  }

  double* rho = 0;
  for (size_t cone = 0; cone < nc; cone++) {
    if (options->solverId == SICONOS_ONECONE_NSN ||
        options->solverId == SICONOS_ONECONE_NSN_GP) {
      rho = &options->dWork[cone];
    } else if (options->solverId == SICONOS_ONECONE_NSN_GP_HYBRID) {
      options->dWork[cone] = 1.0;  // for PLI algorithm.
      rho = &options->dWork[cone + nc];
    }
    numerics_printf_verbose(2,
                            "rfc3d_onecone_nonsmooth_Newton_initialize"
                            " rfc3d_compute rho for cone = %i",
                            cone);

    switch (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY]) {
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM: {
        rolling_fc3d_local_problem_fill_M(problem, localproblem, cone);
        compute_rho_spectral_norm(localproblem, rho);
        break;
      }
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT: {
        rho[0] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
        break;
      }
      case SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE: {
        numerics_warning("rfc3d_onecone_nonsmooth_Newton_initialize",
                         "Adaptive strategy for computing rho not yet implemented");
      }
      default: {
        numerics_printf_verbose(2, "rfc3d_onecone_nonsmooth_Newton_initialize",
                                "strategy for computing rho is unknown");
        numerics_printf_verbose(2, "rfc3d_onecone_nonsmooth_Newton_initialize",
                                "switch back to RHO_STRATEGY_SPECTRAL_NORM");
        rolling_fc3d_local_problem_fill_M(problem, localproblem, cone);
        compute_rho_spectral_norm(localproblem, rho);
      }
    }
    double avg_rho = 0.0;

    if (verbose > 0) {
      avg_rho += rho[0];

      double m_row_norm = 0.0, sum;
      for (int i = 0; i < 5; i++) {
        sum = 0.0;
        for (int j = 0; j < 5; j++) {
          sum += fabs(localproblem->M->matrix0[i + j * 5]);
        }
        m_row_norm = max(sum, m_row_norm);
      }
      // to be updated
      numerics_printf_verbose(2,
                              "rfc3d_onecone_nonsmooth_Newton_initialize : "
                              "cone = %i, rho[0] = %4.2e\t,"
                              "||M||^-1 = %4.2e, ||M||_row^-1 = %4.2e ",
                              cone, rho[0], 1.0 / hypot9(localproblem->M->matrix0),
                              1.0 / m_row_norm);
    }
    numerics_printf(
        "rfc3d_onecone_nonsmooth_Newton_initialize"
        " Avg. rho value = %e\t",
        avg_rho / nc);
    DEBUG_EXPR(NM_display(localproblem->M););
  }
}

static void rfc3d_free(RollingFrictionContactProblem* problem,
                       RollingFrictionContactProblem* localproblem,
                       SolverOptions* localsolver_options) {
  F = NULL;
  jacobianF = NULL;
  free(localsolver_options->dWork);
  localsolver_options->dWork = NULL;
}

void rfc3d_onecone_nonsmooth_Newton_solvers_initialize(
    RollingFrictionContactProblem* problem, RollingFrictionContactProblem* localproblem,
    SolverOptions* localsolver_options) {
  /* Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen
   * formulation. */

  if (localsolver_options->solverId == SICONOS_ONECONE_NSN) {
    rfc3d_onecone_nonsmooth_Newton_initialize(problem, localproblem, localsolver_options);
  } else {
    numerics_error("rfc3d_onecone_nonsmooth_Newton_solvers_initialize",
                   "Unknown formulation type.");
  }
}

int rfc3d_onecone_nonsmooth_Newton_solvers_solve(RollingFrictionContactProblem* localproblem,
                                                 double* r_local, SolverOptions* options) {
  /* numerics_printf_verbose( */
  /*     2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve starts for cone %i",
   */
  /*     options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]); */

  int info = 1;

  /*  check trivial solution */

  double* q = localproblem->q;
  if (q[0] > -DBL_EPSILON) {
    r_local[0] = 0.0;
    r_local[1] = 0.0;
    r_local[2] = 0.0;
    r_local[3] = 0.0;
    r_local[4] = 0.0;
    numerics_printf_verbose(
        2, "rfc3d_nsn  take off, trivial solution reaction = 0, velocity = q.");
    info = 0;
    SET_SOLVER_ITER_DONE(options, 0);
    SET_SOLVER_RESIDUAL(options, 0.0);
    /* numerics_printf_verbose( */
    /*     2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve ends"); */
    return info;
  }
  // to be updated
  if (options->solverId == SICONOS_ONECONE_NSN) {
    info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct(localproblem, r_local, options);
  } /* else if (options->solverId == SICONOS_ONECONE_NSN_GP) { */
  /*   info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped(localproblem, local_reaction,
   */
  /*                                                              options); */
  /* } else if (options->solverId == SICONOS_ONECONE_NSN_GP_HYBRID) { */
  /*   if (options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] == */
  /*           SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP || */
  /*       options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] == */
  /*           SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP) { */
  /*     info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid(localproblem,
     local_reaction, */
  /*                                                                options); */
  /*   } */
  /* else { */
  /*   numerics_error("rfc3d_onecone_nonsmooth_Newton_solvers_solve", */
  /* 		   "Unknown local nsn hybrid solver"); */
  /* } */
  /* } */
  else {
    info = nonSmoothDirectNewton(Fsize, r_local, &F, &jacobianF, options);
  }
  if (info > 0) {
    if (verbose > 0) {
      if (SOLVER_MAX_ITER(options) == SOLVER_ITER_DONE(options)) {
        numerics_warning(
            "rfc3d_onecone_nonsmooth_Newton_solvers_solve",
            "reached max. number of iterations (%i) for cone %i. Residual = %12.8e",
            SOLVER_MAX_ITER(options), SOLVER_CURRENT_BLOCK(options), SOLVER_RESIDUAL(options));
      } else {
        numerics_warning("rfc3d_onecone_nonsmooth_Newton_solvers_solve",
                         "no convergence for cone %i with error = %12.8e",
                         SOLVER_MAX_ITER(options), SOLVER_CURRENT_BLOCK(options),
                         SOLVER_RESIDUAL(options));
      }
      /* note : exit on failure should be done in DefaultCheckSolverOutput */
    }
  }
  /* numerics_printf_verbose(2, */
  /*                         "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve
   * ends"); */
  return info;
  /*  (*postSolver)(cone,reaction); */
}

void rfc3d_onecone_nonsmooth_Newton_solvers_free(RollingFrictionContactProblem* problem,
                                                 RollingFrictionContactProblem* localproblem,
                                                 SolverOptions* localsolver_options) {
  F = NULL;
  jacobianF = NULL;
  if (localsolver_options->solverId == SICONOS_ONECONE_NSN ||
      localsolver_options->solverId == SICONOS_ONECONE_NSN_GP ||
      localsolver_options->solverId == SICONOS_ONECONE_NSN_GP_HYBRID) {
    rfc3d_free(problem, localproblem, localsolver_options);
  } else {
    numerics_error("rfc3d_onecone_nonsmooth_Newton_solvers_initialize",
                   "Unknown formulation type.");
  }
}

void rfc3d_onecone_nonsmooth_Newton_solvers_computeError(int n, double* velocity,
                                                         double* reaction, double* error) {}

void rfc3d_onecone_nonsmooth_Newton_update(int cone, RollingFrictionContactProblem* problem,
                                           RollingFrictionContactProblem* localproblem,
                                           double* reaction, SolverOptions* options) {
  /* Build a local problem for a specific cone
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific cone
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  rolling_fc3d_local_problem_fill_M(problem, localproblem, cone);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.reactionBlock, excluding the block corresponding to the current cone. ****/
  rolling_fc3d_local_problem_compute_q(problem, localproblem, reaction, cone);

  /*  coefficient for current block*/
  localproblem->mu[0] = problem->mu[cone];
  localproblem->mu_r[0] = problem->mu_r[cone];
}

int rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct(
    RollingFrictionContactProblem* localproblem, double* r_local, SolverOptions* options) {
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* numerics_printf_verbose( */
  /*     2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct starts"); */

  double mu = localproblem->mu[0];
  double mu_r = localproblem->mu_r[0];
  double* qLocal = localproblem->q;

  double norm_qLocal =
      sqrt(qLocal[0] * qLocal[0] + qLocal[1] * qLocal[1] + qLocal[2] * qLocal[2] +
           qLocal[3] * qLocal[3] + qLocal[4] * qLocal[4]);
  double norm_relative = 1.0;
  if (norm_qLocal > DBL_EPSILON) {
    norm_relative /= norm_qLocal;
  }

  double* MLocal = localproblem->M->matrix0;

  /* store the increment */
  double dR[5] = {0., 0., 0., 0., 0.};
  /* store the value fo the function */
  double F[5] = {0., 0., 0., 0., 0.};
  /* Store the (sub)-gradient of the function */
  double A[25];
  zero5x5(A);
  double B[25];
  zero5x5(B);
  /* Value of AW+B */
  double AWplusB[25];
  zero5x5(AWplusB);

  /* retrieve value of rho */
  double* rho;
  if (options->solverId == SICONOS_ONECONE_NSN ||
      options->solverId == SICONOS_ONECONE_NSN_GP) {
    rho = &options->dWork[iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]];
  } else {
    int nc = options->dWorkSize / 4;
    rho = &options->dWork[iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] + nc];
  }
  /* compute the velocity */
  double v_local[5] = {0., 0., 0., 0., 0.};
  cpy5(qLocal, v_local);
  mvp5x5(MLocal, r_local, v_local);

  int itermax = SOLVER_MAX_ITER(options);

  /* Newton iteration */
  int inew = 0;
  int info_solv5x5;

  /* compute first residue */
  Function(r_local, v_local, mu, mu_r, rho, F, NULL, NULL);
  /* SET_SOLVER_RESIDUAL(options, */
  /*     0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]+ F[3] * F[3] + F[4] * F[4]) *
   * norm_relative); */
  SET_SOLVER_RESIDUAL(
      options, sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2] + F[3] * F[3] + F[4] * F[4]) *
                   norm_relative);

  // double error = rfc3d_compute_local_error(localproblem, R);
  // printf("New local error = %e\n", error);
  if (iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] == 0)
    numerics_printf_verbose(2, "%-12s %-6s %-4s %-12s %-12s %-12s %-12s", "solver", "cone",
                            "it", "|dR|", "residu", "error", "tol");
  numerics_printf_verbose(2, "%-12s %-6d %-4d %-12.4e %-12.4e %-12.4e %-12.4e", "rfc3d_nsn",
                          iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew,
                          sqrt(dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2] * dR[2] +
                               dR[3] * dR[3] + dR[4] * dR[4]),
                          SOLVER_RESIDUAL(options),
                          rfc3d_compute_local_error(localproblem, r_local),
                          SOLVER_TOL(options));
  if (SOLVER_RESIDUAL(options) < SOLVER_TOL(options)) {
    SET_SOLVER_ITER_DONE(options, inew);
    return 0;
  }

  for (inew = 1; inew < itermax; ++inew) {
    /* Update function and gradient */
    Function(r_local, v_local, mu, mu_r, rho, F, A, B);
    DEBUG_EXPR(printf("F"); display5(F));
    /* compute -(A MLocal +B) */
    mm5x5(A, MLocal, AWplusB);
    add5x5(B, AWplusB);
    scal5x5(-1., AWplusB);
    DEBUG_EXPR(printf("AWplusB"); display5x5(AWplusB));

    /* Solve the linear system */
    cpy5(F, dR);
    solve_5x5_gepp_for_loop(AWplusB, dR);

    DEBUG_EXPR(printf("dR"); display5(dR));
    DEBUG_EXPR(printf("F"); display5(F));
    /* update iterates */
    add5(dR, r_local);

    /* compute new residue */
    cpy5(qLocal, v_local);
    mvp5x5(MLocal, r_local, v_local);
    Function(r_local, v_local, mu, mu_r, rho, F, NULL, NULL);

    //    dparam[SICONOS_DPARAM_RESIDU] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]+ F[3]
    //    * F[3] + F[4] * F[4]) *
    //                               norm_relative;  // improve with relative tolerance

    // VA 01/2025 : should be better but prevents convergence due to rounding error
    SET_SOLVER_RESIDUAL(
        options,
        sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2] * F[2] + F[3] * F[3] + F[4] * F[4]) *
            norm_relative);  // improve with relative tolerance

    // double error = rfc3d_compute_local_error(localproblem, R);
    // printf("New local error = %e\n", error);

    /* if (info_solv3x3) { */
    /*   if (verbose > 0) */
    /*     numerics_warning("rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct", */
    /*                      "cone %i 3x3 linear system is irregular # iteration = %" */
    /*                      "merit  = %.10e", */
    /*                      iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew, */
    /*                      dparam[SICONOS_DPARAM_RESIDU]); */
    /*   break; */
    /* } */

    numerics_printf_verbose(2, "%-12s %-6d %-4d %-12.4e %-12.4e %-12.4e %-12.4e", "rfc3d_nsn",
                            iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew,
                            sqrt(dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2] * dR[2] +
                                 dR[3] * dR[3] + dR[4] * dR[4]),
                            SOLVER_RESIDUAL(options),
                            rfc3d_compute_local_error(localproblem, r_local),
                            SOLVER_TOL(options));

    if (SOLVER_RESIDUAL(options) < SOLVER_TOL(options)) {
      SET_SOLVER_ITER_DONE(options, inew);
      return 0;
    }
    if (sqrt(dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2] * dR[2] + dR[3] * dR[3] +
             dR[4] * dR[4]) < 10 * DBL_EPSILON) {
      return 1;
    }
  }  // End of the Newton iteration

  SET_SOLVER_ITER_DONE(options, inew);
  return 1;
}

#undef DEBUG_MESSAGES
static int LineSearchGP(RollingFrictionContactProblem* localproblem,
                        computeNonsmoothFunction Function, double* t_opt, double R[3],
                        double dR[3], double* rho, int LSitermax, double* F, double* A,
                        double* B, double* velocity) {
  // to be implemented
  assert(0);
  int info = 0;
  return -1;
}
#define DEBUG_MESSAGES
int rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped(
    RollingFrictionContactProblem* localproblem, double* R, SolverOptions* options) {
  // to be implemented
  assert(0);
  int info = 0;
  return 1;
}

static void keep_or_discard_solution(RollingFrictionContactProblem* localproblem,
                                     double* local_reaction, double* local_reaction_backup,
                                     SolverOptions* options, double* current_error) {
  // to be implemented
  assert(0);
}

int rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid(
    RollingFrictionContactProblem* localproblem, double* local_reaction,
    SolverOptions* options) {
  // to be implemented
  assert(0);
  int info = 0;
  return info;
}

void rfc3d_onecone_nsn_set_default(SolverOptions* options) {
  /* Value of rho parameter */
  options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND;

  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.0;

  /* Choice of formulation */
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_NATURALMAP;

  /* Choice of line -search method */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] = SICONOS_FRICTION_3D_NSN_LINESEARCH_NO;

  /* parameters for hybrid solvers */
  options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] =
      SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP;

  options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP] = 1;
  options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 10;
}

void rfc3d_onecone_nsn_gp_set_default(SolverOptions* options) {
  /* Value of rho parameter */
  options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM;

  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.0;

  /* Choice of formulation */
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_NATURALMAP;
  /* Choice of line -search method */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] =
      SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 10;

  /* parameters for hybrid solvers */
  options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] =
      SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP;
  options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP] = 1;
  options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 100;
}
