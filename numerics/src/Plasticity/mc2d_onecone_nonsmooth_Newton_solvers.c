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
#include "mc2d_onecone_nonsmooth_Newton_solvers.h"  // for computeNonsmoo...

#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs, isinf
#include <stdio.h>   // for NULL, printf
#include <stdlib.h>  // for calloc, realloc

#include "MohrCoulomb2DProblem.h"
#include "NSSTools.h"         // for max
#include "NonSmoothNewton.h"  // for nonSmoothDirec...
#include "NumericsFwd.h"      // for SolverOptions
#include "NumericsMatrix.h"   // for NumericsMatrix
#include "Plasticity_cst.h"
#include "SiconosBlas.h"                  // for cblas_ddot
#include "SolverOptions.h"                // for SolverOptions
#include "mc2d_AlartCurnier_functions.h"  // for compute_rho_sp...
#include "mc2d_compute_error.h"           // for mc2d_unitary_c...
#include "mc2d_local_problem_tools.h"     // for mc2d_local_pro...
#include "mc2d_naturalmap_functions.h"    // for mc2d_computeNaturalMap
#include "mc2d_projection.h"              // for mc2d_projectio...
#include "mc2d_solvers.h"                 // for FreeSolverNSGSPtr
#include "numerics_verbose.h"             // for numerics_print...
#include "op3x3.h"                        // for cpy3, mvp3x3

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

static double mc2d_compute_local_error(MohrCoulomb2DProblem* localproblem,
                                       double* local_reaction) {
  double* local_q = localproblem->q;
  double norm_q = cblas_dnrm2(3, localproblem->q, 1);
  double* local_M = localproblem->M->matrix0;
  double worktmp[3];
  double local_velocity[3] = {0., 0., 0.};
  cpy3(local_q, local_velocity);
  mvp3x3(local_M, local_reaction, local_velocity);
  DEBUG_EXPR(NV_display(local_q, 3););
  DEBUG_EXPR(NV_display(local_velocity, 3););
  DEBUG_EXPR(NV_display(local_reaction, 3););

  double current_error = 0.0;

  mc2d_unitary_compute_and_add_error(local_reaction, local_velocity, localproblem->eta[0],
                                     localproblem->theta[0], &current_error, worktmp);
  current_error = sqrt(current_error);
  DEBUG_PRINTF("absolute local error = %e", current_error);
  if (fabs(norm_q) > DBL_EPSILON) current_error /= norm_q;
  DEBUG_PRINTF("relative local error = %e", current_error);
  return current_error;
}

static void mc2d_AC_initialize(MohrCoulomb2DProblem* problem,
                               MohrCoulomb2DProblem* localproblem, SolverOptions* options) {
  /** In initialize, these operators are "connected" to their corresponding static variables,
   * that will be used to build local problem for each considered cone.
   * Local problem is built during call to update (which depends on the storage type for M).
   */

  DEBUG_PRINTF(
      "mc2d_AC_initialize starts with options->iparam[PLASTICITY_NSN_FORMULATION] = "
      "%i\n",
      options->iparam[PLASTICITY_NSN_FORMULATION]);

  if (options->iparam[PLASTICITY_NSN_FORMULATION] ==
      PLASTICITY_NSN_FORMULATION_ALARTCURNIER_STD) {
    numerics_printf(
        "PLASTICITY_NSN_FORMULATION =  PLASTICITY_NSN_FORMULATION_ALARTCURNIER_STD");
    Function = &(mc2d_computeAlartCurnierSTD);
  } else if (options->iparam[PLASTICITY_NSN_FORMULATION] ==
             PLASTICITY_NSN_FORMULATION_NATURALMAP) {
    numerics_printf("PLASTICITY_NSN_FORMULATION =  PLASTICITY_NSN_FORMULATION_NATURALMAP");
    Function = &(mc2d_computeNaturalMap);
  } else if (options->iparam[PLASTICITY_NSN_FORMULATION] == PLASTICITY_NSN_FORMULATION_NULL) {
    Function = NULL;
  }

  /* Compute and store default value of rho value */
  size_t nc = problem->numberOfCones;

  double avg_rho[3] = {0.0, 0.0, 0.0};

  if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN ||
      options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP) {
    if (!options->dWork || options->dWorkSize < 3 * nc) {
      options->dWork = (double*)realloc(options->dWork, 3 * nc * sizeof(double));
      options->dWorkSize = 3 * nc;
    }
  } else if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP_HYBRID) {
    if (!options->dWork || options->dWorkSize < 4 * nc) {
      options->dWork = (double*)realloc(options->dWork, 4 * nc * sizeof(double));
      options->dWorkSize = 4 * nc;
    }
  }

  double* rho = 0;
  for (size_t cone = 0; cone < nc; cone++) {
    if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN ||
        options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP) {
      rho = &options->dWork[3 * cone];
    } else if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP_HYBRID) {
      options->dWork[cone] = 1.0;  // for PLI algorithm.
      rho = &options->dWork[3 * cone + nc];
    }
    numerics_printf(
        "mc2d_AC_initialize "
        " mc2d_compute rho for cone = %i",
        cone);

    if (options->iparam[PLASTICITY_NSN_RHO_STRATEGY] ==
        PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND) {
      mc2d_local_problem_fill_M(problem, localproblem, cone);
      mc2d_compute_rho_split_spectral_norm_cond(localproblem, rho);
    } else if (options->iparam[PLASTICITY_NSN_RHO_STRATEGY] ==
               PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM) {
      mc2d_local_problem_fill_M(problem, localproblem, cone);
      mc2d_compute_rho_split_spectral_norm(localproblem, rho);
    } else if (options->iparam[PLASTICITY_NSN_RHO_STRATEGY] ==
               PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM) {
      mc2d_local_problem_fill_M(problem, localproblem, cone);
      mc2d_compute_rho_spectral_norm(localproblem, rho);
    } else if (options->iparam[PLASTICITY_NSN_RHO_STRATEGY] ==
               PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_CONSTANT) {
      rho[0] = options->dparam[PLASTICITY_NSN_RHO];
      rho[1] = options->dparam[PLASTICITY_NSN_RHO];
      rho[2] = options->dparam[PLASTICITY_NSN_RHO];
    } else if (options->iparam[PLASTICITY_NSN_RHO_STRATEGY] ==
               PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE) {
      numerics_error("mc2d_AC_initialize",
                     "Adaptive strategy for computing rho not yet implemented");
    } else
      numerics_error("mc2d_AC_initialize", "unknown strategy for computing rho");

    if (verbose > 0) {
      avg_rho[0] += rho[0];
      avg_rho[1] += rho[1];
      avg_rho[2] += rho[2];
    }
    numerics_printf(
        "mc2d_AC_initialize"
        "cone = %i, rho[0] = %4.2e, rho[1] = %4.2e, rho[2] = %4.2e",
        cone, rho[0], rho[1], rho[2]);

    mc2d_local_problem_fill_M(problem, localproblem, cone);
    double m_row_norm = 0.0, sum;
    for (int i = 0; i < 3; i++) {
      sum = 0.0;
      for (int j = 0; j < 3; j++) {
        sum += fabs(localproblem->M->matrix0[i + j * 3]);
      }
      m_row_norm = max(sum, m_row_norm);
    }
    numerics_printf(
        "mc2d_AC_initialize"
        " inverse of norm of M = %e",
        1.0 / hypot9(localproblem->M->matrix0));
    numerics_printf(
        "mc2d_AC_initialize"
        " inverse of row norm of M = %e",
        1.0 / m_row_norm);

    DEBUG_EXPR(NM_display(localproblem->M););
  }
  numerics_printf(
      "mc2d_AC_initialize"
      " Avg. rho value = %e\t%e\t%e\t",
      avg_rho[0] / nc, avg_rho[1] / nc, avg_rho[2] / nc);
}

static void mc2d_AC_free(MohrCoulomb2DProblem* problem, MohrCoulomb2DProblem* localproblem,
                         SolverOptions* localsolver_options) {
  F = NULL;
  jacobianF = NULL;
  free(localsolver_options->dWork);
  localsolver_options->dWork = NULL;
}

void mc2d_onecone_nonsmooth_Newton_solvers_initialize(MohrCoulomb2DProblem* problem,
                                                      MohrCoulomb2DProblem* localproblem,
                                                      SolverOptions* localsolver_options) {
  /* Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen
   * formulation. */

  /* Alart-Curnier formulation */
  if (localsolver_options->solverId == MOHR_COULOMB_2D_ONECONE_NSN ||
      localsolver_options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP ||
      localsolver_options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP_HYBRID) {
    mc2d_AC_initialize(problem, localproblem, localsolver_options);
  } else {
    numerics_error("mc2d_onecone_nonsmooth_Newton_solvers_initialize",
                   "Unknown formulation type.");
  }
}

int mc2d_onecone_nonsmooth_Newton_solvers_solve(MohrCoulomb2DProblem* localproblem,
                                                double* local_reaction,
                                                SolverOptions* options) {
  numerics_printf_verbose(
      2, "--------------- mc2d_onecone_nonsmooth_Newton_solvers_solve starts for cone %i",
      options->iparam[PLASTICITY_CURRENT_CONE_NUMBER]);

  int info = 1;

  /*  check trivial solution */

  double* q = localproblem->q;
  if (q[0] > -DBL_EPSILON) {
    local_reaction[0] = 0.0;
    local_reaction[1] = 0.0;
    local_reaction[2] = 0.0;
    numerics_printf_verbose(2, "take off, trivial solution reaction = 0, velocity = q.");
    info = 0;
    options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
    numerics_printf_verbose(
        2, "--------------- mc2d_onecone_nonsmooth_Newton_solvers_solve ends");
    return info;
  }

  if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN) {
    info = mc2d_onecone_nonsmooth_Newton_solvers_solve_direct(localproblem, local_reaction,
                                                              options);
  } else if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP) {
    info = mc2d_onecone_nonsmooth_Newton_solvers_solve_damped(localproblem, local_reaction,
                                                              options);
  } else if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP_HYBRID) {
    if (options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] ==
            PLASTICITY_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP ||
        options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] ==
            PLASTICITY_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP) {
      info = mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid(localproblem, local_reaction,
                                                                options);
    } else {
      numerics_error("mc2d_onecone_nonsmooth_Newton_solvers_solve",
                     "Unknown local nsn hybrid solver");
    }
  } else {
    info = nonSmoothDirectNewton(Fsize, local_reaction, &F, &jacobianF, options);
  }
  if (info > 0) {
    if (verbose > 0) {
      if (options->iparam[SICONOS_IPARAM_MAX_ITER] ==
          options->iparam[SICONOS_IPARAM_ITER_DONE]) {
        numerics_warning(
            "mc2d_onecone_nonsmooth_Newton_solvers_solve",
            "reached max. number of iterations (%i) for cone %i. Residual = %12.8e",
            options->iparam[SICONOS_IPARAM_MAX_ITER],
            options->iparam[PLASTICITY_CURRENT_CONE_NUMBER],
            options->dparam[SICONOS_DPARAM_RESIDU]);
      } else {
        numerics_warning("mc2d_onecone_nonsmooth_Newton_solvers_solve",
                         "no convergence for cone %i with error = %12.8e",
                         options->iparam[SICONOS_IPARAM_MAX_ITER],
                         options->iparam[PLASTICITY_CURRENT_CONE_NUMBER],
                         options->dparam[SICONOS_DPARAM_RESIDU]);
      }
      /* note : exit on failure should be done in DefaultCheckSolverOutput */
    }
  }
  numerics_printf_verbose(2,
                          "--------------- mc2d_onecone_nonsmooth_Newton_solvers_solve ends");
  return info;
  /*  (*postSolver)(cone,reaction); */
}

void mc2d_onecone_nonsmooth_Newton_solvers_free(MohrCoulomb2DProblem* problem,
                                                MohrCoulomb2DProblem* localproblem,
                                                SolverOptions* localsolver_options) {
  F = NULL;
  jacobianF = NULL;
  if (localsolver_options->solverId == MOHR_COULOMB_2D_ONECONE_NSN ||
      localsolver_options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP ||
      localsolver_options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP_HYBRID) {
    mc2d_AC_free(problem, localproblem, localsolver_options);
  } else {
    numerics_error("mc2d_onecone_nonsmooth_Newton_solvers_initialize",
                   "Unknown formulation type.");
  }
}

void mc2d_onecone_nonsmooth_Newton_solvers_computeError(int n, double* velocity,
                                                        double* reaction, double* error) {
  /*   int numberOfCones = n/3; */
  /*   int sizeGlobal = numberOfCones*FSize; */
  /*   //  double * FGlobal = (double*)malloc(sizeGlobal*sizeof(*FGlobal));  */
  /*   (*computeFGlobal)(reaction,velocity); */
  /*   int i; */
  /*   double Fz; */
  /*   *error = 0; */
  /*   for(i=0;i<sizeGlobal;++i) */
  /*     { */
  /*       Fz = velocity[i]*reaction[i]; */
  /*       if(Fz>0) */
  /*  *error+=Fz; */
  /*       if(reaction[i]<0) */
  /*  *error+=reaction[i]; */
  /*       if(velocity[i]<0) */
  /*  *error+=velocity[i]; */
  /*     } */

  /*   // (*computeVelocity)(FGlobal); */

  /*   free(FGlobal); */
}
#ifdef DEBUG_CHECK
static int mc2d_onecone_nonsmooth_Newton_AC_debug(double* R, double* velocity, double theta,
                                                  double* rho, double* MLocal, double* F,
                                                  double* A, double* B, double* AWplusB,
                                                  int* iparam) {
  double AWpB[9];
  if (iparam[PLASTICITY_NSN_FORMULATION] != PLASTICITY_NSN_FORMULATION_JEANMOREAU_GENERATED &&
      iparam[PLASTICITY_NSN_FORMULATION] != PLASTICITY_NSN_FORMULATION_NULL) {
    double Fg[3] = {0., 0., 0.};
    double Ag[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double Bg[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    assert(*rho > 0. && *(rho + 1) > 0. && *(rho + 2) > 0.);
    if (iparam[PLASTICITY_NSN_FORMULATION] == PLASTICITY_NSN_FORMULATION_ALARTCURNIER_STD) {
      mc2d_AlartCurnierFunctionGenerated(R, velocity, theta, rho, Fg, Ag, Bg);
    }

    if (iparam[PLASTICITY_NSN_FORMULATION] == PLASTICITY_NSN_FORMULATION_JEANMOREAU_STD) {
      mc2d_AlartCurnierJeanMoreauFunctionGenerated(R, velocity, theta, rho, Fg, Ag, Bg);
    }

    sub3(F, Fg);
    sub3x3(A, Ag);
    sub3x3(B, Bg);

    assert(hypot3(Fg) <= 1e-7);
    assert(hypot9(Ag) <= 1e-7);
    assert(hypot9(Bg) <= 1e-7);
    cpy3x3(A, Ag);
    cpy3x3(B, Bg);
    mm3x3(A, MLocal, AWpB);
    add3x3(B, AWpB);

    scal3x3(-1., AWpB);
    sub3x3(AWplusB, AWpB);
    assert(hypot9(AWpB) <= 1e-7);
  }

  return 0;
}
#endif

void mc2d_onecone_nonsmooth_Newton_AC_update(int cone, MohrCoulomb2DProblem* problem,
                                             MohrCoulomb2DProblem* localproblem,
                                             double* reaction, SolverOptions* options) {
  /* Build a local problem for a specific cone
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific cone
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  mc2d_local_problem_fill_M(problem, localproblem, cone);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.reactionBlock, excluding the block corresponding to the current cone. ****/
  mc2d_local_problem_compute_q(problem, localproblem, reaction, cone);

  /*  coefficient for current block*/
  localproblem->eta[0] = problem->eta[cone];
  localproblem->theta[0] = problem->theta[cone];
}

int mc2d_onecone_nonsmooth_Newton_solvers_solve_direct(MohrCoulomb2DProblem* localproblem,
                                                       double* R, SolverOptions* options) {
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* numerics_printf_verbose( */
  /*     2, "--------------- mc2d_onecone_nonsmooth_Newton_solvers_solve_direct starts"); */

  double theta = localproblem->theta[0];
  double eta = localproblem->eta[0];
  double* qLocal = localproblem->q;

  double norm_qLocal =
      sqrt(qLocal[0] * qLocal[0] + qLocal[1] * qLocal[1] + qLocal[2] * qLocal[2]);
  double norm_relative = 1.0;
  if (norm_qLocal > DBL_EPSILON) {
    norm_relative /= norm_qLocal;
  }

  double* MLocal = localproblem->M->matrix0;

  /* store the increment */
  double dR[3] = {0., 0., 0.};
  /* store the value fo the function */
  double F[3] = {0., 0., 0.};
  /* Store the (sub)-gradient of the function */
  double A[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double B[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  /* Value of AW+B */
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  /* retrieve value of rho */
  double* rho;
  if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN ||
      options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP) {
    rho = &options->dWork[3 * iparam[PLASTICITY_CURRENT_CONE_NUMBER]];
  } else {
    int nc = options->dWorkSize / 4;
    rho = &options->dWork[3 * iparam[PLASTICITY_CURRENT_CONE_NUMBER] + nc];
  }

  /* compute the velocity */
  double velocity[3] = {0., 0., 0.};
  cpy3(qLocal, velocity);
  mvp3x3(MLocal, R, velocity);

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];

  /* Newton iteration */
  int inew = 0;
  int info_solv3x3;

  /* compute first residue */
  Function(R, velocity, eta, theta, rho, F, NULL, NULL);
  dparam[SICONOS_DPARAM_RESIDU] =
      0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) * norm_relative;
  dparam[SICONOS_DPARAM_RESIDU] =
      sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) * norm_relative;

  numerics_printf_verbose(
      2,
      "---------------    mc2d_onecone_nonsmooth_Newton_solvers_solve_direct  -- "
      "cone %i # iteration = %i  merit = %.10e",
      iparam[PLASTICITY_CURRENT_CONE_NUMBER], inew, dparam[SICONOS_DPARAM_RESIDU]);
  if (dparam[SICONOS_DPARAM_RESIDU] < dparam[SICONOS_DPARAM_TOL]) {
    iparam[SICONOS_IPARAM_ITER_DONE] = inew;
    return 0;
  }

  for (inew = 1; inew < itermax; ++inew) {
    /* Update function and gradient */
    Function(R, velocity, eta, theta, rho, F, A, B);

    /* compute -(A MLocal +B) */
    mm3x3(A, MLocal, AWplusB);
    add3x3(B, AWplusB);
    scal3x3(-1., AWplusB);

#ifdef DEBUG_CHECK
    mc2d_onecone_nonsmooth_Newton_AC_debug(R, velocity, theta, rho, MLocal, F, A, B, AWplusB,
                                           iparam);
#endif

    /* Solve the linear system */
    cpy3(F, dR);
    info_solv3x3 = solve_3x3_gepp(AWplusB, dR);

    /* if determinant is zero, replace dR=NaN with zero (i.e. don't modify R) and return early
     */
    if (info_solv3x3) {
      dR[0] = 0;
      dR[1] = 0;
      dR[2] = 0;
    }

    /* update iterates */
    add3(dR, R);

    /* compute new residue */
    cpy3(qLocal, velocity);
    mvp3x3(MLocal, R, velocity);
    Function(R, velocity, eta, theta, rho, F, NULL, NULL);
    /* dparam[SICONOS_DPARAM_RESIDU] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) * */
    /*                                 norm_relative;  // improve with relative tolerance */
    dparam[SICONOS_DPARAM_RESIDU] = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) *
                                    norm_relative;  // improve with relative tolerance
    /* double error = mc2d_compute_local_error(localproblem, R); */
    /* printf("New local error = %e\n", error); */
    if (info_solv3x3) {
      if (verbose > 0)
        numerics_warning("mc2d_onecone_nonsmooth_Newton_solvers_solve_direct",
                         "cone %i 3x3 linear system is irregular # iteration = %"
                         "merit  = %.10e",
                         iparam[PLASTICITY_CURRENT_CONE_NUMBER], inew,
                         dparam[SICONOS_DPARAM_RESIDU]);
      break;
    }

    numerics_printf_verbose(
        2,
        "---------------    mc2d_onecone_nonsmooth_Newton_solvers_solve_direct  -- "
        "cone %i # iteration = %i  merit = %.10e",
        iparam[PLASTICITY_CURRENT_CONE_NUMBER], inew, dparam[SICONOS_DPARAM_RESIDU]);

    if (dparam[SICONOS_DPARAM_RESIDU] < dparam[SICONOS_DPARAM_TOL]) {
      iparam[SICONOS_IPARAM_ITER_DONE] = inew;
      return 0;
    }
  }  // End of the Newton iteration

  iparam[SICONOS_IPARAM_ITER_DONE] = inew;
  return 1;
}

#undef DEBUG_MESSAGES
static int LineSearchGP(MohrCoulomb2DProblem* localproblem, computeNonsmoothFunction Function,
                        double* t_opt, double R[3], double dR[3], double* rho, int LSitermax,
                        double* F, double* A, double* B, double* velocity) {
  numerics_printf_verbose(2, "-- LineSearchGP starts");

  double alpha = *t_opt;

  double inf = 1e20;

  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.1, m2 = 0.9;

  /*     // store the value fo the function */
  /*     double F[3]={0.,0.,0.}; */

  /*     // Store the (sub)-gradient of the function */
  /*     double A[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.}; */
  /*     double B[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.}; */

  /*     double velocity[3]={0.,0.,0.}; */

  double eta = localproblem->eta[0];
  double theta = localproblem->theta[0];
  double* qLocal = localproblem->q;
  double* MLocal = localproblem->M->matrix0;

  /*     for (int i=0; i<3; i++) velocity[i] = MLocal[i+0*3]*R[0] + qLocal[i] */
  /*          + MLocal[i+1*3]*R[1] + */
  /*          + MLocal[i+2*3]*R[2] ; */

  /*     Function(R,velocity,theta,rho,F,A,B); */

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * cblas_ddot(3, F, 1, F, 1);

  double tmp[3] = {0., 0., 0.};

  // Value of AW+B
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // compute A MLocal +B
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      AWplusB[i + 3 * j] = 0.0;
      for (int k = 0; k < 3; k++) {
        AWplusB[i + 3 * j] += A[i + 3 * k] * MLocal[k + j * 3];
      }
      AWplusB[i + 3 * j] += B[i + 3 * j];
    }
  }

#ifdef DEBUG_MESSAGES
  for (int l = 0; l < 3; l++) {
    for (int k = 0; k < 3; k++) {
      printf("AWplusB[%i+3*%i] = %le\t", l, k, AWplusB[l + 3 * k]);
    }
    printf("\n");
  }
#endif

  for (int i = 0; i < 3; i++) {
    tmp[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      tmp[i] += AWplusB[i + 3 * j] * dR[j];
    }
  }

  double dqdt0 = 0.0;
  for (int i = 0; i < 3; i++) {
    dqdt0 += F[i] * tmp[i];
  }
#ifdef DEBUG_MESSAGES
  printf("q0 = %12.8e \n", q0);
  printf("dqdt0 = %12.8e \n", dqdt0);
  for (int i = 0; i < 3; i++) {
    printf("tmp[%i] = %12.8e \t", i, tmp[i]);
  }
  printf("\n");
  for (int i = 0; i < 3; i++) {
    printf("dR[%i] = %12.8e \t", i, dR[i]);
  }
  printf("\n");
#endif

  for (int iter = 0; iter < LSitermax; iter++) {
    for (int i = 0; i < 3; i++) tmp[i] = R[i] + alpha * dR[i];

    for (int i = 0; i < 3; i++)
      velocity[i] = MLocal[i + 0 * 3] * tmp[0] + qLocal[i] + MLocal[i + 1 * 3] * tmp[1] +
                    +MLocal[i + 2 * 3] * tmp[2];

    Function(tmp, velocity, eta, theta, rho, F, NULL, NULL);

    double q = 0.5 * cblas_ddot(3, F, 1, F, 1);

    double slope = (q - q0) / alpha;

#ifdef DEBUG_MESSAGES
    printf("q = %12.8e \n", q);
    printf("slope = %12.8e \n", slope);
#endif

    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if (C1 && C2) {
      DEBUG_PRINTF("Success in LS: alpha = %12.8e\n", alpha);
      *t_opt = alpha;

      numerics_printf_verbose(
          2, "-- LineSearchGP success number of iteration = %i  alpha = %.10e", iter, alpha);
      return 0;

    } else if (!C1) {
#ifdef DEBUG_MESSAGES
      printf("LS: alpha too small = %12.8e\t, slope =%12.8e\n", alpha, slope);
      printf(" m1*dqdt0 =%12.8e\t, m2*dqdt0 =%12.8e\n ", m1 * dqdt0, m2 * dqdt0);
#endif
      // std::cout << "t = " << t << " is too small : slope = " << slope << ", m2*qp0 = " <<
      // m2*qp0 << std::endl;
      alphamin = alpha;
    } else  // not(C2)
    {
#ifdef DEBUG_MESSAGES
      printf("LS: alpha too big = %12.8e\t, slope =%12.8e\n", alpha, slope);
      printf(" m1*dqdt0 =%12.8e\t, m2*dqdt0 =%12.8e\n ", m1 * dqdt0, m2 * dqdt0);
#endif
      // std::cout << "t = " << t << " is too big : slope = " << slope << ", m1*qp0 = " <<
      // m1*qp0 << std::endl;
      alphamax = alpha;
    }
    if (alpha < inf) {
      alpha = 0.5 * (alphamin + alphamax);
    } else {
      alpha = 10 * alpha;
    }
  }
  numerics_printf_verbose(
      2,
      "-- LineSearchGP not succeed max number of iteration reached  = %i  with alpha = %.10e",
      LSitermax, alpha);
  *t_opt = alpha;
  return -1;
}
#define DEBUG_MESSAGES
int mc2d_onecone_nonsmooth_Newton_solvers_solve_damped(MohrCoulomb2DProblem* localproblem,
                                                       double* R, SolverOptions* options) {
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  numerics_printf_verbose(
      2, "--------------- mc2d_onecone_nonsmooth_Newton_solvers_solve_damped starts");

  double eta = localproblem->eta[0];
  double theta = localproblem->theta[0];
  double* qLocal = localproblem->q;

  double norm_qLocal =
      sqrt(qLocal[0] * qLocal[0] + qLocal[1] * qLocal[1] + qLocal[2] * qLocal[2]);
  double norm_relative = 1.0;
  if (norm_qLocal > DBL_EPSILON) {
    norm_relative /= norm_qLocal;
  }

  double* MLocal = localproblem->M->matrix0;

  /* store the increment */
  double dR[3] = {0., 0., 0.};
  /* store the value fo the function */
  double F[3] = {0., 0., 0.};
  /* Store the (sub)-gradient of the function */
  double A[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double B[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  /* Value of AW+B */
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  /* retrieve value of rho */
  double* rho;
  if (options->solverId == MOHR_COULOMB_2D_ONECONE_NSN ||
      options->solverId == MOHR_COULOMB_2D_ONECONE_NSN_GP) {
    rho = &options->dWork[3 * iparam[PLASTICITY_CURRENT_CONE_NUMBER]];
  } else {
    int nc = options->dWorkSize / 4;
    rho = &options->dWork[3 * iparam[PLASTICITY_CURRENT_CONE_NUMBER] + nc];
  }

  /* compute the velocity */
  double velocity[3] = {0., 0., 0.};
  cpy3(qLocal, velocity);
  mvp3x3(MLocal, R, velocity);

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];

  /* Newton iteration */
  int inew;
  int info_solv3x3;
  double t = 1.;
  double t_opt = 1.;
  double t_init = 1.;

  int LSitermax = iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER];
  for (inew = 0; inew < itermax; ++inew) {
    /* Update function and gradient */
    Function(R, velocity, eta, theta, rho, F, A, B);

    /* compute -(A MLocal +B) */
    mm3x3(A, MLocal, AWplusB);
    add3x3(B, AWplusB);
    scal3x3(-1., AWplusB);

#ifdef DEBUG_CHECK
    mc2d_onecone_nonsmooth_Newton_AC_debug(R, velocity, theta, rho, MLocal, F, A, B, AWplusB,
                                           iparam);
#endif

    /* Solve the linear system */
    cpy3(F, dR);
    info_solv3x3 = solve_3x3_gepp(AWplusB, dR);

    /* if determinant is zero, replace dR=NaN with zero (i.e. don't modify R) and return early
     */
    if (info_solv3x3) {
      dR[0] = 0;
      dR[1] = 0;
      dR[2] = 0;
    }
    if (!info_solv3x3) {
      /* Perform Line Search */

      t_opt = t_init;
      LineSearchGP(localproblem, Function, &t_opt, R, dR, rho, LSitermax, F, A, B, velocity);
      t = t_opt;
    }

    /* update iterates */
    R[0] = R[0] + t * dR[0];
    R[1] = R[1] + t * dR[1];
    R[2] = R[2] + t * dR[2];

    /* compute new residue */
    cpy3(qLocal, velocity);
    mvp3x3(MLocal, R, velocity);
    Function(R, velocity, eta, theta, rho, F, NULL, NULL);
    /* dparam[SICONOS_DPARAM_RESIDU] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) * */
    /*                                 norm_relative;  // improve with relative tolerance */
    dparam[SICONOS_DPARAM_RESIDU] = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) *
                                    norm_relative;  // improve with relative tolerance

    if (info_solv3x3) {
      if (verbose > 0)
        numerics_warning("mc2d_onecone_nonsmooth_Newton_solvers_solve_damped",
                         "cone %i 3x3 linear system is irregular # iteration = %i"
                         " error = %.10e \n",
                         iparam[PLASTICITY_CURRENT_CONE_NUMBER], inew,
                         dparam[SICONOS_DPARAM_RESIDU]);
      break;
    }
    numerics_printf_verbose(2, "-- cone %i # iteration = %i  error = %.10e",
                            iparam[PLASTICITY_CURRENT_CONE_NUMBER], inew,
                            dparam[SICONOS_DPARAM_RESIDU]);

    if (dparam[SICONOS_DPARAM_RESIDU] < dparam[SICONOS_DPARAM_TOL]) {
      iparam[SICONOS_IPARAM_ITER_DONE] = inew;
      return 0;
    }
  }  // End of the Newton iteration

  iparam[SICONOS_IPARAM_ITER_DONE] = inew;
  return 1;
}

static void keep_or_discard_solution(MohrCoulomb2DProblem* localproblem,
                                     double* local_reaction, double* local_reaction_backup,
                                     SolverOptions* options, double* current_error) {
  double error = 0.0;
  /* error = options->dparam[SICONOS_DPARAM_RESIDU]; */

  error = mc2d_compute_local_error(localproblem, local_reaction);
  DEBUG_PRINTF("New local error = %e\n", error);
  int nan = isnan(options->dparam[SICONOS_DPARAM_RESIDU]) ||
            isinf(options->dparam[SICONOS_DPARAM_RESIDU]);
  if (nan) {
    DEBUG_PRINT("Residu is equal to nan or inf\n");
    DEBUG_EXPR(NM_display(localproblem->M));
    DEBUG_EXPR(NM_vector_display(localproblem->q, 3));
    DEBUG_EXPR(NM_vector_display(local_reaction, 3));
    DEBUG_EXPR(NM_vector_display(local_reaction_backup, 3));

    /* DEBUG_PRINTF("No hope for cone %d, setting to zero.\n", */
    /*              options->iparam[PLASTICITY_CURRENT_CONE_NUMBER]); */
    /* local_reaction[0] = 0; */
    /* local_reaction[1] = 0; */
    /* local_reaction[2] = 0; */
    DEBUG_PRINTF("Discard the new local solution with error = %e\n", error);
    DEBUG_PRINTF("Get back to the local backup solution = %e\n", *current_error);
    cpy3(local_reaction_backup, local_reaction);
    // memcpy(local_reaction, local_reaction_backup, sizeof(double)*3);
  } else {
    if (error <= options->dparam[SICONOS_DPARAM_TOL] || error <= *current_error) {
      DEBUG_PRINTF("Keep the new local solution with error = %e\n", error);
      *current_error = error;
      cpy3(local_reaction, local_reaction_backup);
      // memcpy(local_reaction_backup, local_reaction, sizeof(double)*3);
    } else {
      DEBUG_PRINTF("Discard the new local solution with error = %e\n", error);
      DEBUG_PRINTF("Get back to the local backup solution = %e\n", *current_error);
      cpy3(local_reaction_backup, local_reaction);
      // memcpy(local_reaction, local_reaction_backup, sizeof(double)*3);
    }
  }
}

int mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid(MohrCoulomb2DProblem* localproblem,
                                                       double* local_reaction,
                                                       SolverOptions* options) {
  numerics_printf_verbose(
      2, "--------------- mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid starts");

  int info = -1;
  double local_reaction_backup[3] = {local_reaction[0], local_reaction[1], local_reaction[2]};

  int max_loop = options->iparam[PLASTICITY_NSN_HYBRID_MAX_LOOP];

  int newton_iteration_number = options->iparam[SICONOS_IPARAM_MAX_ITER];
  int pli_iteration_number = options->iparam[PLASTICITY_NSN_HYBRID_MAX_ITER];

  /* compute current_error */
  double current_error = mc2d_compute_local_error(localproblem, local_reaction);

  DEBUG_PRINTF("mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid current_error= %e\n",
               current_error);

  if (!(options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] ==
            PLASTICITY_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP ||
        options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] ==
            PLASTICITY_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP)) {
    numerics_error("mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid",
                   "Unknown local nsn hybrid solver");
  }

  /* 0 - Perform a first call to NSN solver to see if it succeeds quickly */

  if (options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] ==
      PLASTICITY_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP) {
    options->iparam[SICONOS_IPARAM_MAX_ITER] = newton_iteration_number;
    info = mc2d_onecone_nonsmooth_Newton_solvers_solve_damped(localproblem, local_reaction,
                                                              options);

    DEBUG_PRINTF("NSN solver ended with residual = %e\n",
                 options->dparam[SICONOS_DPARAM_RESIDU]);

    keep_or_discard_solution(localproblem, local_reaction, local_reaction_backup, options,
                             &current_error);

    if (current_error <= options->dparam[SICONOS_DPARAM_TOL]) {
      options->dparam[SICONOS_DPARAM_RESIDU] = current_error;
      DEBUG_PRINTF("First call of NSN solver ends with current_error= %e\n", current_error);
      DEBUG_END("mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid\n");
      return info;
    }
  }

  if (options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] ==
          PLASTICITY_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP ||
      options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] ==
          PLASTICITY_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP) {
    int loop = 0;
    while (loop < max_loop && current_error >= options->dparam[SICONOS_DPARAM_TOL]) {
      loop++;
      DEBUG_PRINTF("PLASTICITY_ONECONE_NSN_GP_HYBRID:  loop = %i\n", loop);

      /*  1 - fixed point projection solver */

      options->iparam[SICONOS_IPARAM_MAX_ITER] = pli_iteration_number;
      info =
          mc2d_projectionOnConeWithLocalIteration_solve(localproblem, local_reaction, options);
      DEBUG_PRINTF("PLI solver ended with residual = %e\n",
                   options->dparam[SICONOS_DPARAM_RESIDU]);

      keep_or_discard_solution(localproblem, local_reaction, local_reaction_backup, options,
                               &current_error);
      if (current_error < options->dparam[SICONOS_DPARAM_TOL]) {
        options->iparam[SICONOS_IPARAM_MAX_ITER] = newton_iteration_number;
        break;
      }

      /* 2 - nonsmooth Newton solver */

      options->iparam[SICONOS_IPARAM_MAX_ITER] = newton_iteration_number;
      info = mc2d_onecone_nonsmooth_Newton_solvers_solve_damped(localproblem, local_reaction,
                                                                options);
      DEBUG_PRINTF("NSN solver  ended with residual = %e\n",
                   options->dparam[SICONOS_DPARAM_RESIDU]);

      keep_or_discard_solution(localproblem, local_reaction, local_reaction_backup, options,
                               &current_error);
    }

    if (loop == max_loop && max_loop != 1) {
      if (verbose > 0) {
        printf(
            "Maximum number of loop (%i) in PLASTICITY_ONECONE_NSN_GP_HYBRID has "
            "been reached for cone %i with error = %e \n",
            max_loop, options->iparam[PLASTICITY_CURRENT_CONE_NUMBER], current_error);
      }
    }
  }

  options->dparam[SICONOS_DPARAM_RESIDU] = current_error;
  DEBUG_PRINTF("mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid END current_error= %e\n",
               current_error);
  DEBUG_END("mc2d_onecone_nonsmooth_Newton_solvers_solve_hybrid\n")
  return info;
}

void mc2d_onecone_nsn_set_default(SolverOptions* options) {
  /* Value of rho parameter */
  options->iparam[PLASTICITY_NSN_RHO_STRATEGY] =
      PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND;
  options->dparam[PLASTICITY_NSN_RHO] = 1.0;

  /* Choice of formulation */
  options->iparam[PLASTICITY_NSN_FORMULATION] = PLASTICITY_NSN_FORMULATION_NATURALMAP;

  /* Choice of line -search method */
  options->iparam[PLASTICITY_NSN_LINESEARCH] = PLASTICITY_NSN_LINESEARCH_NO;

  /* parameters for hybrid solvers */
  options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] =
      PLASTICITY_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP;

  options->iparam[PLASTICITY_NSN_HYBRID_MAX_LOOP] = 1;
  options->iparam[PLASTICITY_NSN_HYBRID_MAX_ITER] = 10;
}

void mc2d_onecone_nsn_gp_set_default(SolverOptions* options) {
  /* Value of rho parameter */
  options->iparam[PLASTICITY_NSN_RHO_STRATEGY] =
      PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM;
  options->dparam[PLASTICITY_NSN_RHO] = 1.0;

  /* Choice of formulation */
  options->iparam[PLASTICITY_NSN_FORMULATION] = PLASTICITY_NSN_FORMULATION_NATURALMAP;
  /* Choice of line -search method */
  options->iparam[PLASTICITY_NSN_LINESEARCH] = PLASTICITY_NSN_LINESEARCH_GOLDSTEINPRICE;
  options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER] = 10;

  /* parameters for hybrid solvers */
  options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] =
      PLASTICITY_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP;
  options->iparam[PLASTICITY_NSN_HYBRID_MAX_LOOP] = 1;
  options->iparam[PLASTICITY_NSN_HYBRID_MAX_ITER] = 100;
}
