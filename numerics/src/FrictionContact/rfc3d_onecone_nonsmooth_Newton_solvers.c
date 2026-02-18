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
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs, isinf
#include <stdio.h>   // for NULL, printf
#include <stdlib.h>  // for calloc, realloc

#include "Friction_cst.h"
#include "NSSTools.h"         // for max
#include "NonSmoothNewton.h"  // for nonSmoothDirec...
#include "NumericsFwd.h"      // for SolverOptions
#include "NumericsMatrix.h"   // for NumericsMatrix
#include "RollingFrictionContactProblem.h"
#include "SiconosBlas.h"                             // for cblas_ddot
#include "SolverOptions.h"                           // for SolverOptions
#include "numerics_verbose.h"                        // for numerics_print...
#include "op5x5.h"                                   // for cpy3, mvp3x3
#include "op3x3.h"                                   // for cpy3, mvp3x3
#include "rfc3d_onecone_nonsmooth_Newton_solvers.h"  // for computeNonsmoo...
#include "rolling_fc2d_projection.h"
#include "rolling_fc3d_compute_error.h"
#include "rolling_fc3d_local_problem_tools.h"
#include "rolling_fc3d_projection.h"
#include "rolling_naturalmap_functions.h"

/* #define DEBUG_CHECK */
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#define DEBUG_STDOUT
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
  double norm_q = cblas_dnrm2(3, localproblem->q, 1);
  double* local_M = localproblem->M->matrix0;
  double worktmp[3];
  double local_velocity[3] = {0., 0., 0.};
  cpy5(local_q, local_velocity);
  mvp5x5(local_M, local_reaction, local_velocity);
  DEBUG_EXPR(NV_display(local_q, 3););
  DEBUG_EXPR(NV_display(local_velocity, 3););
  DEBUG_EXPR(NV_display(local_reaction, 3););

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

static void compute_rho_split_spectral_norm_cond(RollingFrictionContactProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->matrix0;
  assert(MLocal[0 + 0 * 5] > 0);

  DEBUG_EXPR(NM_dense_display(MLocal, 5, 5, 5););
  double sw = MLocal[1 + 1 * 5] + MLocal[2 + 2 * 5];
  double dw = sw * sw - 4.0 * (sw - MLocal[2 + 1 * 5] + MLocal[1 + 2 * 5]);
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
  
  sw = MLocal[3 + 3 * 5] + MLocal[4 + 4 * 5];
  dw = sw * sw - 4.0 * (sw - MLocal[4 + 1 * 5] + MLocal[3 + 2 * 5]);
  DEBUG_PRINTF("dw = %e\n", dw);
  if (dw > 0.0)
    dw = sqrt(dw);
  else
    dw = 0.0;

  rho[3] = 2.0 * (sw - dw) / ((sw + dw) * (sw + dw));
  rho[4] = rho[3];

  assert(rho[0] > 0);
  assert(rho[1] > 0);
  assert(rho[2] > 0);
  assert(rho[3] > 0);
  assert(rho[4] > 0);


  DEBUG_PRINTF("sw=%le\t  ", sw);
  DEBUG_PRINTF("dw=%le\n ", dw);
  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\t", rho[2]);
  DEBUG_PRINTF("rho[3]=%le\t", rho[3]);
  DEBUG_PRINTF("rho[4]=%le\n", rho[4]);
}

static void compute_rho_split_spectral_norm(RollingFrictionContactProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->matrix0;
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
  sw = MLocal[3 + 3 * 3] + MLocal[4 + 4 * 3];
  dw = sw * sw - 4.0 * (sw - MLocal[4 + 1 * 3] + MLocal[3 + 2 * 3]);
  DEBUG_PRINTF("dw = %e\n", dw);
  if (dw > 0.0)
    dw = sqrt(dw);
  else
    dw = 0.0;

  rho[3] = 2.0 / (sw + dw);
  rho[4] = rho[3];
  assert(rho[3] > 0);
  assert(rho[4] > 0);

  
  DEBUG_PRINTF("sw=%le\t  ", sw);
  DEBUG_PRINTF("dw=%le\n ", dw);
  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\t", rho[2]);
  DEBUG_PRINTF("rho[3]=%le\t", rho[3]);
  DEBUG_PRINTF("rho[4]=%le\n", rho[4]);
}

static void compute_rho_spectral_norm(RollingFrictionContactProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->matrix0;

  SET5X5(MLocal);  

  double Mlocal3x3[9] = {*MLocal00, *MLocal10, *MLocal20, *MLocal01, *MLocal11, *MLocal21, *MLocal03, *MLocal13, *MLocal23};

  
  double worktmp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eig[3] = {0.0, 0.0, 0.0};
  if (eig_3x3(MLocal, worktmp, eig)) numerics_printf("compute_rho_spectral_norm : failed");
  DEBUG_PRINTF("eig[0] = %4.2e, eig[1] = %4.2e, eig[2] = %4.2e", eig[0], eig[1], eig[2]);
  DEBUG_PRINTF("1/eig[0] = %4.2e, 1/eig[1] = %4.2e, 1/eig[2] = %4.2e", 1.0 / eig[0],
               1.0 / eig[1], 1.0 / eig[2]);
  rho[0] = 1.0 / eig[0];
  rho[1] = rho[0];
  rho[2] = rho[0];
  rho[3] = rho[0];
  rho[4] = rho[0];

 
  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\t", rho[2]);
  DEBUG_PRINTF("rho[3]=%le\t", rho[3]);
  DEBUG_PRINTF("rho[4]=%le\n", rho[4]);
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

  double avg_rho[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  // to be updated
  if (options->solverId == SICONOS_ONECONE_NSN ||
      options->solverId == SICONOS_ONECONE_NSN_GP) {
    if (!options->dWork || options->dWorkSize < 5 * nc) {
      options->dWork = (double*)realloc(options->dWork, 5 * nc * sizeof(double));
      options->dWorkSize = 5 * nc;
    }
  } else if (options->solverId == SICONOS_ONECONE_NSN_GP_HYBRID) {
    if (!options->dWork || options->dWorkSize < 6 * nc) {
      options->dWork = (double*)realloc(options->dWork, 4 * nc * sizeof(double));
      options->dWorkSize = 6 * nc;
    }
  }

  double* rho = 0;
  for (size_t cone = 0; cone < nc; cone++) {
    if (options->solverId == SICONOS_ONECONE_NSN ||
        options->solverId == SICONOS_ONECONE_NSN_GP) {
      rho = &options->dWork[5 * cone];
    } else if (options->solverId == SICONOS_ONECONE_NSN_GP_HYBRID) {
      options->dWork[cone] = 1.0;  // for PLI algorithm.
      rho = &options->dWork[5 * cone + nc];
    }
    numerics_printf_verbose(2,
                            "rfc3d_onecone_nonsmooth_Newton_initialize"
                            " rfc3d_compute rho for cone = %i",
                            cone);

    if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] ==
        SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND) {
      rolling_fc3d_local_problem_fill_M(problem, localproblem, cone);
      compute_rho_split_spectral_norm_cond(localproblem, rho);
    } else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] ==
               SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM) {
      rolling_fc3d_local_problem_fill_M(problem, localproblem, cone);
      compute_rho_split_spectral_norm(localproblem, rho);
    } else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] ==
               SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM) {
      rolling_fc3d_local_problem_fill_M(problem, localproblem, cone);
      compute_rho_spectral_norm(localproblem, rho);
    } else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] ==
               SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT) {
      rho[0] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
      rho[1] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
      rho[2] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
      rho[3] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
      rho[4] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
    } else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] ==
               SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE) {
      numerics_error("rfc3d_onecone_nonsmooth_Newton_initialize",
                     "Adaptive strategy for computing rho not yet implemented");
    } else
      numerics_error("rfc3d_onecone_nonsmooth_Newton_initialize",
                     "unknown strategy for computing rho");

    rolling_fc3d_local_problem_fill_M(problem, localproblem, cone);

    if (verbose > 0) {
      avg_rho[0] += rho[0];
      avg_rho[1] += rho[1];
      avg_rho[2] += rho[2];
      avg_rho[3] += rho[3];
      avg_rho[4] += rho[4];

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
                              "cone = %i, rho[0] = %4.2e, rho[1] = %4.2e, rho[2] = %4.2e, "
                              "||M||^-1 = %4.2e, ||M||_row^-1 = %4.2e ",
                              cone, rho[0], rho[1], rho[2],
                              1.0 / hypot9(localproblem->M->matrix0), 1.0 / m_row_norm);
    }
    DEBUG_EXPR(NM_display(localproblem->M););
  }
  numerics_printf(
      "rfc3d_onecone_nonsmooth_Newton_initialize"
      " Avg. rho value = %e\t%e\t%e\t%e\t%e\t",
      avg_rho[0] / nc, avg_rho[1] / nc, avg_rho[2] / nc, avg_rho[3] / nc, avg_rho[4] / nc);
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
                                                 double* local_reaction,
                                                 SolverOptions* options) {
  verbose=3;
  numerics_printf_verbose(
      2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve starts for cone %i",
      options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]);

  int info = 1;

  /*  check trivial solution */

  double* q = localproblem->q;
  if (q[0] > -DBL_EPSILON) {
    local_reaction[0] = 0.0;
    local_reaction[1] = 0.0;
    local_reaction[2] = 0.0;
    local_reaction[3] = 0.0;
    local_reaction[4] = 0.0;
    numerics_printf_verbose(2, "take off, trivial solution reaction = 0, velocity = q.");
    info = 0;
    options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
    numerics_printf_verbose(
        2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve ends");
    return info;
  }
  // to be updated
  if (options->solverId == SICONOS_ONECONE_NSN) {
    info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct(localproblem, local_reaction,
                                                               options);
  } else if (options->solverId == SICONOS_ONECONE_NSN_GP) {
    info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped(localproblem, local_reaction,
                                                               options);
  } else if (options->solverId == SICONOS_ONECONE_NSN_GP_HYBRID) {
    if (options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==
            SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP ||
        options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==
            SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP) {
      info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid(localproblem, local_reaction,
                                                                 options);
    } else {
      numerics_error("rfc3d_onecone_nonsmooth_Newton_solvers_solve",
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
            "rfc3d_onecone_nonsmooth_Newton_solvers_solve",
            "reached max. number of iterations (%i) for cone %i. Residual = %12.8e",
            options->iparam[SICONOS_IPARAM_MAX_ITER],
            options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER],
            options->dparam[SICONOS_DPARAM_RESIDU]);
      } else {
        numerics_warning("rfc3d_onecone_nonsmooth_Newton_solvers_solve",
                         "no convergence for cone %i with error = %12.8e",
                         options->iparam[SICONOS_IPARAM_MAX_ITER],
                         options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER],
                         options->dparam[SICONOS_DPARAM_RESIDU]);
      }
      /* note : exit on failure should be done in DefaultCheckSolverOutput */
    }
  }
  numerics_printf_verbose(2,
                          "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve ends");
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
    RollingFrictionContactProblem* localproblem, double* R, SolverOptions* options) {
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* numerics_printf_verbose( */
  /*     2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct starts"); */

  double mu = localproblem->mu[0];
  double mu_r = localproblem->mu_r[0];
  double* qLocal = localproblem->q;

  double norm_qLocal = sqrt(qLocal[0] * qLocal[0] + qLocal[1] * qLocal[1] +
                            qLocal[2] * qLocal[2] + qLocal[3] * qLocal[3] + qLocal[4] * qLocal[4]);
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
    rho = &options->dWork[5 * iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]];
  } else {
    int nc = options->dWorkSize / 4;
    rho = &options->dWork[5 * iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] + nc];
  }

  /* compute the velocity */
  double velocity[5] = {0., 0., 0., 0., 0.};
  cpy5(qLocal, velocity);
  mvp5x5(MLocal, R, velocity);

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];

  /* Newton iteration */
  int inew = 0;
  int info_solv5x5;

  /* compute first residue */
  Function(R, velocity, mu, mu_r, rho, F, NULL, NULL);
  dparam[SICONOS_DPARAM_RESIDU] =
      0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]+ F[3] * F[3] + F[4] * F[4]) * norm_relative;
  /* dparam[SICONOS_DPARAM_RESIDU] = */
  /*     sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) * norm_relative; */

  numerics_printf_verbose(
      2,
      "---------------    rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct  -- "
      "cone %i # iteration = %i  merit = %.10e",
      iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew, dparam[SICONOS_DPARAM_RESIDU]);
  if (dparam[SICONOS_DPARAM_RESIDU] < dparam[SICONOS_DPARAM_TOL]) {
    iparam[SICONOS_IPARAM_ITER_DONE] = inew;
    return 0;
  }

  for (inew = 1; inew < itermax; ++inew) {
    /* Update function and gradient */
    Function(R, velocity, mu, mu_r, rho, F, A, B);
    DEBUG_EXPR(printf("F"); display5(F));
    /* compute -(A MLocal +B) */
    mm5x5(A, MLocal, AWplusB);
    add5x5(B, AWplusB);
    //scal5x5(1., AWplusB);
    scal5(-1., F);    
    DEBUG_EXPR(printf("AWplusB"); display5x5(AWplusB));

    double AWplusB_save[25];
    cpy5x5(AWplusB, AWplusB_save);
    
    double F_save[5];
    cpy5(F, F_save);

    /* Solve the linear system */

    solve_nxn_gepp(5, AWplusB, F, dR);

    scal5(-1., F_save); 
    mvp5x5(AWplusB_save, dR, F_save);
    
    printf("hypot5 = %e \n", hypot5(F_save));
    
    
    DEBUG_EXPR(printf("dR"); display5(dR));
    DEBUG_EXPR(printf("F"); display5(F));
    /* update iterates */
    add5(dR, R);

    /* compute new residue */
    cpy5(qLocal, velocity);
    mvp5x5(MLocal, R, velocity);
    Function(R, velocity, mu, mu_r, rho, F, NULL, NULL);

    dparam[SICONOS_DPARAM_RESIDU] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]+ F[3] * F[3] + F[4] * F[4]) *
                                    norm_relative;  // improve with relative tolerance

    // VA 01/2025 : should be better but prevents convergemce due to rounding error
    /* dparam[SICONOS_DPARAM_RESIDU] = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2] )* */
    /* 					 norm_relative;  // improve with relative tolerance */
    /* double error = rfc3d_compute_local_error(localproblem, R); */
    /* printf("New local error = %e\n", error); */

    /* if (info_solv3x3) { */
    /*   if (verbose > 0) */
    /*     numerics_warning("rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct", */
    /*                      "cone %i 3x3 linear system is irregular # iteration = %" */
    /*                      "merit  = %.10e", */
    /*                      iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew, */
    /*                      dparam[SICONOS_DPARAM_RESIDU]); */
    /*   break; */
    /* } */

    numerics_printf_verbose(
        2,
        "---------------    rfc3d_onecone_nonsmooth_Newton_solvers_solve_direct  -- "
        "cone %i # iteration = %i  merit = %.10e",
        iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew,
        dparam[SICONOS_DPARAM_RESIDU]);

    if (dparam[SICONOS_DPARAM_RESIDU] < dparam[SICONOS_DPARAM_TOL]) {
      iparam[SICONOS_IPARAM_ITER_DONE] = inew;
      return 0;
    }
  }  // End of the Newton iteration

  iparam[SICONOS_IPARAM_ITER_DONE] = inew;
  return 1;
}

#undef DEBUG_MESSAGES
static int LineSearchGP(RollingFrictionContactProblem* localproblem,
                        computeNonsmoothFunction Function, double* t_opt, double R[3],
                        double dR[3], double* rho, int LSitermax, double* F, double* A,
                        double* B, double* velocity) {
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

  double mu = localproblem->mu[0];
  double mu_r = localproblem->mu_r[0];
  double* qLocal = localproblem->q;
  double* MLocal = localproblem->M->matrix0;

  /*     for (int i=0; i<3; i++) velocity[i] = MLocal[i+0*3]*R[0] + qLocal[i] */
  /*          + MLocal[i+1*3]*R[1] + */
  /*          + MLocal[i+2*3]*R[2] ; */

  /*     Function(R,velocity,mu_r,rho,F,A,B); */

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

    Function(tmp, velocity, mu, mu_r, rho, F, NULL, NULL);

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
int rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped(
    RollingFrictionContactProblem* localproblem, double* R, SolverOptions* options) {
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  numerics_printf_verbose(
      2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped starts");

  double mu = localproblem->mu[0];
  double mu_r = localproblem->mu_r[0];
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
  if (options->solverId == SICONOS_ONECONE_NSN ||
      options->solverId == SICONOS_ONECONE_NSN_GP) {
    rho = &options->dWork[3 * iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]];
  } else {
    int nc = options->dWorkSize / 4;
    rho = &options->dWork[3 * iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] + nc];
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

  int LSitermax = iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER];
  for (inew = 0; inew < itermax; ++inew) {
    /* Update function and gradient */
    Function(R, velocity, mu, mu_r, rho, F, A, B);

    /* compute -(A MLocal +B) */
    mm3x3(A, MLocal, AWplusB);
    add3x3(B, AWplusB);
    scal3x3(-1., AWplusB);


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
    Function(R, velocity, mu, mu_r, rho, F, NULL, NULL);
    dparam[SICONOS_DPARAM_RESIDU] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) *
                                    norm_relative;  // improve with relative tolerance
    /* dparam[SICONOS_DPARAM_RESIDU] = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) * */
    /*                                 norm_relative;  // improve with relative tolerance */

    if (info_solv3x3) {
      if (verbose > 0)
        numerics_warning("rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped",
                         "cone %i 3x3 linear system is irregular # iteration = %i"
                         " error = %.10e \n",
                         iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew,
                         dparam[SICONOS_DPARAM_RESIDU]);
      break;
    }
    numerics_printf_verbose(2, "-- cone %i # iteration = %i  error = %.10e",
                            iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER], inew,
                            dparam[SICONOS_DPARAM_RESIDU]);

    if (dparam[SICONOS_DPARAM_RESIDU] < dparam[SICONOS_DPARAM_TOL]) {
      iparam[SICONOS_IPARAM_ITER_DONE] = inew;
      return 0;
    }
  }  // End of the Newton iteration

  iparam[SICONOS_IPARAM_ITER_DONE] = inew;
  return 1;
}

static void keep_or_discard_solution(RollingFrictionContactProblem* localproblem,
                                     double* local_reaction, double* local_reaction_backup,
                                     SolverOptions* options, double* current_error) {
  double error = 0.0;
  /* error = options->dparam[SICONOS_DPARAM_RESIDU]; */

  error = rfc3d_compute_local_error(localproblem, local_reaction);
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
    /*              options->iparam[ SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]); */
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

int rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid(
    RollingFrictionContactProblem* localproblem, double* local_reaction,
    SolverOptions* options) {
  numerics_printf_verbose(
      2, "--------------- rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid starts");

  int info = -1;
  double local_reaction_backup[3] = {local_reaction[0], local_reaction[1], local_reaction[2]};

  int max_loop = options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP];

  int newton_iteration_number = options->iparam[SICONOS_IPARAM_MAX_ITER];
  int pli_iteration_number = options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER];

  /* compute current_error */
  double current_error = rfc3d_compute_local_error(localproblem, local_reaction);

  DEBUG_PRINTF("rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid current_error= %e\n",
               current_error);

  if (!(options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==
            SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP ||
        options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==
            SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP)) {
    numerics_error("rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid",
                   "Unknown local nsn hybrid solver");
  }

  /* 0 - Perform a first call to NSN solver to see if it succeeds quickly */

  if (options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==
      SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP) {
    options->iparam[SICONOS_IPARAM_MAX_ITER] = newton_iteration_number;
    info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped(localproblem, local_reaction,
                                                               options);

    DEBUG_PRINTF("NSN solver ended with residual = %e\n",
                 options->dparam[SICONOS_DPARAM_RESIDU]);

    keep_or_discard_solution(localproblem, local_reaction, local_reaction_backup, options,
                             &current_error);

    if (current_error <= options->dparam[SICONOS_DPARAM_TOL]) {
      options->dparam[SICONOS_DPARAM_RESIDU] = current_error;
      DEBUG_PRINTF("First call of NSN solver ends with current_error= %e\n", current_error);
      DEBUG_END("rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid\n");
      return info;
    }
  }

  if (options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==
          SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP ||
      options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==
          SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP) {
    int loop = 0;
    while (loop < max_loop && current_error >= options->dparam[SICONOS_DPARAM_TOL]) {
      loop++;
      DEBUG_PRINTF(" SICONOS_ONECONE_NSN_GP_HYBRID:  loop = %i\n", loop);

      /*  1 - fixed point projection solver */

      options->iparam[SICONOS_IPARAM_MAX_ITER] = pli_iteration_number;
      info = rolling_fc3d_projectionOnConeWithLocalIteration_solve(localproblem,
                                                                   local_reaction, options);
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
      info = rfc3d_onecone_nonsmooth_Newton_solvers_solve_damped(localproblem, local_reaction,
                                                                 options);
      DEBUG_PRINTF("NSN solver  ended with residual = %e\n",
                   options->dparam[SICONOS_DPARAM_RESIDU]);

      keep_or_discard_solution(localproblem, local_reaction, local_reaction_backup, options,
                               &current_error);
    }

    if (loop == max_loop && max_loop != 1) {
      if (verbose > 0) {
        printf(
            "Maximum number of loop (%i) in  SICONOS_ONECONE_NSN_GP_HYBRID has "
            "been reached for cone %i with error = %e \n",
            max_loop, options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER],
            current_error);
      }
    }
  }

  options->dparam[SICONOS_DPARAM_RESIDU] = current_error;
  DEBUG_PRINTF("rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid END current_error= %e\n",
               current_error);
  DEBUG_END("rfc3d_onecone_nonsmooth_Newton_solvers_solve_hybrid\n")
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
