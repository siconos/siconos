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
#include <assert.h>  // for assert
#include <stdio.h>   // for printf, fclose, fopen
#include <stdlib.h>  // for malloc, free, exit

#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContac...
#include "NumericsFwd.h"                   // for NumericsMatrix, Fric...
#include "NumericsVector.h"
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_d...
#include "SolverOptions.h"                       // for SICONOS_DPARAM_TOL
#include "fc3d_Solvers.h"                        // for fc3d_DeSaxceFixedPoint
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"  // for fc3d_nonsmooth_Newto...
#include "fc3d_short_names.h"
#include "gfc3d_Solvers.h"                       // for gfc3d_DeSaxceFixedPo...
#include "gfc3d_compute_error.h"
#include "numerics_verbose.h"
#include "utils/numerics_errors.h"
/* Solver registration system */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */

#include <float.h>
#include <string.h>
#include <time.h>

#include "gfc3d_ipm.h"      // for primalResidual, dualResidual ...
#include "siconos_debug.h"  // for DEBUG_EXPR, DEBUG_P...

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void gfc3d_nsgs_wr(GlobalFrictionContactProblem* problem, double* reaction, double* velocity,
                   double* globalVelocity, int* info, SolverOptions* options) {
  /* verbose = 1; */

  DEBUG_BEGIN("gfc3d_nsgs_wr\n");
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation
    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);

    DEBUG_EXPR(frictionContact_display(reduced_problem););
    if (verbose) {
      printf("Call to the fc3d solver ...\n");
    }
    // call nsgs solver for the local problem
    // long clk_tck = CLOCKS_PER_SEC;
    // clock_t t1 = clock();
    fc3d_nsgs(reduced_problem, reaction, velocity, info, options);
    // clock_t t2 = clock();
    // printf("\nTIME = %10.4f\n", (double)(t2 - t1) / (double)clk_tck);

    // // printf("\n\n Compute v of NSGS: \n");
    // FILE *sol_file = fopen("sol_data.res", "r");
    // for (int i=0; i < problem->numberOfContacts*3; i++)
    // {
    //   fscanf(sol_file, "%lf ", reaction+i);
    // }
    // fscanf(sol_file, "\n");
    // fclose(sol_file);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    // printf("\n\nNSGS v = "); printBlockVec(globalVelocity, problem->M->size0, 3, 0);
    // printf("\n\n");

    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;

    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
  DEBUG_END("gfc3d_nsgs_wr\n");
}

REGISTER_SOLVER(GFC3D_NSGS_WR,
                "GFC3D_NSGS_WR",
                "Non-smooth Gauss-Seidel for Global 3D Friction Contact with reduction",
                NULL,
                NULL,
                NULL,
                NULL,  /* error function */
                fc3d_nsgs_set_default,  /* set_default */
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */)

void gfc3d_admm_wr(GlobalFrictionContactProblem* problem, double* reaction, double* velocity,
                   double* globalVelocity, int* info, SolverOptions* options) {
  DEBUG_BEGIN("gfc3d_admm_wr\n");
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation

    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(reduced_problem););

    if (verbose) {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_admm(reduced_problem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;
    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
  DEBUG_END("gfc3d_admm_wr\n");
}
REGISTER_SOLVER(GFC3D_ADMM_WR,
                "GFC3D_ADMM_WR",
                "ADMM for Global 3D Friction Contact with reduction",
                NULL,
                NULL,
                NULL,
                NULL,  /* error function */
                fc3d_admm_set_default,  /* set_default */
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */)

void gfc3d_nonsmooth_Newton_AlartCurnier_wr(GlobalFrictionContactProblem* problem,
                                            double* reaction, double* velocity,
                                            double* globalVelocity, int* info,
                                            SolverOptions* options) {
  DEBUG_BEGIN("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n");
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation
    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(reduced_problem););

    numerics_printf("gfc3d_nonsmooth_Newton_AlartCurnier_wr - Call to the fc3d solver ...\n");

    fc3d_nonsmooth_Newton_AlartCurnier(reduced_problem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;

    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }

  DEBUG_END("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n")
}

REGISTER_SOLVER(GFC3D_NSN_AC_WR,
                "GFC3D_NSN_AC_WR",
                "Non-smooth Newton for Global 3D Friction Contact with reduction", NULL, NULL,
                NULL, NULL,           /* error function */
                fc3d_nsn_ac_set_default, /* set_default */
                1000,                 /* default_max_iter */
                1e-4,                 /* default_tol */
                0 /* is_local_solver */)

void gfc3d_nonsmooth_Newton_AlartCurnier_new_wr(GlobalFrictionContactProblem* problem,
                                                double* reaction, double* velocity,
                                                double* globalVelocity, int* info,
                                                SolverOptions* options) {
  DEBUG_BEGIN("gfc3d_nonsmooth_Newton_AlartCurnier_new_wr(...)\n");
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation
    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ... "
                            "this make take a while");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(reduced_problem););

    numerics_printf(
        "gfc3d_nonsmooth_Newton_AlartCurnier_new_wr - Call to the fc3d solver ...\n");

    fc3d_nonsmooth_Newton_AlartCurnier_new(reduced_problem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;
    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }

  DEBUG_END("gfc3d_nonsmooth_Newton_AlartCurnier_wr(...)\n")
}

void gfc3d_nsgs_velocity_wr(GlobalFrictionContactProblem* problem, double* reaction,
                            double* velocity, double* globalVelocity, int* info,
                            SolverOptions* options) {
  NumericsMatrix* H = problem->H;

  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation

    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ... "
                            "this make take a while");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);

    DEBUG_EXPR(frictionContact_display(reduced_problem););
    if (verbose) {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_nsgs_velocity(reduced_problem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;
    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
}
REGISTER_SOLVER(GFC3D_NSGSV_WR, "GFC3D_NSGSV_WR",
                "Global 3D Friction Contact with reduction", NULL, NULL,
                NULL, NULL,              /* error function */
                fc3d_nsgs_velocity_set_default, /* set_default */
                1000,                    /* default_max_iter */
                1e-4,                    /* default_tol */
                0 /* is_local_solver */)

void gfc3d_proximal_wr(GlobalFrictionContactProblem* problem, double* reaction,
                       double* velocity, double* globalVelocity, int* info,
                       SolverOptions* options) {
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation

    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);

    DEBUG_EXPR(frictionContact_display(reduced_problem););
    DEBUG_EXPR(frictionContact_display(reduced_problem););
    if (verbose) {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_proximal(reduced_problem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;
    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
}
REGISTER_SOLVER(GFC3D_PROX_WR, "GFC3D_PROX_WR",
                "Global 3D Friction Contact with reduction", NULL, NULL,
                NULL, NULL,              /* error function */
                fc3d_proximal_set_default, /* set_default */
                1000,                    /* default_max_iter */
                1e-4,                    /* default_tol */
                0 /* is_local_solver */)


void gfc3d_DeSaxceFixedPoint_wr(GlobalFrictionContactProblem* problem, double* reaction,
                                double* velocity, double* globalVelocity, int* info,
                                SolverOptions* options) {
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation

    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ... "
                            "this make take a while");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);

    DEBUG_EXPR(frictionContact_display(reduced_problem););

    if (verbose) {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_DeSaxceFixedPoint(reduced_problem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;
    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
}
REGISTER_SOLVER(GFC3D_DSFP_WR, "GFC3D_DSFP_WR",
                "Global 3D Friction Contact with reduction", NULL, NULL,
                NULL, NULL,              /* error function */
                fc3d_dsfp_set_default, /* set_default */
                1000,                    /* default_max_iter */
                1e-4,                    /* default_tol */
                0 /* is_local_solver */)

void gfc3d_TrescaFixedPoint_wr(GlobalFrictionContactProblem* problem, double* reaction,
                               double* velocity, double* globalVelocity, int* info,
                               SolverOptions* options) {
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation
    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ... "
                            "this make take a while");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);

    DEBUG_EXPR(frictionContact_display(reduced_problem););

    if (verbose) {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_TrescaFixedPoint(reduced_problem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;
    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
}
REGISTER_SOLVER(GFC3D_TFP_WR, "GFC3D_TFP_WR",
                "Global 3D Friction Contact with reduction", NULL, NULL,
                NULL, NULL,              /* error function */
                fc3d_tfp_set_default, /* set_default */
                1000,                    /* default_max_iter */
                1e-4,                    /* default_tol */
                0 /* is_local_solver */)
void gfc3d_ipm_snm_wr(GlobalFrictionContactProblem* problem, double* reaction,
                      double* velocity, double* globalVelocity, int* info,
                      SolverOptions* options) {
  verbose = 1;
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / 3);
  if (H->size1 > 0) {
    // Reformulation
    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ... "
                            "this make take a while");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reformulation_FrictionContact(problem);
    DEBUG_EXPR(frictionContact_display(reduced_problem););

    if (verbose) {
      printf("Call to the fc3d solver ...\n");
    }
    fc3d_IPM_SNM(reduced_problem, reaction, velocity, info, options);
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = 3 * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;
    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
}
REGISTER_SOLVER(GFC3D_IPM_SNM_WR, "GFC3D_DSFP_WRGFC3D_IPM_SNM_WR",
                "Global 3D Friction Contact with reduction", NULL, NULL,
                NULL, NULL,              /* error function */
                fc3d_ipm_snm_set_default, /* set_default */
                1000,                    /* default_max_iter */
                1e-4,                    /* default_tol */
                0 /* is_local_solver */)
