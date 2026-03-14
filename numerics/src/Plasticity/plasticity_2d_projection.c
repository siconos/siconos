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

#include "plasticity_2d_projection.h"  // for plasticity_2d_projectionOnConeWithDiago...

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt
#include <stdio.h>   // for fprintf, printf, NULL, stderr
#include <stdlib.h>  // for calloc, free, exit, EXIT_FAILURE

#include "PlasticityProblem.h"      // for PlasticityProblem
#include "NumericsFwd.h"               // for SolverOptions, MohrCoulomb2D...
#include "NumericsMatrix.h"            // for NumericsMatrix, RawNumericsMatrix
#include "Plasticity_options.h"            // for PLASTICITY_NSGS_LOCAL...
#include "SiconosBlas.h"               // for cblas_ddot
#include "SolverOptions.h"             // for SolverOptions, solver_options_...
#include "SparseBlockMatrix.h"         // for SBM_row_prod
#include "plasticity_2d_compute_error.h"        // for plasticity_2d_Tresca_unitary_compute_an...
#include "plasticity_2d_local_problem_tools.h"  // for plasticity_2d_local_problem_compute_q
#include "plasticity_2d_solvers.h"
#include "plasticity_2d_vonmises_projection.h"  // for Von Mises projection
#include "numerics_verbose.h"

/* Solver registration system */
#include "solver_registry.h"
#include "numerics_errors.h"
#include "projectionOnCone.h"      // for projectionOnCone
#include "projectionOnCylinder.h"  // for projectionOnCylinder
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"  // for DEBUG_PRINTF, DEBUG_EXPR, DEBU...
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

/* Static variables */

/* The global problem of size n= 3*nc, nc being the number of contacts, is locally saved in
 * MGlobal and qGlobal */
/* note that either MGlobal or MBGlobal is used, depending on the chosen storage */
/* static int n=0; */
/* static const NumericsMatrix* MGlobal = NULL; */
/* static const double* qGlobal = NULL; */
/* static const double* theta = NULL; */

/* Local problem operators */
/* static const int nLocal = 3; */
/* static double* MLocal; */
/* static int isMAllocatedIn = 0; /\* True if a malloc is done for MLocal, else false *\/ */
/* static double qLocal[3]; */
/* static double theta_i = 0.0; */

void plasticity_2d_projection_initialize(PlasticityProblem* problem,
                                PlasticityProblem* localproblem) {}

void plasticity_2d_projection_update(int contact, PlasticityProblem* problem,
                            PlasticityProblem* localproblem, double* stress,
                            SolverOptions* options) {
  /* Build a local problem for a specific contact
     stress corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  plasticity_2d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.stressBlock, excluding the block corresponding to the current contact. ****/
  plasticity_2d_local_problem_compute_q(problem, localproblem, stress, contact);

  /* coefficient for current block - handle both model types */
  if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
    localproblem->model.drucker_prager->eta[0] = problem->model.drucker_prager->eta[contact];
    localproblem->model.drucker_prager->theta[0] = problem->model.drucker_prager->theta[contact];
  } else if (problem->model_type == PLASTICITY_MODEL_VON_MISES) {
    localproblem->model.von_mises->sigma_y[0] = problem->model.von_mises->sigma_y[contact];
  }
}

void plasticity_2d_projectionOnConeWithLocalIteration_initialize(PlasticityProblem* problem,
                                                        PlasticityProblem* localproblem,
                                                        SolverOptions* localsolver_options) {
  size_t nc = problem->numberOfCones;
  /* printf("plasticity_2d_projectionOnConeWithLocalIteration_initialize. Allocation of dwork\n"); */
  if (!localsolver_options->dWork || localsolver_options->dWorkSize < nc) {
    localsolver_options->dWork =
        (double*)realloc(localsolver_options->dWork, nc * sizeof(double));
    localsolver_options->dWorkSize = nc;
  }
  for (size_t i = 0; i < nc; i++) {
    localsolver_options->dWork[i] = 1.0;
  }
}
void plasticity_2d_projectionOnConeWithLocalIteration_free(PlasticityProblem* problem,
                                                  PlasticityProblem* localproblem,
                                                  SolverOptions* localsolver_options) {
  free(localsolver_options->dWork);
  localsolver_options->dWork = NULL;
}

int plasticity_2d_projectionOnConeWithLocalIteration_solve(PlasticityProblem* localproblem,
                                                  double* stress_local, SolverOptions* options) {
  DEBUG_BEGIN("plasticity_2d_projectionOnConeWithLocalIteration_solve(...)\n");

  DEBUG_EXPR(plasticity2DProblem_display(localproblem););

  double* MLocal = localproblem->M->matrix0;
  double* qLocal = localproblem->q;
  double eta_i = localproblem->model.drucker_prager->eta[0];
  double theta_i = localproblem->model.drucker_prager->theta[0];
  /* int nLocal = 3; */

  /*   /\* Builds local problem for the current contact *\/ */
  /*   plasticity_2d_projection_update(localproblem, stress_local); */

  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] +
   * MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  /* double an = 1. / (MLocal[0]); */

  /* double at = 1.0 / (MLocal[4] + theta_i); */
  /* double as = 1.0 / (MLocal[8] + theta_i); */
  /* at = an; */
  /* as = an; */
  double rho = options->dWork[options->iparam[PLASTICITY_CURRENT_CONE_NUMBER]], rho_k;
  DEBUG_PRINTF(" Contact options->iparam[PLASTICITY_CURRENT_CONTACT_NUMBER] = %i\n",
               options->iparam[PLASTICITY_CURRENT_CONTACT_NUMBER]);
  DEBUG_PRINTF("saved rho = %14.7e\n", rho);
  assert(rho > 0);

  /* int incx = 1, incy = 1; */
  int i;

  double strainrate_local[3], strainrate_k[3], stress_k[3], worktmp[3];
  double normUT;
  double localerror = 1.0;
  // printf ("localerror = %14.7e\n",localerror );
  int localiter = 0;
  double localtolerance = SOLVER_TOL(options);

  /* Variable for Line_search */
  double a1, a2;
  int success = 0;
  double localerror_k;
  int ls_iter = 0;
  int ls_itermax = 10;

  double tau = 2.0 / 3.0, tauinv = 3.0 / 2.0, L = 0.9, Lmin = 0.3;

  numerics_printf_verbose(2, "--  plasticity_2d_projectionOnConeWithLocalIteration_solve contact = %i",
                          options->iparam[PLASTICITY_CURRENT_CONE_NUMBER]);
  numerics_printf_verbose(2,
                          "--  plasticity_2d_projectionOnConeWithLocalIteration_solve | localiter \t| "
                          "rho \t\t\t| error\t\t\t|");
  numerics_printf_verbose(
      2, "--                                                | %i \t\t| %.10e\t| %.10e\t|",
      localiter, rho, localerror);

  /*     printf ("localtolerance = %14.7e\n",localtolerance ); */
  while ((localerror > localtolerance) && (localiter < SOLVER_MAX_ITER(options))) {
    DEBUG_PRINT("\n Local iteration starts \n");
    localiter++;

    /*    printf ("stress_local[0] = %14.7e\n",stress_local[0]); */
    /*    printf ("stress_local[1] = %14.7e\n",stress_local[1]); */
    /*    printf ("stress_local[2] = %14.7e\n",stress_local[2]); */

    /* Store the error */
    localerror_k = localerror;

    /* store the stress at the beginning of the iteration */
    /* cblas_dcopy(nLocal , stress_local , 1 , stress_k, 1); */

    stress_k[0] = stress_local[0];
    stress_k[1] = stress_local[1];
    stress_k[2] = stress_local[2];
    DEBUG_EXPR(NV_display(stress_k, 3););
    /* /\* strainrate_k <- q  *\/ */
    /* cblas_dcopy_msan(nLocal , qLocal , 1 , strainrate_k, 1); */
    /* /\* strainrate_k <- q + M * stress_local  *\/ */
    /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, stress_local,
     * incx, 1.0, strainrate_k, incy); */
    for (i = 0; i < 3; i++)
      strainrate_k[i] = MLocal[i + 0 * 3] * stress_local[0] + qLocal[i] +
                      MLocal[i + 1 * 3] * stress_local[1] + +MLocal[i + 2 * 3] * stress_local[2];
    DEBUG_EXPR(NV_display(strainrate_k, 3););
    ls_iter = 0;
    success = 0;
    rho_k = rho / tau;

    normUT = sqrt(strainrate_k[1] * strainrate_k[1] + strainrate_k[2] * strainrate_k[2]);
    while (!success && (ls_iter < ls_itermax)) {
      rho_k = rho_k * tau;
      DEBUG_PRINTF("rho_k =%f\n", rho_k);
      stress_local[0] = stress_k[0] - rho_k * (strainrate_k[0] + theta_i * normUT);
      stress_local[1] = stress_k[1] - rho_k * strainrate_k[1];
      stress_local[2] = stress_k[2] - rho_k * strainrate_k[2];
      DEBUG_PRINT("r-rho tilde v before projection")
      DEBUG_EXPR(NV_display(stress_local, 3););

      projectionOnCone(&stress_local[0], eta_i);

      /* strainrate <- q  */
      /* cblas_dcopy(nLocal , qLocal , 1 , strainrate, 1); */
      /* strainrate <- q + M * stress_local  */
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, stress_local,
       * incx, 1.0, strainrate, incy); */

      for (i = 0; i < 3; i++)
        strainrate_local[i] = MLocal[i + 0 * 3] * stress_local[0] + qLocal[i] +
                      MLocal[i + 1 * 3] * stress_local[1] + +MLocal[i + 2 * 3] * stress_local[2];

      a1 = sqrt((strainrate_k[0] - strainrate_local[0]) * (strainrate_k[0] - strainrate_local[0]) +
                (strainrate_k[1] - strainrate_local[1]) * (strainrate_k[1] - strainrate_local[1]) +
                (strainrate_k[2] - strainrate_local[2]) * (strainrate_k[2] - strainrate_local[2]));

      a2 = sqrt((stress_k[0] - stress_local[0]) * (stress_k[0] - stress_local[0]) +
                (stress_k[1] - stress_local[1]) * (stress_k[1] - stress_local[1]) +
                (stress_k[2] - stress_local[2]) * (stress_k[2] - stress_local[2]));

      success = (rho_k * a1 <= L * a2) ? 1 : 0;

      DEBUG_PRINTF("rho_k = %12.8e\t", rho_k);
      DEBUG_PRINTF("a1 = %12.8e\t", a1);
      DEBUG_PRINTF("a2 = %12.8e\t", a2);
      DEBUG_PRINTF("norm stress = %12.8e\t",
                   sqrt(stress_local[0] * stress_local[0] + stress_local[1] * stress_local[1] +
                        stress_local[2] * stress_local[2]));
      DEBUG_PRINTF("success = %i\n", success);

      ls_iter++;
    }

    /* printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho,
     * localerror); */

    /* compute local error */
    localerror = 0.0;
    plasticity_2d_unitary_compute_and_add_error(stress_local, strainrate_local, eta_i, theta_i, &localerror, worktmp);

    /*Update rho*/
    if ((rho_k * a1 < Lmin * a2) && (localerror < localerror_k)) {
      rho = rho_k * tauinv;
    } else
      rho = rho_k;

    numerics_printf_verbose(
        2, "--                                                | %i \t\t| %.10e\t| %.10e\t|",
        localiter, rho, localerror);
  }
  options->dWork[options->iparam[PLASTICITY_CURRENT_CONE_NUMBER]] = rho;
  SET_SOLVER_RESIDUAL(options, localerror);
  DEBUG_PRINTF("final rho  =%e\n", rho);

  DEBUG_END("plasticity_2d_projectionOnConeWithLocalIteration_solve(...)\n");
  if (localerror > localtolerance) return 1;
  return 0;
}

int plasticity_2d_projectionOnCone_solve(PlasticityProblem* localproblem, double* stress_local,
                                SolverOptions* options) {
  /*  /\* Builds local problem for the current contact *\/ */
  /*   plasticity_2d_projection_update(contact, stress_local); */

  double* MLocal = localproblem->M->matrix0;
  double* qLocal = localproblem->q;
  double eta_i = localproblem->model.drucker_prager->eta[0];
  double theta_i = localproblem->model.drucker_prager->theta[0];
  /* int nLocal = 3; */

  /* this part is critical for the success of the projection */
  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] +
   * MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  // double an = 1./(MLocal[0]+theta_i);
  double an = 1. / (MLocal[0]);

  /* int incx = 1, incy = 1; */
  double worktmp[3];
  double normUT;
  /* cblas_dcopy_msan(nLocal , qLocal, incx , worktmp , incy); */
  /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, stress_local,
   * incx, 1.0, worktmp, incy); */

  for (int i = 0; i < 3; i++)
    worktmp[i] = MLocal[i + 0 * 3] * stress_local[0] + qLocal[i] +
                 MLocal[i + 1 * 3] * stress_local[1] + +MLocal[i + 2 * 3] * stress_local[2];

  normUT = sqrt(worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
  stress_local[0] -= an * (worktmp[0] + theta_i * normUT);
  stress_local[1] -= an * worktmp[1];
  stress_local[2] -= an * worktmp[2];

  projectionOnCone(stress_local, eta_i);
  return 0;
}

void plasticity_2d_projection_free(PlasticityProblem* problem, PlasticityProblem* localproblem,
                          SolverOptions* localsolver_options) {}



void plasticity_2d_poc_set_default(SolverOptions* options) {
  options->iparam[PLASTICITY_CURRENT_CONE_NUMBER] = 0;  // this will be set by external solver
  options->dparam[PLASTICITY_NSN_RHO] =
      0.;  // Used only for ProjectionOnConeWithRegularization
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers local one-cone projection solvers in the global solver registry
 */

/* Wrapper for projection on cone (local solver) - Note: local solvers use 3 args */
static void plasticity_2d_projectionOnCone_set_default(SolverOptions* options) {
  /* No specific defaults */
}

static int plasticity_2d_projectionOnCone_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_projectionOnCone_solve_wrap(void* localproblem, double* reaction,
                                            double* velocity, SolverOptions* options) {
  (void)velocity;  /* Local solvers don't use velocity parameter */
  return plasticity_2d_projectionOnCone_solve((PlasticityProblem*)localproblem, reaction, options);
}

/* Wrapper for projection on cone with local iteration (local solver) */
static void plasticity_2d_projectionOnConeWithLocalIteration_set_default(SolverOptions* options) {
  /* No specific defaults */
}

static int plasticity_2d_projectionOnConeWithLocalIteration_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_projectionOnConeWithLocalIteration_solve_wrap(void* localproblem,
                                                               double* reaction,
                                                               double* velocity,
                                                               SolverOptions* options) {
  (void)velocity;  /* Local solvers don't use velocity parameter */
  return plasticity_2d_projectionOnConeWithLocalIteration_solve((PlasticityProblem*)localproblem,
                                                        reaction, options);
}

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_ProjectionOnCone,
                "PLASTICITY_2D_ONECONE_PROJECTION",
                "Projection on cone for 2D Mohr Coulomb (one cone)",
                plasticity_2d_projectionOnCone_init_wrap,
                plasticity_2d_projectionOnCone_solve_wrap,
                NULL,
                NULL,
                plasticity_2d_projectionOnCone_set_default,
                100,   /* default_max_iter */
                1e-14, /* default_tol */
                1);    /* is_local */

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration,
                "PLASTICITY_2D_ONECONE_PROJECTION_LI",
                "Projection on cone with local iteration for 2D Mohr Coulomb (one cone)",
                plasticity_2d_projectionOnConeWithLocalIteration_init_wrap,
                plasticity_2d_projectionOnConeWithLocalIteration_solve_wrap,
                NULL,
                NULL,
                plasticity_2d_projectionOnConeWithLocalIteration_set_default,
                100,   /* default_max_iter */
                1e-14, /* default_tol */
                1);    /* is_local */

/* ===========================================================================
 * Von Mises Projection Solver
 * ===========================================================================
 * Simple radial return algorithm for Von Mises plasticity
 */

int plasticity_2d_projectionOnVonMises_solve(PlasticityProblem* localproblem, double* stress_local,
                                 SolverOptions* options) {
  (void)options;  /* Unused for now - simple projection doesn't need options */
  
  double* MLocal = localproblem->M->matrix0;
  double* qLocal = localproblem->q;
  double sigma_y = localproblem->model.von_mises->sigma_y[0];
  
  /* Compute trial stress: sigma_trial = MLocal * stress_local + qLocal */
  double trial_stress[3];
  for (int i = 0; i < 3; i++) {
    trial_stress[i] = MLocal[i + 0 * 3] * stress_local[0] + qLocal[i] +
                      MLocal[i + 1 * 3] * stress_local[1] +
                      MLocal[i + 2 * 3] * stress_local[2];
  }
  
  /* Project onto Von Mises yield surface */
  plasticity_2d_projectionOnVonMises(trial_stress, sigma_y);
  
  /* Update stress_local */
  stress_local[0] = trial_stress[0];
  stress_local[1] = trial_stress[1];
  stress_local[2] = trial_stress[2];
  
  return 0;
}

void plasticity_2d_vonmises_set_default(SolverOptions* options) {
  (void)options;
}

/* Wrapper for Von Mises projection (local solver) */
static int plasticity_2d_projectionOnVonMises_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_projectionOnVonMises_solve_wrap(void* localproblem,
                                                         double* reaction,
                                                         double* velocity,
                                                         SolverOptions* options) {
  (void)velocity;  /* Local solvers don't use velocity parameter */
  return plasticity_2d_projectionOnVonMises_solve((PlasticityProblem*)localproblem,
                                                  reaction, options);
}

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_VONMISES,
                "PLASTICITY_2D_ONECONE_VONMISES",
                "Radial return projection for Von Mises plasticity (one cone)",
                plasticity_2d_projectionOnVonMises_init_wrap,
                plasticity_2d_projectionOnVonMises_solve_wrap,
                NULL,
                NULL,
                plasticity_2d_vonmises_set_default,
                100,   /* default_max_iter */
                1e-14, /* default_tol */
                1);    /* is_local */
