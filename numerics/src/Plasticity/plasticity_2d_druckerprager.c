// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "plasticity_2d_druckerprager.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt
#include <stdio.h>   // for fprintf, printf, NULL, stderr
#include <stdlib.h>  // for calloc, free, exit, EXIT_FAILURE

#include "NumericsFwd.h"         // for SolverOptions, MohrCoulomb2D...
#include "NumericsMatrix.h"      // for NumericsMatrix, RawNumericsMatrix
#include "PlasticityProblem.h"   // for PlasticityProblem
#include "Plasticity_options.h"  // for PLASTICITY_NSGS_LOCAL...
#include "naming_conventions.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "op3x3.h"                        // for solve_3x3_gepp
#include "plasticity_2d_compute_error.h"  // for plasticity_2d_Tresca_unitary_compute_an...
#include "plasticity_2d_local_newton.h"
#include "plasticity_2d_local_problem_tools.h"  // for plasticity_2d_local_problem_compute_q
#include "plasticity_2d_solvers.h"
#include "projection_soc.h"  // for soc_projection_3d
#include "quartic.h"         // for BIQUADROOTS
#include "solver_registry.h"
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
    localproblem->model.drucker_prager->theta[0] =
        problem->model.drucker_prager->theta[contact];
  } else if (problem->model_type == PLASTICITY_MODEL_VON_MISES) {
    localproblem->model.von_mises->sigma_y[0] = problem->model.von_mises->sigma_y[contact];
  }
}

void plasticity_2d_projectionOnConeWithLocalIteration_initialize(
    PlasticityProblem* problem, PlasticityProblem* localproblem,
    SolverOptions* localsolver_options) {
  size_t nc = problem->numberOfCones;
  /* printf("plasticity_2d_projectionOnConeWithLocalIteration_initialize. Allocation of
   * dwork\n"); */
  if (!localsolver_options->dWork || localsolver_options->dWorkSize < nc) {
    localsolver_options->dWork =
        (double*)realloc(localsolver_options->dWork, nc * sizeof(double));
    localsolver_options->dWorkSize = nc;
  }
  for (size_t i = 0; i < nc; i++) {
    localsolver_options->dWork[i] = 1.0;
  }
}
void plasticity_2d_projectionOnConeWithLocalIteration_free(
    PlasticityProblem* problem, PlasticityProblem* localproblem,
    SolverOptions* localsolver_options) {
  free(localsolver_options->dWork);
  localsolver_options->dWork = NULL;
}

int plasticity_2d_projectionOnConeWithLocalIteration_solve(PlasticityProblem* localproblem,
                                                           double* stress_local,
                                                           SolverOptions* options) {
  DEBUG_BEGIN("plasticity_2d_projectionOnConeWithLocalIteration_solve(...)\n");

  DEBUG_EXPR(plasticity2DProblem_display(localproblem););

  double* MLocal = localproblem->M->dense;
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

  double strainrate_local[3], strainrate_k[3], stress_k[3];
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

  numerics_printf_verbose(
      2, "--  plasticity_2d_projectionOnConeWithLocalIteration_solve contact = %i",
      options->iparam[PLASTICITY_CURRENT_CONE_NUMBER]);
  numerics_printf_verbose(
      2,
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
                        MLocal[i + 1 * 3] * stress_local[1] +
                        +MLocal[i + 2 * 3] * stress_local[2];
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

      soc_projection_3d(&stress_local[0], eta_i);

      /* strainrate <- q  */
      /* cblas_dcopy(nLocal , qLocal , 1 , strainrate, 1); */
      /* strainrate <- q + M * stress_local  */
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, stress_local,
       * incx, 1.0, strainrate, incy); */

      for (i = 0; i < 3; i++)
        strainrate_local[i] = MLocal[i + 0 * 3] * stress_local[0] + qLocal[i] +
                              MLocal[i + 1 * 3] * stress_local[1] +
                              +MLocal[i + 2 * 3] * stress_local[2];

      a1 = sqrt(
          (strainrate_k[0] - strainrate_local[0]) * (strainrate_k[0] - strainrate_local[0]) +
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
    /* scaled as the global error (it was the squared absolute residual) */
    localerror = plasticity_2d_DruckerPrager_local_scaled_error(MLocal, qLocal, eta_i, theta_i,
                                                                stress_local);

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

  double* MLocal = localproblem->M->dense;
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

  soc_projection_3d(stress_local, eta_i);
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
static void plasticity_2d_projection_on_cone_set_default(SolverOptions* options) {
  /* No specific defaults */
}

static int plasticity_2d_projection_on_cone_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_projection_on_cone_solve_wrap(void* localproblem, double* reaction,
                                                       double* velocity,
                                                       SolverOptions* options) {
  (void)velocity; /* Local solvers don't use velocity parameter */
  return plasticity_2d_projectionOnCone_solve((PlasticityProblem*)localproblem, reaction,
                                              options);
}

/* Wrapper for projection on cone with local iteration (local solver) */
static void plasticity_2d_projection_on_cone_with_local_iteration_set_default(
    SolverOptions* options) {
  /* No specific defaults */
}

static int plasticity_2d_projection_on_cone_with_local_iteration_init_wrap(
    void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_projection_on_cone_with_local_iteration_solve_wrap(
    void* localproblem, double* reaction, double* velocity, SolverOptions* options) {
  (void)velocity; /* Local solvers don't use velocity parameter */
  return plasticity_2d_projectionOnConeWithLocalIteration_solve(
      (PlasticityProblem*)localproblem, reaction, options);
}

/* ===========================================================================
 * Drucker-Prager, one cone: natural map, return mapping, semi-smooth Newton
 *   sigma in K, -(w + theta ||w_T|| e_N) in N_K(sigma), K = {||sigma_T|| <= eta sigma_N}
 * =========================================================================== */

static double dp_norm3(const double a[3]) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/* y = M x */
static void dp_mv(const double M[9], const double x[3], double y[3]) {
  for (int i = 0; i < 3; i++) y[i] = M[i] * x[0] + M[i + 3] * x[1] + M[i + 6] * x[2];
}

/* w = M s + q */
static void dp_w(const double M[9], const double q[3], const double s[3], double w[3]) {
  for (int i = 0; i < 3; i++) w[i] = q[i] + M[i] * s[0] + M[i + 3] * s[1] + M[i + 6] * s[2];
}

/* generalized Jacobian (column-major) of the projection on K at y, same cases as
 * soc_projection_3d: polar cone (0), inside (I), boundary */
void plasticity_2d_DruckerPrager_projection_jacobian(const double y[3], double eta,
                                                     double J[9]) {
  soc_projection_jacobian_3d(J, y, eta);
}

/* Phi = sigma - P_K(sigma - rho v), v = w + theta ||w_T|| e_N; returns ||Phi||, and the
 * generalized Jacobian H = I - J_P(y)(I - rho J_v) if H != NULL, with
 * J_v = (I + theta e_N [0, n^T]) M, n = w_T/||w_T|| (0 if w_T = 0) */
static double dp_natural_map(const double M[9], const double q[3], double eta, double theta,
                             double rho, const double s[3], double Phi[3], double H[9]) {
  double w[3];
  dp_w(M, q, s, w);
  double nwT = sqrt(w[1] * w[1] + w[2] * w[2]);
  double y[3] = {s[0] - rho * (w[0] + theta * nwT), s[1] - rho * w[1], s[2] - rho * w[2]};
  double P[3] = {y[0], y[1], y[2]};
  soc_projection_3d(P, eta);
  for (int i = 0; i < 3; i++) Phi[i] = s[i] - P[i];
  if (H) {
    double Jv[9], J[9];
    for (int i = 0; i < 9; i++) Jv[i] = M[i];
    if (nwT > 0.0) {
      double n1 = w[1] / nwT, n2 = w[2] / nwT;
      for (int j = 0; j < 3; j++) Jv[3 * j] += theta * (n1 * M[1 + 3 * j] + n2 * M[2 + 3 * j]);
    }
    plasticity_2d_DruckerPrager_projection_jacobian(y, eta, J);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double sum = 0.0;
        for (int k = 0; k < 3; k++)
          sum += J[i + 3 * k] * ((k == j ? 1.0 : 0.0) - rho * Jv[k + 3 * j]);
        H[i + 3 * j] = (i == j ? 1.0 : 0.0) - sum;
      }
  }
  return dp_norm3(Phi);
}

double plasticity_2d_DruckerPrager_natural_map_residual(const double M[9], const double q[3],
                                                        double eta, double theta, double rho,
                                                        const double stress[3]) {
  double Phi[3];
  return dp_natural_map(M, q, eta, theta, rho, stress, Phi, NULL);
}

static double dp_det3(const double a[3], const double b[3], const double c[3]) {
  return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
         a[2] * (b[0] * c[1] - b[1] * c[0]);
}

/* sliding: sigma = sigma_N a(phi), w = lambda b(phi), a = (1, eta c, eta s),
 * b = (eta - theta, -c, -s), c = cos(phi), s = sin(phi). F = sigma_N M a - lambda b + q = 0.
 * From an angle phi: least-squares (sigma_N, lambda), then Newton on (sigma_N, lambda, phi).
 * Returns the stress. */
static void dp_sliding_candidate(const double M[9], const double q[3], double eta,
                                 double theta, double phi, double stress[3]) {
  double z[3] = {0.0, 0.0, phi}; /* sigma_N, lambda, phi */
  double scale = dp_norm3(q) + DBL_MIN;
  for (int it = 0; it < 30; it++) {
    double c = cos(z[2]), s = sin(z[2]);
    double a[3] = {1.0, eta * c, eta * s}, b[3] = {eta - theta, -c, -s};
    double Ma[3];
    dp_mv(M, a, Ma);
    if (it == 0) {
      /* least squares for (sigma_N, lambda): [Ma, -b] z = -q */
      double g11 = Ma[0] * Ma[0] + Ma[1] * Ma[1] + Ma[2] * Ma[2];
      double g12 = -(Ma[0] * b[0] + Ma[1] * b[1] + Ma[2] * b[2]);
      double g22 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
      double r1 = -(Ma[0] * q[0] + Ma[1] * q[1] + Ma[2] * q[2]);
      double r2 = b[0] * q[0] + b[1] * q[1] + b[2] * q[2];
      double det = g11 * g22 - g12 * g12;
      if (!(fabs(det) > DBL_MIN)) break;
      z[0] = (g22 * r1 - g12 * r2) / det;
      z[1] = (g11 * r2 - g12 * r1) / det;
    }
    double F[3], da[3] = {0.0, -eta * s, eta * c}, db[3] = {0.0, s, -c}, Mda[3];
    dp_mv(M, da, Mda);
    for (int i = 0; i < 3; i++) F[i] = z[0] * Ma[i] - z[1] * b[i] + q[i];
    if (dp_norm3(F) <= 1e-15 * scale) break;
    double J[9];
    for (int i = 0; i < 3; i++) {
      J[i] = Ma[i];
      J[i + 3] = -b[i];
      J[i + 6] = z[0] * Mda[i] - z[1] * db[i];
    }
    double dz[3] = {-F[0], -F[1], -F[2]};
    if (solve_3x3_gepp(J, dz)) break;
    for (int i = 0; i < 3; i++) z[i] += dz[i];
  }
  stress[0] = z[0];
  stress[1] = z[0] * eta * cos(z[2]);
  stress[2] = z[0] * eta * sin(z[2]);
}

/* angles phi of the possible sliding solutions: det[M a(phi), b(phi), q] = 0. With
 * u = (1, cos, sin), det = sum_ij D_ij u_i u_j, D_ij = det[A_i, B_j, q], A = (M e_0,
 * eta M e_1, eta M e_2), B = ((eta - theta) e_0, -e_1, -e_2); tau = tan(phi/2) gives a
 * quartic: (1 + tau^2)^2 det = sum_ij D_ij U_i U_j, U = (1 + tau^2, 1 - tau^2, 2 tau). */
static int dp_sliding_angles(const double M[9], const double q[3], double eta, double theta,
                             double phis[8]) {
  const double A[3][3] = {{M[0], M[1], M[2]},
                          {eta * M[3], eta * M[4], eta * M[5]},
                          {eta * M[6], eta * M[7], eta * M[8]}};
  const double B[3][3] = {{eta - theta, 0., 0.}, {0., -1., 0.}, {0., 0., -1.}};
  const double U[3][3] = {{1., 0., 1.}, {1., 0., -1.}, {0., 2., 0.}}; /* powers 0..2 */
  double P[5] = {0., 0., 0., 0., 0.};                                 /* power k */
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double D = dp_det3(A[i], B[j], q);
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) P[a + b] += D * U[i][a] * U[j][b];
    }
  int n = 0;
  phis[n++] = M_PI; /* tau = infinity */
  double pmax = 0.0;
  for (int k = 0; k < 5; k++) pmax = fmax(pmax, fabs(P[k]));
  if (pmax == 0.0) {
    /* degenerate: try a few directions */
    for (int k = 0; k < 6; k++) phis[n++] = 2.0 * M_PI * k / 7.0;
    return n;
  }
  int deg = 4;
  while (deg > 0 && fabs(P[deg]) <= 1e-13 * pmax) deg--;
  double p[5], r[3][5];
  for (int k = 0; k <= deg; k++) p[k] = P[deg - k]; /* leading coefficient first */
  int nroots = 0;
  if (deg == 4) {
    BIQUADROOTS(p, r);
    nroots = 4;
  } else if (deg == 3) {
    CUBICROOTS(p, r);
    nroots = 3;
  } else if (deg == 2) {
    QUADROOTS(p, r);
    nroots = 2;
  } else if (deg == 1) {
    r[1][1] = -p[1] / p[0];
    r[2][1] = 0.0;
    nroots = 1;
  }
  for (int k = 1; k <= nroots; k++)
    if (fabs(r[2][k]) <= 1e-6 * (1.0 + fabs(r[1][k]))) phis[n++] = 2.0 * atan(r[1][k]);
  return n;
}

int plasticity_2d_DruckerPrager_return_mapping(const double M[9], const double q[3],
                                               double eta, double theta, double stress[3]) {
  double n2 = plasticity_2d_local_spectral_norm(M);
  double rho0 = (n2 > 0.0) ? 1.0 / n2 : 1.0;
  double candidates[10][3];
  int nc = 0;
  /* elastic: w = 0 */
  candidates[nc][0] = -q[0];
  candidates[nc][1] = -q[1];
  candidates[nc][2] = -q[2];
  if (solve_3x3_gepp(M, candidates[nc]) == 0) nc++;
  /* apex */
  candidates[nc][0] = candidates[nc][1] = candidates[nc][2] = 0.0;
  nc++;
  /* sliding */
  double phis[8];
  int nphi = dp_sliding_angles(M, q, eta, theta, phis);
  for (int k = 0; k < nphi; k++)
    dp_sliding_candidate(M, q, eta, theta, phis[k], candidates[nc++]);

  /* solutions: zero natural-map residual (up to round-off); the closest to the input stress
   * is chosen (continuity in the Gauss-Seidel iterations when the solution is not unique,
   * theta != eta or non-symmetric M) */
  int best = -1, best_valid = -1;
  double best_res = INFINITY, best_dist = INFINITY;
  for (int k = 0; k < nc; k++) {
    double Phi[3];
    double res = dp_natural_map(M, q, eta, theta, rho0, candidates[k], Phi, NULL);
    double scale = 1.0 + dp_norm3(candidates[k]) + rho0 * dp_norm3(q);
    if (res < best_res) {
      best_res = res;
      best = k;
    }
    if (res <= 1e-10 * scale) {
      double d[3] = {candidates[k][0] - stress[0], candidates[k][1] - stress[1],
                     candidates[k][2] - stress[2]};
      double dist = dp_norm3(d);
      if (dist < best_dist) {
        best_dist = dist;
        best_valid = k;
      }
    }
  }
  int chosen = (best_valid >= 0) ? best_valid : best;
  for (int i = 0; i < 3; i++) stress[i] = candidates[chosen][i];
  return (best_valid >= 0) ? 0 : 1;
}

/* semi-smooth Newton (shared implementation, plasticity_2d_local_newton.h) */
typedef struct {
  const double* M;
  const double* q;
  double eta;
  double theta;
} DruckerPragerData;

static double dp_newton_map(void* data, double rho, const double s[3], double Phi[3],
                            double H[9]) {
  DruckerPragerData* d = (DruckerPragerData*)data;
  return dp_natural_map(d->M, d->q, d->eta, d->theta, rho, s, Phi, H);
}

static int dp_newton_exact(void* data, double s[3]) {
  DruckerPragerData* d = (DruckerPragerData*)data;
  return plasticity_2d_DruckerPrager_return_mapping(d->M, d->q, d->eta, d->theta, s);
}

int plasticity_2d_DruckerPrager_semismooth_newton(const double M[9], const double q[3],
                                                  double eta, double theta, double stress[3],
                                                  Plasticity2DLocalParameters* p) {
  DruckerPragerData data = {M, q, eta, theta};
  return plasticity_2d_local_semismooth_newton(dp_newton_map, dp_newton_exact, &data, M,
                                               stress, p);
}

/* SolverOptions interface */
int plasticity_2d_DruckerPrager_return_mapping_solve(PlasticityProblem* localproblem,
                                                     double* stress_local,
                                                     SolverOptions* options) {
  int status = plasticity_2d_DruckerPrager_return_mapping(
      localproblem->M->dense, localproblem->q, localproblem->model.drucker_prager->eta[0],
      localproblem->model.drucker_prager->theta[0], stress_local);
  if (options) SET_SOLVER_ITER_DONE(options, 1);
  return status;
}

int plasticity_2d_DruckerPrager_nsn_solve(PlasticityProblem* localproblem,
                                          double* stress_local, SolverOptions* options) {
  Plasticity2DLocalParameters p = {
      .max_iter = SOLVER_MAX_ITER(options),
      .tol = SOLVER_TOL(options),
      .rho_strategy = options->iparam[PLASTICITY_NSN_RHO_STRATEGY],
      .rho = options->dparam[PLASTICITY_NSN_RHO],
      .linesearch = options->iparam[PLASTICITY_NSN_LINESEARCH],
      .linesearch_max_iter = options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER],
      .hybrid = (options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] !=
                 PLASTICITY_NSN_HYBRID_STRATEGY_NO)};
  int info = plasticity_2d_DruckerPrager_semismooth_newton(
      localproblem->M->dense, localproblem->q, localproblem->model.drucker_prager->eta[0],
      localproblem->model.drucker_prager->theta[0], stress_local, &p);
  SET_SOLVER_ITER_DONE(options, p.iter_done);
  SET_SOLVER_RESIDUAL(options, p.residual);
  return info;
}

static void plasticity_2d_druckerprager_return_mapping_set_default(SolverOptions* options) {
  (void)options;
}

void plasticity_2d_druckerprager_nsn_set_default(SolverOptions* options) {
  options->iparam[PLASTICITY_NSN_RHO_STRATEGY] =
      PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM;
  options->iparam[PLASTICITY_NSN_LINESEARCH] = PLASTICITY_NSN_LINESEARCH_ARMIJO;
  options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER] = 20;
  /* fall back to the return mapping if the Newton method fails */
  options->iparam[PLASTICITY_NSN_HYBRID_STRATEGY] =
      PLASTICITY_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP;
  options->dparam[PLASTICITY_NSN_RHO] = 0.0;
}

static int plasticity_2d_druckerprager_local_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_druckerprager_return_mapping_wrap(void* localproblem,
                                                           double* reaction, double* velocity,
                                                           SolverOptions* options) {
  (void)velocity;
  return plasticity_2d_DruckerPrager_return_mapping_solve((PlasticityProblem*)localproblem,
                                                          reaction, options);
}

static int plasticity_2d_druckerprager_nsn_wrap(void* localproblem, double* reaction,
                                                double* velocity, SolverOptions* options) {
  (void)velocity;
  return plasticity_2d_DruckerPrager_nsn_solve((PlasticityProblem*)localproblem, reaction,
                                               options);
}

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING,
                "PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING",
                "Drucker-Prager, one cone: return mapping (exact enumeration: elastic, apex, "
                "sliding)",
                plasticity_2d_druckerprager_local_init_wrap,
                plasticity_2d_druckerprager_return_mapping_wrap, NULL, NULL,
                plasticity_2d_druckerprager_return_mapping_set_default, 1, /* max_iter */
                1e-14,                                                     /* tol */
                1);                                                        /* is_local */

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN,
                "PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN",
                "Drucker-Prager, one cone: semi-smooth Newton on the natural map",
                plasticity_2d_druckerprager_local_init_wrap,
                plasticity_2d_druckerprager_nsn_wrap, NULL, NULL,
                plasticity_2d_druckerprager_nsn_set_default, 50, /* max_iter */
                1e-14,                                           /* tol */
                1);                                              /* is_local */

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_ProjectionOnCone, "PLASTICITY_2D_ONECONE_PROJECTION",
                "Projection on cone for 2D Mohr Coulomb (one cone)",
                plasticity_2d_projection_on_cone_init_wrap,
                plasticity_2d_projection_on_cone_solve_wrap, NULL, NULL,
                plasticity_2d_projection_on_cone_set_default, 100, /* default_max_iter */
                1e-14,                                             /* default_tol */
                1);                                                /* is_local */

REGISTER_SOLVER(PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration,
                "PLASTICITY_2D_ONECONE_PROJECTION_LI",
                "Projection on cone with local iteration for 2D Mohr Coulomb (one cone)",
                plasticity_2d_projection_on_cone_with_local_iteration_init_wrap,
                plasticity_2d_projection_on_cone_with_local_iteration_solve_wrap, NULL, NULL,
                plasticity_2d_projection_on_cone_with_local_iteration_set_default,
                100,   /* default_max_iter */
                1e-14, /* default_tol */
                1);    /* is_local */
