/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <openblas_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FrictionContact_options.h"
#include "Friction_tools.h"  // for ComputeErrorPtr
#include "JordanAlgebra.h"
#include "NumericsFwd.h"
#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include "SiconosLapack.h"  // IWYU pragma: keep
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "fc3d_short_names.h"  // for ComputeErrorPtr
#include "gfc3d_ipm.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "safe_casts.h"
#include "solver_registry.h"

/* ------------------------- Helper functions implementation ------------------------------ */
/* Compute the primal constraint vector for local fricprob: out = Wr + q + Es - u
   and the relative 2-norm of this vector: rnorm = |out|/max{|Wr|, |q|, |velocity|, |Es|} */
static void primalResidual_s_for_local_friction(const double* velocity, const double* reaction,
                                                NumericsMatrix* W, const double* q,
                                                const double* s, double* out, double* rnorm,
                                                const double tol) {
  const size_t nd = to_size_t(W->size0);
  double rn;
  const blasint blasnd = to_blasint(nd);

  NM_gemv(1.0, W, reaction, 0.0, out);
  rn = cblas_dnrm2(blasnd, out, 1);
  cblas_daxpy(blasnd, 1.0, q, 1, out, 1);
  cblas_daxpy(blasnd, -1.0, velocity, 1, out, 1);

  for (size_t i = 0; i < nd; i += 3) out[i] += s[i / 3];
  rn = fmax(rn, cblas_dnrm2(blasnd, velocity, 1));
  rn = fmax(rn, cblas_dnrm2(blasnd, q, 1));
  rn = fmax(rn, cblas_dnrm2(blasnd / 3, s, 1));
  *rnorm = (rn > tol ? cblas_dnrm2(blasnd, out, 1) : cblas_dnrm2(blasnd, out, 1));
}

/* --------------------------- Interior-point method implementation
 * ------------------------------ */
/*
 * Implementation contains the following functions:
 *  - fc3d_IPM_SNM_init - initialize solver (allocate memory)
 *  - fc3d_IPM_SNM_free - deallocate memory
 *  - fc3d_ipm_snm_set_default - setup default solver parameters
 *  - fc3d_IPM_SNM - optimization method
 *
 * NOTE: Friction contact (local) has the same structure of velocity and reaction with Global
 * friction contact. So, it is not necessary to create a separate struct
 */
void fc3d_IPM_SNM_init(FrictionContactProblem* problem, SolverOptions* options) {
  size_t n = to_size_t(problem->numberOfContacts);
  // size_t m = 3 * n;
  size_t d = to_size_t(problem->dimension);
  size_t nd = n * d;

  const size_t worksize = 2 * nd + n;
  if (!options->dWork || options->dWorkSize < worksize) {
    double* p = (double*)realloc(options->dWork, worksize * sizeof(double));
    if (!p) {
      numerics_error("fc3d_IP_SNM_init", "bad alloc");
    }
    options->dWork = p;
    options->dWorkSize = worksize;
  }

  /* ------------- initialize starting point ------------- */
  Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)calloc(1, sizeof(Gfc3d_IPM_init_data));
  options->solverData = data;

  /* --------- allocate memory for tmp point ----------- */
  data->tmp_point = (IPM_tmp_point*)malloc(sizeof(IPM_tmp_point));
  data->tmp_point->t_velocity = (double*)calloc(nd, sizeof(double));
  data->tmp_point->t_reaction = (double*)calloc(nd, sizeof(double));

  /* NOTE: Assigning values for starting points is done in fc3d_IPM_SNM for convenience */
  // /* u */
  // data->starting_point->velocity = (double*)calloc(nd, sizeof(double));
  // for(size_t i = 0; i < nd; ++ i)
  // {
  //   data->starting_point->velocity[i] = 0.01;
  //   if(i % d == 0)
  //     data->starting_point->velocity[i] = 0.1;
  // }

  // /* r */
  // data->starting_point->reaction = (double*)calloc(nd, sizeof(double));
  // for(size_t i = 0; i < nd; ++ i)
  // {
  //   data->starting_point->reaction[i] = 0.01; //0.0351;
  //   if(i % d == 0)
  //     data->starting_point->reaction[i] = 0.1; //0.2056;
  // }

  /* ------ initialize the change of variable matrix P_mu ------- */
  data->P_mu = (IPM_change_of_variable*)calloc(1, sizeof(IPM_change_of_variable));
  int ndint = to_int(nd);  // TEMP ...
  data->P_mu->mat = NM_create(NM_SPARSE, ndint, ndint);
  NM_triplet_alloc(data->P_mu->mat, ndint);
  for (size_t i = 0; i < nd; ++i)
    if (i % d == 0)
      NM_entry(data->P_mu->mat, i, i, 1.);
    else
      NM_entry(data->P_mu->mat, i, i, problem->mu[(size_t)(i / d)]);

  /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
  data->P_mu->inv_mat = NM_create(NM_SPARSE, ndint, ndint);
  NM_triplet_alloc(data->P_mu->inv_mat, ndint);
  for (size_t i = 0; i < nd; ++i)
    if (i % d == 0)
      NM_entry(data->P_mu->inv_mat, i, i, 1.);
    else
      NM_entry(data->P_mu->inv_mat, i, i, 1.0 / problem->mu[(int)(i / d)]);

  /* ------ initial parameters initialization ---------- */
  data->internal_params = (IPM_internal_params*)malloc(sizeof(IPM_internal_params));
  data->internal_params->alpha_primal = 1.0;
  data->internal_params->alpha_dual = 1.0;
  data->internal_params->sigma = 0.1;
  data->internal_params->barr_param = 1.0;

  /* ----- temporary vaults initialization ------- */
  data->tmp_vault_nd = (double**)malloc(10 * sizeof(double*));
  for (size_t i = 0; i < 10; ++i) data->tmp_vault_nd[i] = (double*)calloc(nd, sizeof(double));
}

static void fc3d_IPM_SNM_free(FrictionContactProblem* problem, Gfc3d_IPM_init_data* data) {
  if (!data) return;

  if (data->P_mu) {
    data->P_mu->mat = NM_free(data->P_mu->mat);
    data->P_mu->inv_mat = NM_free(data->P_mu->inv_mat);
    free(data->P_mu);
  }
  data->P_mu = NULL;

  if (data->tmp_vault_nd) {
    for (size_t i = 0; i < 10; ++i) free(data->tmp_vault_nd[i]);
    free(data->tmp_vault_nd);
  }
  data->tmp_vault_nd = NULL;

  if (data->tmp_point) {
    free(data->tmp_point->t_velocity);
    data->tmp_point->t_velocity = NULL;
    free(data->tmp_point->t_reaction);
    data->tmp_point->t_reaction = NULL;
    free(data->tmp_point);
  }
  data->tmp_point = NULL;

  if (data->internal_params) {
    free(data->internal_params);
  }
  data->internal_params = NULL;

  free(data);
}

void fc3d_IPM_SNM(FrictionContactProblem* restrict problem, double* restrict reaction,
                  double* restrict velocity, int* restrict info,
                  SolverOptions* restrict options) {
  // the size of the problem detection
  const size_t n = to_size_t(problem->numberOfContacts);
  const blasint blas_n = to_blasint(n);
  const size_t m = 3 * n;
  const size_t d = to_size_t(problem->dimension);
  const size_t nd = n * d;
  const blasint blasnd = to_blasint(nd);
  size_t no_nd = 0;
  const blasint size_blas = 2 * to_blasint(nd) + to_blasint(n);
  const size_t size_sol = 2 * nd + n;

  // initialize solver if it is not set
  int internal_allocation = 0;
  if (!options->dWork || (options->dWorkSize != size_sol)) {
    fc3d_IPM_SNM_init(problem, options);
    internal_allocation = 1;
  }

  Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)options->solverData;
  NumericsMatrix* P_mu = data->P_mu->mat;
  NumericsMatrix* P_mu_inv = data->P_mu->inv_mat;

  NumericsMatrix* W_ori = problem->M;
  double* q_ori = problem->q;
  double* q = data->tmp_vault_nd[no_nd++];

  // change of variable to eliminate the friction coefficients: W_ori --> W and q_ori --> q
  NumericsMatrix* WP = NM_multiply(W_ori, P_mu);
  NumericsMatrix* W = NM_multiply(P_mu, WP);
  NM_gemv(1.0, P_mu, q_ori, 0.0, q);

  size_t W_nzmax = NM_nnz(W);
  if (WP) {
    WP = NM_free(WP);
  }  // Free buffer var to avoid memory leaking

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;
  double barr_param_a = 0., e = 0.;  // for Mehrotra

  /* COMPUTATION OF A NEW STARTING POINT */
  for (size_t i = 0; i < nd; i++)
    if (i % d == 0)
      reaction[i] = 0.1;
    else
      reaction[i] = 0.01;

  for (size_t i = 0; i < nd; i++)
    if (i % d == 0)
      velocity[i] = 0.1;
    else
      velocity[i] = 0.01;

  double* s = (double*)calloc(n, sizeof(double));
  for (size_t i = 0; i < n; i++) s[i] = 0.014;

  double* d_velocity = (double*)calloc(nd, sizeof(double));
  double* d_reaction = (double*)calloc(nd, sizeof(double));
  double* d_s = (double*)calloc(n, sizeof(double));
  double* u_plus_du = data->tmp_vault_nd[no_nd++];   // for Mehrotra
  double* r_plus_dr = data->tmp_vault_nd[no_nd++];   // for Mehrotra
  double* dudr_jprod = data->tmp_vault_nd[no_nd++];  // for Mehrotra

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];
  double gmm = gmmp1 + gmmp2;

  int hasNotConverged = 1;
  int iteration = 0;
  double pinfeas = 1e300;
  // double complem = 1e300;
  // double dualgap = 1e300;
  double udotr = 1e300, uor_mu = 1e300;
  double projerr = 1e300;
  double totalresidual = 1e300, totalresidual_mu = 1e300;
  double diff_fixp = 1e300;

  double* primalConstraint = data->tmp_vault_nd[no_nd++];
  double* complemConstraint = data->tmp_vault_nd[no_nd++];
  double* fixpConstraint = (double*)calloc(n, sizeof(double));

  double kappa_mu = 1.;  // reduction of barparam
  double max_uor_2mu = 0., tmp_uor_2mu = 0.;
  double scale_sub_diff = 0.9;

  double* rhs = options->dWork;
  double* rhs_2 = (double*)calloc(size_sol, sizeof(double));
  double* sol = (double*)calloc(size_sol, sizeof(double));

  char fws = ' '; /* finish without scaling */

  /* norm of the residuals of teh second linear system */
  double LS_norm_p = 0.;  // primal feasibility
  double LS_norm_c = 0.;  // complementarity
  double LS_norm_f = 0.;  // fixed point

  long J_nzmax;

  NumericsMatrix* minus_eye_nd = NM_eye(nd);
  NM_scal(-1.0, minus_eye_nd);
  NumericsMatrix* eye_n = NM_eye(n);

  NumericsMatrix* mat_E = NM_create(NM_SPARSE, to_int(nd), to_int(n));
  int mat_E_nzmax = to_int(n);
  NM_triplet_alloc(mat_E, mat_E_nzmax);
  // NM_fill(mat_E, NM_SPARSE, nd, n, mat_E->matrix2); // Useless?
  for (size_t i = 0; i < n; ++i) {
    NM_entry(mat_E, i * d, i, 1.);
  }

  /* ---- IPM iterations ---- */
  numerics_printf_verbose(-1, "problem dimensions n, nd x m: %1i, %6i x %-6i", n, nd, m);
  switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL: {
      numerics_printf_verbose(-1, "Friction contact problem - LS solution: 3x3 no scaling\n");
      break;
    }
    default: {
      printf("ERROR\n");
    }
  }

  numerics_printf_verbose(
      -1,
      "| it  | pinfeas |  |s-ub| | |uor-mu||2max|uor-mu||   4*mu  |  u'r/n  | prj err | "
      "barpram |  alpha  |  |du|   |  |dr|   |  |ds|   | ls prim | ls comp | ls fixP |");
  numerics_printf_verbose(
      -1,
      "---------------------------------------------------------------------------------------"
      "-------------------------------------------------------------------------");

  FILE* iterates = NULL;
  //  FILE* matrixH = NULL;
  FILE* iterates_2 = NULL;
  FILE* sol_file = NULL;

  // Read problem names
  FILE* f_problem_name = fopen("problem_names.res", "r");
  char problem_name[100];
  if (!f_problem_name)
    printf("\n\nProblem names data file is not available!!! \n\n");
  else {
    fscanf(f_problem_name, "%s\n", problem_name);
  }

  char* str = (char*)malloc(200);
  strcpy(str, problem_name);
  const char* separators = "/";
  char* strToken = strtok(str, separators);
  for (int i = 0; i < 5; i++) {
    if (strToken != NULL) strToken = strtok(NULL, separators);
  }

  strToken = strtok(strToken, ".");
  // for(int i=0; i<strlen(strToken); i++)
  // {
  //   if(strToken[i] == '-') strToken[i] = '_';
  // }

  char matlab_name[100];  //, probName[100];

  // int count=0; for (int i = 13; problem_name[i] != '.'; i++) {probName[count] =
  // problem_name[i]; count++;} probName[count] = '\0'; sprintf(matlab_name, "%s.m",probName);

  // sprintf(matlab_name, "%s.m",strToken);
  sprintf(matlab_name, "iterates_local_Aqueduc.m");

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
    // iterates_2 = fopen("box466_s_u0_ub_noSubDiff.m", "w");
    iterates = fopen(matlab_name, "a+");
    // iterates = fopen(matlab_name, "w");
    fprintf(iterates, "%% data = struct;\n");
    fprintf(iterates, "data(end+1).name = \"%s\";\n", strToken);
    fprintf(iterates, "data(end).val = [\n");
    // printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates_2);
  }
  free(str);
  /* check the full criterion */
  ComputeErrorPtr computeError = &fc3d_compute_error;

  int load_starting_point = 0, save_sol_point = 0;
  if (load_starting_point) {
    sol_file = fopen("sol_data.res", "r");
    if (!sol_file)
      printf("\n\nSolution data file is not available!!! \n\n");
    else {
      // load u
      for (size_t i = 0; i < nd; i++) {
        fscanf(sol_file, "%lf ", velocity + i);
      }
      fscanf(sol_file, "\n");

      // load r
      for (size_t i = 0; i < nd; i++) {
        fscanf(sol_file, "%lf ", reaction + i);
      }
      fscanf(sol_file, "\n");

      // load s
      for (size_t i = 0; i < n; i++) {
        fscanf(sol_file, "%lf ", s + i);
      }
      fscanf(sol_file, "\n");
    }

    fclose(sol_file);

    // Sol perturbation
    for (size_t i = 0; i < nd; i++) {
      if (i % d == 0) {
        velocity[i] *= 1.2;
        reaction[i] *= 1.1;
      } else {
        velocity[i] *= 0.9;
        reaction[i] *= 0.8;
      }
    }

    for (size_t i = 0; i < n; i++) {
      s[i] *= 1.05;
    }
  }

  int jacobian_is_nan = 0;
  while (iteration < max_iter) {
    jacobian_is_nan = 0;
    switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
      /*  Build the Jacobian matrix with reducing the linear system for vector globalVelocity.
       *
       *  In the case where the NT scaling is not used, the matrix to factorize is the
       following:
       *
       *          nd           nd         n
       *      |    W           -I         E | nd
       *      |                             |
       *  J = |  Arw(u)      Arw(r)       0 | nd
       *      |                             |
       *      |    0     [0 -ub'/|ub|]    I | n
       *
       *
       * Case of a non-reduced system. Building the right-hand side related to the first system

         without NT scaling
         rhs = -
         [     Wr + q + Es - u      ]  nd        primalConstraint
         [      u o r - 2 mu e      ]  nd        complemConstraint
         [         s - |ub|         ]  n         fixpConstraint
      */
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL: {
        NumericsMatrix* Jmat = NM_create(NM_SPARSE, to_int(size_sol), to_int(size_sol));
        J_nzmax = to_csint(W_nzmax + nd + n + 2 * (d * 3 - 2) * n + 2 * n + n);
        NM_triplet_alloc(Jmat, J_nzmax);

        NumericsMatrix* arrow_r = Arrow_repr(reaction, nd, n);
        NumericsMatrix* arrow_u = Arrow_repr(velocity, nd, n);

        // Create subdiff_u
        NumericsMatrix* subdiff_u = NM_create(NM_SPARSE, to_int(n), to_int(nd));
        CS_INT subdiff_u_nzmax = to_csint(2 * n);
        NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);
        // NM_fill(subdiff_u, NM_SPARSE, n, nd, subdiff_u->matrix2); // useless?

        /* Matrix filling */
        size_t pos;
        double ub;
        scale_sub_diff = 0.95;
        for (size_t i = 0; i < n; ++i) {
          pos = i * d;
          ub = sqrt(velocity[pos + 1] * velocity[pos + 1] +
                    velocity[pos + 2] * velocity[pos + 2]);
          NM_entry(subdiff_u, i, pos + 1, -1. * scale_sub_diff * velocity[pos + 1] / ub);
          NM_entry(subdiff_u, i, pos + 2, -1. * scale_sub_diff * velocity[pos + 2] / ub);

          fixpConstraint[i] = s[i] - ub;  // fixpConstraint = s - |u_bar|
        }
        NM_insert(Jmat, W, 0, 0);
        NM_insert(Jmat, arrow_u, nd, 0);

        NM_insert(Jmat, minus_eye_nd, 0, nd);
        NM_insert(Jmat, arrow_r, nd, nd);
        NM_insert(Jmat, subdiff_u, 2 * nd, nd);

        NM_insert(Jmat, mat_E, 0, 2 * nd);
        NM_insert(Jmat, eye_n, 2 * nd, 2 * nd);

        subdiff_u = NM_free(subdiff_u);

        if (arrow_r) {
          arrow_r = NM_free(arrow_r);
        }
        if (arrow_u) {
          arrow_u = NM_free(arrow_u);
        }

        jacobian_is_nan = NM_isnan(Jmat);
        if (jacobian_is_nan) {
          numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
          Jmat = NM_free(Jmat);
          break;
        }
        // MEHROTRA scheme
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_YES) {
          primalResidual_s_for_local_friction(velocity, reaction, W, q, s, primalConstraint,
                                              &pinfeas, tol);
          JA_prod(velocity, reaction, nd, n, complemConstraint);

          cblas_dcopy(blasnd, primalConstraint, 1, rhs, 1);
          cblas_dcopy(blasnd, complemConstraint, 1, rhs + nd, 1);
          cblas_dcopy(blas_n, fixpConstraint, 1, rhs + 2 * nd, 1);

          cblas_dscal(size_blas, -1.0, rhs, 1);
          cblas_dcopy(size_blas, rhs, 1, rhs_2, 1);  // rhs_2 = old rhs

          double* rhs_tmp = (double*)calloc(size_sol, sizeof(double));
          cblas_dcopy(size_blas, rhs_2, 1, rhs_tmp, 1);
          for (size_t k = 0; k < size_sol; sol[k] = 0., k++);  // reset sol

          int max_refine = 1;
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
              SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
            max_refine = 10;

          for (int itr = 0; itr < max_refine; itr++) {
            NM_LU_solve(Jmat, rhs_tmp, 1);  // rhs_tmp = d = solution of J*d = rhs
            cblas_daxpy(size_blas, 1.0, rhs_tmp, 1, sol, 1);  // sol = x_+ = x + d
            cblas_dcopy(size_blas, rhs_2, 1, rhs_tmp, 1);     // rhs_tmp = old rhs = b
            NM_gemv(-1.0, Jmat, sol, 1.0, rhs_tmp);           // rhs_tmp = b - J*x_+
            // printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(2*nd, rhs_tmp, 1));
            if (cblas_dnrm2(size_blas, rhs_tmp, 1) <= 1e-14) {
              // printf("\nrefinement iterations = %d %8.2e\n", itr + 1, cblas_dnrm2(2 * nd +
              // n, rhs_tmp, 1));
              break;
            }
          }
          if (rhs_tmp) {
            free(rhs_tmp);
            rhs_tmp = NULL;
          }

          cblas_dcopy(blasnd, sol, 1, d_reaction, 1);
          cblas_dcopy(blasnd, sol + nd, 1, d_velocity, 1);
          cblas_dcopy(blas_n, sol + 2 * nd, 1, d_s, 1);

          alpha_primal = getStepLength(velocity, d_velocity, nd, n, 1.);
          alpha_dual = getStepLength(reaction, d_reaction, nd, n, 1.);

          if (alpha_primal < alpha_dual)
            alpha_dual = alpha_primal;
          else
            alpha_primal = alpha_dual;

          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* ----- Corrector step of Mehrotra ----- */
          cblas_dcopy(blasnd, velocity, 1, u_plus_du, 1);
          cblas_dcopy(blasnd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(blasnd, alpha_primal, d_velocity, 1, u_plus_du, 1);
          cblas_daxpy(blasnd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          barr_param = cblas_ddot(blasnd, velocity, 1, reaction, 1) / (double)n;
          barr_param_a = cblas_ddot(blasnd, u_plus_du, 1, r_plus_dr, 1) / (double)n;
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(alpha_primal, 2)) : sgmp3;
          // e = 3.;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e));
          // sigma = 0.4;

          cblas_dcopy(size_blas, rhs_2, 1, rhs, 1);  // Get back value for rhs
          JA_prod(d_velocity, d_reaction, nd, n, dudr_jprod);
          cblas_daxpy(blasnd, -1.0, dudr_jprod, 1, rhs + nd, 1);
          for (size_t k = 0; k < nd; rhs[nd + k] += 2 * sigma * barr_param, k += d);

          cblas_dcopy(size_blas, rhs, 1, rhs_2, 1);  // rhs_2 = old rhs
          // printf("\n barr_param_a = %e, e = %e, sigma = %e, gmm 1st = %e, ", barr_param_a,
          // e, sigma, gmm);
        }

        else  // Not MEHROTRA
        {
          if (iteration == 0) {
            primalResidual_s_for_local_friction(velocity, reaction, W, q, s, primalConstraint,
                                                &pinfeas, tol);
            JA_prod(velocity, reaction, nd, n, complemConstraint);
            for (size_t k = 0; k < nd; complemConstraint[k] -= 2 * barr_param, k += d);
          }

          cblas_dcopy(blasnd, primalConstraint, 1, rhs, 1);
          cblas_dcopy(blasnd, complemConstraint, 1, rhs + nd, 1);
          cblas_dcopy(blas_n, fixpConstraint, 1, rhs + 2 * nd, 1);

          cblas_dscal(size_blas, -1.0, rhs, 1);
          cblas_dcopy(size_blas, rhs, 1, rhs_2, 1);  // rhs_2 = old rhs
        }

        // SOLVE
        // NM_LU_solve(J, rhs, 1);

        double* rhs_tmp = (double*)calloc(size_sol, sizeof(double));
        cblas_dcopy(size_blas, rhs_2, 1, rhs_tmp, 1);
        for (size_t k = 0; k < size_sol; sol[k] = 0., k++);  // reset sol

        int max_refine = 1;
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
          max_refine = 10;

        for (int itr = 0; itr < max_refine; itr++) {
          NM_LU_solve(Jmat, rhs_tmp, 1);  // rhs_tmp = d = solution of J*d = rhs
          cblas_daxpy(size_blas, 1.0, rhs_tmp, 1, sol, 1);  // sol = x_+ = x + d
          cblas_dcopy(size_blas, rhs_2, 1, rhs_tmp, 1);     // rhs_tmp = old rhs = b
          NM_gemv(-1.0, Jmat, sol, 1.0, rhs_tmp);           // rhs_tmp = b - J*x_+
          // printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(2*nd, rhs_tmp, 1));
          if (cblas_dnrm2(size_blas, rhs_tmp, 1) <= 1e-14) {
            // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(2*nd+n,
            // rhs_tmp, 1));
            if (rhs_tmp) {
              free(rhs_tmp);
              rhs_tmp = NULL;
            }
            break;
          }
        }
        if (rhs_tmp) {
          free(rhs_tmp);
          rhs_tmp = NULL;
        }

        // cblas_dcopy(2*nd+n, rhs, 1, sol, 1);
        // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
        // SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
        // {
        //   double residu;
        //   NM_LU_refine(J, sol, tol/1000, 10,  &residu);
        // }
        // else
        //   NM_LU_solve(J, sol, 1);

        NM_gemv(1.0, Jmat, sol, -1.0, rhs_2);

        LS_norm_p = cblas_dnrm2(blasnd, rhs_2, 1);
        LS_norm_c = cblas_dnrm2(blasnd, rhs_2 + nd, 1);
        LS_norm_f = cblas_dnrm2(blas_n, rhs_2 + 2 * nd, 1);

        cblas_dcopy(blasnd, sol, 1, d_reaction, 1);
        cblas_dcopy(blasnd, sol + nd, 1, d_velocity, 1);
        cblas_dcopy(blas_n, sol + 2 * nd, 1, d_s, 1);

        if (Jmat) Jmat = NM_free(Jmat);

        break;
      }

      default: {
        printf("ERROR\n");
      }
    }

    if (jacobian_is_nan) {
      hasNotConverged = 2;
      break;
    }

    /* computing the affine step-length */
    // alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.99);
    // alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.99);
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * alpha_primal;
    // printf("gmm 2nd = %e\n", gmm);

    /* ----- Update variables ----- */
    if (NV_isnan(d_velocity, nd) | NV_isnan(d_reaction, nd) | NV_isnan(d_s, n)) {
      hasNotConverged = 2;
      break;
    }

    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s,
    //   d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);
    // printIteresProbMatlabFile(iteration, pinfeas, dinfeas, udotr, diff_fixp, d, n, m,
    // iterates);

    cblas_daxpy(blasnd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(blasnd, alpha_dual, d_reaction, 1, reaction, 1);
    cblas_daxpy(blas_n, alpha_primal, d_s, 1, s, 1);
    if (NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n)) {
      hasNotConverged = 2;
      break;
    }

    // if (iteration > 20)
    // {
    // printf("\n");
    // for(size_t i=0; i<nd; i+=d)
    // {
    // if (s[i/d] < 0) printf("\ns[%i] = %8.20e,\t|ub[%i]| = %8.20e\n", i/d, s[i/d], i/d,
    // cblas_dnrm2(2, velocity+i+1,1)); printf("\ns[%i] = %8.20e,\t|ub[%i]| = %8.20e\n",
    // i/d, s[i/d], i/d, cblas_dnrm2(2, velocity+i+1,1)); printf("s[%i] = %e,\t|ub[%i]| =
    // %e, \ts-ub = %e\n", i/d, s[i/d], i/d, cblas_dnrm2(2, velocity+i+1,1),
    // s[i/d]-cblas_dnrm2(2, velocity+i+1,1)); if (s[i/d] < 0) s[i/d] = cblas_dnrm2(2,
    // velocity+i+1,1); s[i/d] = cblas_dnrm2(2, velocity+i+1,1);
    // }
    // printf("\n");
    // }

    /* Computation of the values of
     - primal residual: Wr + q + Es - u
     - duality gap: u'*r / n
     - complementarity: u o r
     - projection error: r - proj(r-u)
    */
    primalResidual_s_for_local_friction(velocity, reaction, W, q, s, primalConstraint,
                                        &pinfeas, tol);
    // complem =
    complemResidualNorm(velocity, reaction, nd, n);
    udotr = cblas_ddot(blasnd, velocity, 1, reaction, 1) / (double)n;

    JA_prod(velocity, reaction, nd, n, complemConstraint);
    max_uor_2mu = 0.0;
    for (size_t k = 0; k < nd; k += d) {
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_YES) {
        // complemConstraint[k] -= 2*sigma*barr_param;
        // complemConstraint[k] += dudr_jprod[k];
        // complemConstraint[k+1] += dudr_jprod[k+1];
        // complemConstraint[k+2] += dudr_jprod[k+2];
      } else {
        complemConstraint[k] -= 2 * barr_param;
      }

      tmp_uor_2mu = cblas_dnrm2(3, complemConstraint + k, 1);

      if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
    }

    uor_mu = cblas_dnrm2(blasnd, complemConstraint, 1);

    diff_fixp = 0.;
    for (size_t i = 0; i < nd; i += d) {
      diff_fixp += (s[i / d] - cblas_dnrm2(2, velocity + i + 1, 1)) *
                   (s[i / d] - cblas_dnrm2(2, velocity + i + 1, 1));
    }
    // diff_fixp = sqrt(diff_fixp);

    // for (size_t i = 0; i<nd; i+=d)
    // {
    //   printf("s[%i] = %8.20e,\t |ub[%i]| = %8.20e\n", i/d, s[i/d], i/d, cblas_dnrm2(2,
    //   velocity+i+1, 1));
    // }
    // printf("\n");

    // projerr = projectionError(velocity, reaction, n, tol);
    NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
    NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
    double norm_q = cblas_dnrm2(blasnd, problem->q, 1);
    (*computeError)(problem, data->tmp_point->t_reaction, data->tmp_point->t_velocity, tol,
                    options, norm_q, &projerr);

    totalresidual = fmax(fmax(fmax(pinfeas, diff_fixp), 2. * max_uor_2mu), 4. * barr_param);
    totalresidual_mu = fmax(fmax(pinfeas, uor_mu), diff_fixp);

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
      fprintf(iterates,
              "%d %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e "
              "%.1e;\n",
              iteration, pinfeas, diff_fixp, uor_mu, 2. * max_uor_2mu, 4. * barr_param, udotr,
              projerr, barr_param, alpha_primal,
              fabs(d_velocity[cblas_idamax(blasnd, d_velocity, 1)]),
              fabs(d_reaction[cblas_idamax(blasnd, d_reaction, 1)]),
              fabs(d_s[cblas_idamax(blas_n, d_s, 1)]), LS_norm_p, LS_norm_c, LS_norm_f);

      if (iterates_2) {
        // printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s,
        // d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates_2);

        // store s & u0 & ub in matlab file for inspection
        fprintf(iterates_2, "s(%3i,:) = [", iteration + 1);
        for (size_t i = 0; i < n; i++) {
          fprintf(iterates_2, "%8.20e, ", s[i]);
        }
        fprintf(iterates_2, "];\n");

        fprintf(iterates_2, "u0(%3i,:) = [", iteration + 1);
        for (size_t i = 0; i < nd; i += d) {
          fprintf(iterates_2, "%8.20e, ", velocity[i]);
        }
        fprintf(iterates_2, "];\n");

        fprintf(iterates_2, "ub(%3i,:) = [", iteration + 1);
        for (size_t i = 0; i < n; i++) {
          fprintf(iterates_2, "%8.20e, ", cblas_dnrm2(2, velocity + i * d + 1, 1));
        }
        fprintf(iterates_2, "];\n");
      }
    }

    numerics_printf_verbose(-1,
                            "| %3i%c| %.1e | %.1e | %.1e |   %.1e  | %.1e | %.1e | %.1e | "
                            "%.1e | %.1e | %.1e | "
                            "%.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, pinfeas, diff_fixp, uor_mu, 2. * max_uor_2mu,
                            4. * barr_param, udotr, projerr, barr_param, alpha_primal,
                            fabs(d_velocity[cblas_idamax(blasnd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(blasnd, d_reaction, 1)]),
                            fabs(d_s[cblas_idamax(blas_n, d_s, 1)]), LS_norm_p, LS_norm_c,
                            LS_norm_f);

    if (totalresidual <= tol) {
      // numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |
      // %.1e |
      // %.1e |",
      //                         iteration, fws, pinfeas, dinfeas, diff_fixp, uor_mu,
      //                         udotr, complem, projerr, barr_param);

      // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s,
      //   d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);
      // printIteresProbMatlabFile(iteration, pinfeas, dinfeas, udotr, diff_fixp, d, n, m,
      // iterates);

      double unitur;
      for (size_t i = 0; i < n; i++) {
        unitur = cblas_ddot(3, velocity + 3 * i, 1, reaction + 3 * i, 1);
        if (unitur < 0) printf("UR NEGATIF %9.2e\n", unitur);
      }

      hasNotConverged = 0;
      break;
    }

    if (alpha_primal < 1e-8) {
      numerics_printf_verbose(-1, " fc3d_ipm_snm :: failure, step alpha is too small");
      break;
    }

    kappa_mu = 0.7;
    // if (totalresidual_mu <= 1e-8)
    if (totalresidual_mu <= 10 * barr_param &&
        options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_NO)
    // if ((barr_param > 1e-6 || totalresidual_mu < tol) && alpha_primal > 1e-1)
    // if (barr_param > 1e-11 && alpha_primal > 1e-1)
    // if (barr_param > 1e-11)
    {
      barr_param *= kappa_mu;
      // barr_param = udotr / 3;

      printf("\n");
      numerics_printf_verbose(-1,
                              "| it  | pinfeas |  |s-ub| | |uor-mu||2max|uor-mu||   4*mu  "
                              "|  u'r/n  | prj err | "
                              "barpram |  alpha  |  |du|   |  |dr|   |  |ds|   | ls prim "
                              "| ls comp | ls fixP |");
    }

    iteration++;

  }  // while loop

  free(sol);
  sol = NULL;

  /* Checking strict complementarity */
  /* For each cone i from 1 to n, one checks if u+r is in the interior of the Lorentz cone */
  /* One first computes the 3 dimensional vector somme = (u+r)/norm(u+r) */
  /* Then one checks if somme[0] > sqrt(somme[1]^2 + somme[2]^2) + ceps */

  W = NM_free(W);

  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns;
  int nsc = 0;
  int nN = 0;
  int nB = 0;
  int nR = 0;
  for (size_t i = 0; i < n; i++) {
    somme[0] = velocity[3 * i] + reaction[3 * i];
    somme[1] = velocity[3 * i + 1] + reaction[3 * i + 1];
    somme[2] = velocity[3 * i + 2] + reaction[3 * i + 2];
    dur[0] = velocity[3 * i] - reaction[3 * i];
    dur[1] = velocity[3 * i + 1] - reaction[3 * i + 1];
    dur[2] = velocity[3 * i + 2] - reaction[3 * i + 2];

    ns = somme[0] - cblas_dnrm2(2, somme + 1, 1);
    // dur = dur[0] - cblas_dnrm2(2, dur + 1, 1);
    if (ns > cesp * cblas_dnrm2(3, somme, 1)) {
      nsc += 1;
      if (dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        nN += 1;
      else if (-dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        nB += 1;
      else
        nR += 1;
    } else
      // printf("cone %i %9.2e %9.2e\n", i, somme[0], cblas_dnrm2(2,somme+1, 1));
      printf(
          "cone %04zu: (u+r)_0 = %9.2e,  |(u+r)_bar| = %9.2e, (u+r)_0 - |(u+r)_bar| = "
          "%9.40e\n",
          i, somme[0], cblas_dnrm2(2, somme + 1, 1), ns);
  }
  if (nsc < (int)n) {
    printf("Ratio of Strict complementarity solutions: %4i / %4zu = %4.2f,\n", nsc, n,
           (double)nsc / (double)n);
    printf("Provisional classification BNR: %4i %4i %4i over %4zu number of contact.\n", nB,
           nN, nR, n);
  } else
    printf("Strict complementarity satisfied: %4i / %4zu  %9.2e %9.2e  %4i %4i %4i\n", nsc, n,
           ns, ns / cblas_dnrm2(3, somme, 1), nB, nN, nR);

  // Store solution into file
  if (save_sol_point) {
    sol_file = fopen("sol_data.res", "w");

    // store u
    for (size_t i = 0; i < nd; i++) {
      fprintf(sol_file, "%8.20e ", velocity[i]);
    }
    fprintf(sol_file, "\n");

    // store r
    for (size_t i = 0; i < nd; i++) {
      fprintf(sol_file, "%8.20e ", reaction[i]);
    }
    fprintf(sol_file, "\n");

    // store s
    for (size_t i = 0; i < n; i++) {
      fprintf(sol_file, "%8.20e ", s[i]);
    }
    fprintf(sol_file, "\n");

    fclose(sol_file);
  }

  /* ----- return to original variables ------ */
  NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
  cblas_dcopy(blasnd, data->tmp_point->t_velocity, 1, velocity, 1);

  NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
  cblas_dcopy(blasnd, data->tmp_point->t_reaction, 1, reaction, 1);

  options->dparam[SICONOS_DPARAM_RESIDU] = totalresidual;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;

  if (internal_allocation) {
    fc3d_IPM_SNM_free(problem, data);
    options->solverData = NULL;
  }

  if (minus_eye_nd) minus_eye_nd = NM_free(minus_eye_nd);
  if (eye_n) eye_n = NM_free(eye_n);
  if (mat_E) {
    mat_E = NM_free(mat_E);
  }

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
    fprintf(iterates, "];\n\n");
    fclose(iterates);
    if (iterates_2) fclose(iterates_2);
  }

  if (rhs_2) free(rhs_2);
  if (sol) free(sol);

  if (d_velocity) free(d_velocity);
  if (d_reaction) free(d_reaction);
  if (s) free(s);
  if (d_s) free(d_s);
  if (fixpConstraint) free(fixpConstraint);

  *info = hasNotConverged;
}

void fc3d_ipm_snm_set_default(SolverOptions* options) {
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
      SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA] =
      SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_NO;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] =
  // SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 1;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] =
      SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO;

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 500;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 */

static int fc3d_ipm_snm_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int fc3d_ipm_snm_solve_wrap(void* problem, double* reaction, double* velocity,
                                   SolverOptions* options) {
  int info = NUMERICS_OK;
  fc3d_IPM_SNM((FrictionContactProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void fc3d_ipm_snm_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  fc3d_IPM_SNM_free((FrictionContactProblem*)problem,
                    (Gfc3d_IPM_init_data*)options->solverData);
  options->solverData = NULL;
}

REGISTER_SOLVER(FC3D_IPM_SNM, "FC3D_IPM_SNM",
                "Interior Point Method with Smoothing and Newton for 3D Friction Contact",
                fc3d_ipm_snm_init_wrap, fc3d_ipm_snm_solve_wrap, fc3d_ipm_snm_free_wrap, NULL,
                fc3d_ipm_snm_set_default, /* set_default */
                500,                      /* default_max_iter - from set_default */
                1e-10,                    /* default_tol - from set_default */
                0 /* is_local_solver */)
