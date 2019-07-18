/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "fc3d_projection.h"
//#include "gfc3d_projection.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
#include "projectionOnCone.h"
#include "SiconosLapack.h"
#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "sanitizer.h"
#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "float.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
#include "float.h"
#include "JordanAlgebra.h"
#include "cs.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"
const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_STR = "GFC3D IPM";

typedef struct {
  double * reaction_hat;
  double * reaction_k;
  double * u_hat;
  double * u_k;
  double * u;
  double * b;
}
  Gfc3d_IPM_data;


typedef struct {
    /* initial interior points */
    double * globalVelocity_k;
    double * reaction_k;
    double * velocity_k;

    /* change of variable matrix */
    NumericsMatrix* P_mu;
    NumericsMatrix* P_mu_inv;

    /* initial constant solver parameters */
    double alpha_primal0;
    double alpha_dual0;
    double sigma0;
    double mu0;
}
  Gfc3d_IPM_init_data;

/** Returns the step length for variables update in IPM [1, p. 29]
 * \param x is the initial point to update.
 * \param dx is the Newton step.
 * \param vecSize is the size of the vectors x and dx.
 * \param varsCount is the count of variables concatenated into vector x.
 * \param gamma is the safety parameter.
 * \return scalar, the step length
 *
 * \cite 1. K.C. Toh, R.H. Tutuncu, M.J. Todd,
 *          On the implementation and usage of SDPT3 - a Matlab software package
 *          for semidefinite-quadratic-linear programming, version 4.0
 *          Draft, 17 July 2006
 */
static double _getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                             const unsigned int varsCount, const double gamma);

/**
 * Returns the primal constraint vector for global fricprob ( velocity - H @ globalVelocity - w )
 * \param velocity is the vector of relative velocities.
 * \param H is the constraint matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param w is the constraint vector.
 * \param out is the result velocity - H @ globalVelocity - w vector.
 */
static void _getPrimalConstraint(const double * velocity, const NumericsMatrix * H,
                                     const double * globalVelocity, const double * w, double * out);


/**
 * Returns the dual constraint vector for global fricprob ( M @ globalVelocity + f - H @ reaction )
 * \param M is the mass matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param H is the constraint matrix.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param f is the constraint vector (vector of internal and external forces).
 * \param out os the result M @ globalVelocity + f - H @ reaction vector.
 */
static void _getDualConstraint(const NumericsMatrix * M, const double * globalVelocity,
                               const NumericsMatrix * H, const double * reaction, const double * f,
                               double * out);

/* Returns the step length for variables update in IPM */
static double _getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                             const unsigned int varsCount, const double gamma)
{
    unsigned int dimension = (int)(vecSize / varsCount);
    double * alpha_list = (double*)calloc(varsCount, sizeof(double));

    unsigned int pos;
    double ai, bi, ci, di, alpha, min_alpha;
    double *xi, *xi2, *dxi, *dxi2, *xi_dxi;

    dxi2 = (double*)calloc(dimension, sizeof(double));
    xi2 = (double*)calloc(dimension, sizeof(double));
    xi_dxi = (double*)calloc(dimension, sizeof(double));

    for (unsigned int i = 0; i < varsCount; ++i)
    {
        pos = i * dimension;
        xi = x + pos;
        dxi = dx + pos;

        NV_power2(dxi, dimension, dxi2);
        ai = dxi2[0] - NV_reduce((dxi2 + 1), dimension - 1);

        NV_prod(xi, dxi, dimension, xi_dxi);
        bi = xi_dxi[0] - NV_reduce((xi_dxi + 1), dimension - 1);

        NV_power2(xi, dimension, xi2);
        ci = xi2[0] - NV_reduce((xi2 + 1), dimension - 1);

        di = bi * bi - ai * ci;

        if (ai < 0 || (bi < 0 && ai < (bi * bi) / ci))
            alpha = ((-bi - sqrt(di)) / ai);
        else if ((fabs(ai) < 1e-12) && (bi < 0))
            alpha = (-ci / (2 * bi));
        else
            alpha = DBL_MAX;

        if (fabs(alpha) < 1e-12)
            alpha = 0.0;

        alpha_list[i] = alpha;
    }

    min_alpha = NV_min(alpha_list, varsCount);

    free(xi2);
    free(dxi2);
    free(xi_dxi);
    free(alpha_list);

    return gamma * fmin(1.0, min_alpha);
}

/* Returns the primal constraint vector for global fricprob ( velocity - H @ globalVelocity - w ) */
static void _getPrimalConstraint(const double * velocity, const NumericsMatrix * H,
                                     const double * globalVelocity, const double * w, double * out)
{
    double nd = H->size0;

    /* The memory for the result vectors should be allocated using calloc
     * since H is a sparse matrix. In other case the behaviour will be undefined.*/
    double *Hv = (double*)calloc(nd, sizeof(double));
    double *u_minus_Hv = (double*)calloc(nd, sizeof(double));

    // Hv = H @ globalVelocity
    NM_gemv(1.0, H, globalVelocity, 0.0, Hv);

    // u_minus_Hv = velocity - H @ globalVelocity
    NV_sub(velocity, Hv, nd, u_minus_Hv);

    // out = velocity - H @ globalVelocity - w
    NV_sub(u_minus_Hv, w, nd, out);

    // free allocated memory
    free(Hv);
    free(u_minus_Hv);
}

/* Returns the dual constraint vector for global fricprob ( M @ globalVelocity + f - H @ reaction ) */
static void _getDualConstraint(const NumericsMatrix * M, const double * globalVelocity,
                                   const NumericsMatrix * H, const double * reaction, const double * f,
                                   double * out)
{
    double m = H->size1;

    /* The memory for the result vectors should be allocated using calloc
     * since H is a sparse matrix. In other case the behaviour will be undefined.*/
    double *Mv = (double*)calloc(m, sizeof(double));
    double *HTr = (double*)calloc(m, sizeof(double));
    double * Mv_plus_f = (double*)calloc(m, sizeof(double));

    // Mv = M @ globalVelocity
    NM_gemv(1.0, M, globalVelocity, 0.0, Mv);

    // Mv_plus_f = M @ globalVelocity + f
    NV_add(Mv, f, m, Mv_plus_f);

    // HT = H^T
    NumericsMatrix* HT = NM_transpose(H);

    // HTr = H^T @ reaction
    NM_gemv(1.0, HT, reaction, 0.0, HTr);

    // Mv_plus_f_minus_HTr = M @ globalVelocity + f - H^T @ reaction
    NV_sub(Mv_plus_f, HTr, m, out);

    // free allocated memory
    NM_free(HT);

    free(Mv);
    free(HTr);
    free(Mv_plus_f);

    return out;
}

static double _getPrimalInfeasibility(const double * velocity, const NumericsMatrix * H, const double * globalVelocity, const double * w)
{
    double * u_Hv_w = (double*)calloc(H->size0, sizeof(double));
    _getPrimalConstraint(velocity, H, globalVelocity, w, u_Hv_w);
    double norm_inf = NV_norm_inf(u_Hv_w, H->size0);
    free(u_Hv_w);
    printf("DISPLAY w\n");
    NV_display(w, H->size0);
    double norm2 = NV_norm_2(w, H->size0);
    return norm_inf / (1 + norm2);
}

static double _getDualInfeasibility(const NumericsMatrix * M, const double * globalVelocity, const NumericsMatrix * H, const double * reaction, const double * f)
{
    double * Mv_f_HTr = (double*)calloc(H->size1, sizeof(double));
    _getDualConstraint(M, globalVelocity, H, reaction, f, Mv_f_HTr);
    double norm_inf = NV_norm_inf(Mv_f_HTr, H->size1);
    free(Mv_f_HTr);
    return norm_inf / (1 + NV_norm_2(f, H->size1));
}

static double _getComplementarInfeasibility(const double * const velocity, const double * const reaction, const unsigned int vecSize, const unsigned int varsCount)
{
    double * vr_jprod = (double*)calloc(vecSize, sizeof(double));
    JA_prod(velocity, reaction, vecSize, varsCount, vr_jprod);
    double norm2 = NV_norm_2(vr_jprod, vecSize);
    free(vr_jprod);
    return norm2 / (double)varsCount;
}

void _setError(double * error, const double pinfeas, const double dinfeas, const double complem, const double barr_param)
{
    error[0] = pinfeas;
    error[1] = dinfeas;
    error[2] = complem;
    error[3] = barr_param;
}

static int saveMatrix(NumericsMatrix* m, const char * filename)
{
    NumericsMatrix * md = NM_create(NM_DENSE, m->size0, m->size1);
    NM_to_dense(m, md);
    FILE *f;
    f = fopen(filename, "wb");
    if (!f)
        return 1;
    for (int i = 0; i < m->size0; ++i)
        for (int j = 0; j < m->size1; ++j)
            fwrite(&(md->matrix0[i+j*md->size0]), sizeof(double), 1, f);
    fclose(f);
    NM_free(md);
    return 0;
}

static int saveVector(double * vec, const unsigned int vecSize, const char * filename)
{
    FILE *f;
    f = fopen(filename, "wb");
    if (!f)
        return 1;
    for (int i = 0; i < vecSize; ++i)
        fwrite(&(vec[i]), sizeof(double), 1, f);
    fclose(f);
    return 0;
}

void gfc3d_IPM_init(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
    unsigned int m = problem->M->size0;
    unsigned int nd = problem->H->size1;
    unsigned int d = problem->dimension;
    unsigned int n = (int)(nd / d);


    if (!options->dWork || options->dWorkSize != m + nd + nd)
    {
      options->dWork = (double*)calloc(m + nd + nd, sizeof(double));
      options->dWorkSize = m + nd + nd;
    }

    /* ------------- initialize starting point ------------- */
    options->solverData=(Gfc3d_IPM_init_data *)malloc(sizeof(Gfc3d_IPM_init_data));
    Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;

    /* 1. v */
    data->globalVelocity_k = (double*)calloc(m, sizeof(double));
    for (unsigned int i = 0; i < m; ++ i)
        data->globalVelocity_k[i] = 0.01;

    /* 2. u */
    data->velocity_k = (double*)calloc(nd, sizeof(double));
    for (unsigned int i = 0; i < nd; ++ i)
    {
        data->velocity_k[i] = 0.01;
        if (i % d == 0)
            data->velocity_k[i] = 0.1;
    }

    /* 2. r */
    data->reaction_k = (double*)calloc(nd, sizeof(double));
    for (unsigned int i = 0; i < nd; ++ i)
    {
        data->reaction_k[i] = 0.0351;
        if (i % d == 0)
            data->reaction_k[i] = 0.2056;
    }

    /* ------ initialize the change of variable matrix P_mu ------- */
    data->P_mu = NM_create(NM_SPARSE, nd, nd);
    NM_triplet_alloc(data->P_mu, nd);
    data->P_mu->matrix2->origin = NSM_TRIPLET;
    for (unsigned int i = 0; i < nd; ++i)
        if (i % d == 0)
            NM_zentry(data->P_mu, i, i, 1. / problem->mu[(int)(i/d)]);
        else
            NM_zentry(data->P_mu, i, i, 1.);

    /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
    data->P_mu_inv = NM_create(NM_SPARSE, nd, nd);
    NM_triplet_alloc(data->P_mu_inv, nd);
    data->P_mu_inv->matrix2->origin = NSM_TRIPLET;
    for (unsigned int i = 0; i < nd; ++i)
        if (i % d == 0)
            NM_zentry(data->P_mu_inv, i, i, problem->mu[(int)(i/d)]);
        else
            NM_zentry(data->P_mu_inv, i, i, 1.);

    /* ------ initial parameters initialization ---------- */
    data->alpha_primal0 = 1.0;
    data->alpha_dual0 = 1.0;
    data->sigma0 = 0.0;
    data->mu0 = 1.0;
}

void gfc3d_IPM_free(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
    if (options->dWork)
    {
      free(options->dWork);
      options->dWork=NULL;
      options->dWorkSize = 0;
    }
    if (options->solverData)
    {
      Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;
      free(data->globalVelocity_k);
      data->globalVelocity_k = NULL;

      free(data->velocity_k);
      data->velocity_k = NULL;

      free(data->reaction_k);
      data->reaction_k = NULL;

      NM_free(data->P_mu);
      data->P_mu = NULL;

      NM_free(data->P_mu_inv);
      data->P_mu_inv = NULL;
    }

}

int gfc3d_IPM_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the IPM Solver\n");
  }

  solver_options_nullify(options);

  options->solverId = SICONOS_GLOBAL_FRICTION_3D_IPM;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;

  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
          SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-5;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09;

  return 0;
}

void gfc3d_IPM(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
                double* restrict velocity, double* restrict globalVelocity,
                int* restrict info, SolverOptions* restrict options)
{
    unsigned int m = problem->M->size0;
    unsigned int nd = problem->H->size1;
    unsigned int d = problem->dimension;
    unsigned int n = problem->numberOfContacts;

    // initialize solver if it is not set
    int internal_allocation=0;
    if (!options->dWork || (options->dWorkSize != m + nd + nd))
    {
      gfc3d_IPM_init(problem, options);
      internal_allocation = 1;
    }

    Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;
    NumericsMatrix *P_mu = data->P_mu;
    NumericsMatrix *P_mu_inv = data->P_mu_inv;

    NumericsMatrix *M = problem->M;
    NumericsMatrix *H_tilde = NM_transpose(problem->H);
    double *w_tilde = problem->b;
    double *w = (double*)calloc(nd, sizeof(double));
    double *f = problem->q;
    double *iden;

//    char buffer[50];
//    sprintf(buffer, "/home/maksym/Work/INRIA/M_%d_%d_%d.bin", m*n*100, M->size0, M->size0);
//    saveMatrix(M, buffer);
//    memset(buffer, 0, sizeof(buffer));

//    sprintf(buffer, "/home/maksym/Work/INRIA/H_%d_%d_%d.bin", m*n*100, H_tilde->size0, H_tilde->size1);
//    saveMatrix(H_tilde, buffer);
//    memset(buffer, 0, sizeof(buffer));

//    sprintf(buffer, "/home/maksym/Work/INRIA/w_%d_%d.bin", m*n*100, nd);
//    saveVector(w_tilde, nd, buffer);
//    memset(buffer, 0, sizeof(buffer));

//    sprintf(buffer, "/home/maksym/Work/INRIA/f_%d_%d.bin", m*n*100, m);
//    saveVector(f, m, buffer);
//    memset(buffer, 0, sizeof(buffer));

//    sprintf(buffer, "/home/maksym/Work/INRIA/mu_%d_%d.bin", m*n*100, n);
//    saveVector(problem->mu, n, buffer);
//    memset(buffer, 0, sizeof(buffer));


    NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
//    NM_display(H);
    NumericsMatrix *minus_H = NM_create(H->storageType, H->size0, H->size1);
    NM_copy(H, minus_H);
    NM_gemm(-1.0, H, NM_eye(H->size1), 0.0, minus_H);

    printf("DISPLAY w: before NM_gemv\n");
    NV_display(w, nd);
    NM_gemv(1.0, P_mu, w_tilde, 0.0, w);
    printf("DISPLAY w: agter NM_gemv\n");
    NV_display(w, nd);

    double alpha_primal = data->alpha_dual0;
    double alpha_dual = data->alpha_dual0;

    cblas_dcopy(nd, data->reaction_k, 1, reaction, 1);
    cblas_dcopy(nd, data->velocity_k, 1, velocity, 1);
    cblas_dcopy(m, data->globalVelocity_k, 1, globalVelocity, 1);

    double barr_param = data->mu0;
    double sigma = data->sigma0;

    double tol = options->dparam[SICONOS_DPARAM_TOL];
    unsigned int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

    double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
    double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
    double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
    double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
    double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

    int hasNotConverged = 1;
    size_t iteration = 0;
    double pinfeas = -1.;
    double dinfeas = -1.;
    double complem = -1.;
    double error[] = {pinfeas, dinfeas, complem, barr_param};
    double *primalConstraint, *dualConstraint, *complemConstraint;
    double *d_globalVelocity, *d_velocity, *d_reaction;
    double gmm, barr_param_a, e;

    double norm_f = cblas_dnrm2(m , f , 1);
    double norm_w = cblas_dnrm2(nd , w , 1);

    double * rhs = (double*)malloc((m + nd + nd) * sizeof(double));
    double *vr_jprod, *dvdr_jprod, *vr_prod_sub_iden, *v_plus_dv, *r_plus_dr, *gv_plus_dgv;

    vr_jprod = (double*)calloc(nd, sizeof(double));
    v_plus_dv = (double*)calloc(nd, sizeof(double));
    r_plus_dr = (double*)calloc(nd, sizeof(double));
    gv_plus_dgv = (double*)calloc(m, sizeof(double));
    vr_prod_sub_iden = (double*)calloc(nd, sizeof(double));
    dvdr_jprod = (double*)calloc(nd, sizeof(double));

    complemConstraint = (double*)calloc(nd, sizeof(double));
    primalConstraint = (double*)calloc(nd, sizeof(double));
    dualConstraint = (double*)calloc(m, sizeof(double));

    NumericsMatrix *J;
    long H_nzmax, J_nzmax;
    H_nzmax = NM_triplet(H)->nzmax;
    free(H->matrix2->triplet);
    H->matrix2->triplet = NULL;

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] ==
        SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES)
    {
      numerics_printf_verbose(1,"---- GFC3D - IPM - Problem information");
      numerics_printf_verbose(1,"---- GFC3D - IPM - 1-norm of M = %g norm of f = %g ", NM_norm_1(M), norm_f);
      numerics_printf_verbose(1,"---- GFC3D - IPM - inf-norm of M = %g ", NM_norm_inf(M));

      numerics_printf_verbose(1,"---- GFC3D - IPM - 1-norm of H = %g norm of w = %g ", NM_norm_1(problem->H), norm_w);
      numerics_printf_verbose(1,"---- GFC3D - IPM - inf-norm of H = %g ", NM_norm_inf(problem->H));
      numerics_printf_verbose(1,"---- GFC3D - IPM - M is symmetric = %i ", NM_is_symmetric(M));

      numerics_printf_verbose(1,"---- GFC3D - IPM - M size = (%i, %i) ", M->size0, M->size1);
      numerics_printf_verbose(1,"---- GFC3D - IPM - H size = (%i, %i) ", problem->H->size0, problem->H->size1);
    }

    /* ---- IPM iterations ---- */
    numerics_printf_verbose(-1, "---- GFC3D - IPM - | it  |   dgap   | pinfeas  | dinfeas  | complem  | alpha_p  | alpha_d  |  sigma   |");
    numerics_printf_verbose(-1, "---- GFC3D - IPM - ------------------------------------------------------------------------------------");
    while (hasNotConverged && (iteration < max_iter))
    {
        pinfeas = _getPrimalInfeasibility(velocity, H, globalVelocity, w);
        dinfeas = _getDualInfeasibility(M, globalVelocity, H, reaction, f);
        complem = _getComplementarInfeasibility(velocity, reaction, nd, n);

        _setError(&error, pinfeas, dinfeas, complem, barr_param);

        numerics_printf_verbose(-1, "---- GFC3D - IPM - | %3i | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
                                iteration, barr_param, pinfeas, dinfeas, complem, alpha_primal, alpha_dual, sigma);

        // check exit condition
        if (NV_max(&error, 4) < tol)
        {
            hasNotConverged = 0;
            break;
        }

        /*    1. Build the Jacobian matrix
         *
         *         m     nd       nd
         *      |  M     0      -H^T  | m
         *      |                     |
         *  J = |  0   Arw(r)  Arw(u) | nd
         *      |                     |
         *      | -H     I        0   | nd
         *
         */

        J = NM_create(NM_SPARSE, m + nd + nd, m + nd + nd);
        J_nzmax = (d * d) * (m / d) + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd;
        NM_triplet_alloc(J, J_nzmax);
        J->matrix2->origin = NSM_TRIPLET;


        NM_insert(J, M, 0, 0);
        NM_insert(J, NM_transpose(minus_H), 0, m + nd);
        NM_insert(J, Arrow_repr(reaction, nd, n), m, m);
        NM_insert(J, Arrow_repr(velocity, nd, n), m, m + nd);
        NM_insert(J, minus_H, m + nd, 0);
        NM_insert(J, NM_eye(nd), m + nd, m);

        /* 2. ---- Predictor step of Mehrotra ---- */

        /*  2.1 Build predictor right-hand side */

        _getPrimalConstraint(velocity, H, globalVelocity, w, primalConstraint);
        _getDualConstraint(M, globalVelocity, H, reaction, f, dualConstraint);
        JA_prod(velocity, reaction, nd, n, complemConstraint);

        NV_insert(rhs, m + nd + nd, dualConstraint, m, 0);
        NV_insert(rhs, m + nd + nd, complemConstraint, nd, m);
        NV_insert(rhs, m + nd + nd, primalConstraint, nd, m + nd);
        cblas_dscal(m + nd + nd, -1.0, rhs, 1);

        /* Newton system solving */
        NM_gesv_expert(J, rhs, NM_KEEP_FACTORS);

        d_globalVelocity = rhs;
        d_velocity = rhs + m;
        d_reaction = rhs + m + nd;


        cblas_dcopy(nd, d_velocity, 1, data->velocity_k, 1);
        cblas_dcopy(nd, d_reaction, 1, data->reaction_k, 1);

        alpha_primal = _getStepLength(velocity, d_velocity, nd, n, 1.);
        alpha_dual = _getStepLength(reaction, d_reaction, nd, n, 1.);
        gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

        /* ----- Corrector step of Mehrotra ----- */
        cblas_dscal(nd, alpha_primal, data->velocity_k, 1);
        cblas_dscal(nd, alpha_dual, data->reaction_k, 1);

        NV_add(velocity, data->velocity_k, nd, v_plus_dv);
        NV_add(reaction, data->reaction_k, nd, r_plus_dr);

        barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / nd;

        e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * fmin(alpha_primal, alpha_dual)) : sgmp3;
        sigma = fmin(1.0, pow(barr_param_a / barr_param, e));

        iden = JA_iden(nd, n);
        cblas_dscal(nd, 2 * barr_param * sigma, iden, 1);

        JA_prod(velocity, reaction, nd, n, vr_jprod);
        JA_prod(d_velocity, d_reaction, nd, n, dvdr_jprod);
        NV_sub(vr_jprod, iden, nd, vr_prod_sub_iden);

        free(iden);

        NV_add(vr_prod_sub_iden, dvdr_jprod, nd, complemConstraint);

        NV_insert(rhs, m + nd + nd, dualConstraint, m, 0);
        NV_insert(rhs, m + nd + nd, complemConstraint, nd, m);
        NV_insert(rhs, m + nd + nd, primalConstraint, nd, m + nd);
        cblas_dscal(m + nd + nd, -1.0, rhs, 1);

        /* Newton system solving */
        NM_gesv_expert(J, rhs, NM_KEEP_FACTORS);

        NM_free(J);

        d_globalVelocity = rhs;
        d_velocity = rhs + m;
        d_reaction = rhs + m + nd;

        alpha_primal = _getStepLength(velocity, d_velocity, nd, n, gmm);
        alpha_dual = _getStepLength(reaction, d_reaction, nd, n, gmm);


        /* ----- Update variables ----- */
        cblas_dscal(nd, alpha_primal, d_velocity, 1);
        cblas_dscal(nd, alpha_dual, d_reaction, 1);
        cblas_dscal(m, alpha_primal, d_globalVelocity, 1);

        NV_add(velocity, d_velocity, nd, v_plus_dv);
        NV_add(reaction, d_reaction, nd, r_plus_dr);
        NV_add(globalVelocity, d_globalVelocity, m, gv_plus_dgv);

        cblas_dcopy(nd, v_plus_dv, 1, velocity, 1);
        cblas_dcopy(nd, r_plus_dr, 1, reaction, 1);
        cblas_dcopy(m, gv_plus_dgv, 1, globalVelocity, 1);

        barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / nd;

        iteration++;
    }

    free(rhs);
    if (internal_allocation)
    {
      gfc3d_IPM_free(problem,options);
    }
    NM_free(H_tilde);
    NM_free(minus_H);
    NM_free(H);
    free(w);
    free(v_plus_dv);
    free(r_plus_dr);
    free(gv_plus_dgv);
    free(vr_prod_sub_iden);
    free(vr_jprod);
    free(dvdr_jprod);
    free(complemConstraint);
    free(primalConstraint);
    free(dualConstraint);

    *info = hasNotConverged;
}


void gfc3d_IPM_init_EXAMPLE(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;
  if (!options->dWork || options->dWorkSize != m+n)
  {
    options->dWork = (double*)calloc(m+n,sizeof(double));
    options->dWorkSize = m+n;
  }
  if  (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_IPM_ACCELERATION ||
       options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_IPM_ACCELERATION_AND_RESTART )
  {
    options->solverData=(Gfc3d_IPM_data *)malloc(sizeof(Gfc3d_IPM_data));
    Gfc3d_IPM_data * data = (Gfc3d_IPM_data *)options->solverData;
    data->reaction_hat =  (double*)calloc(m,sizeof(double));
    data->reaction_k =  (double*)calloc(m,sizeof(double));
    data->u_hat =  (double*)calloc(m,sizeof(double));
    data->u_k =  (double*)calloc(m,sizeof(double));
    data->u =  (double*)calloc(m,sizeof(double));
    data->b =  (double*)calloc(m,sizeof(double));
  }
}
void gfc3d_IPM_free_EXAMPLE(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  if (options->dWork)
  {
    free(options->dWork);
    options->dWork=NULL;
    options->dWorkSize = 0;
  }
  if (options->solverData)
  {
    Gfc3d_IPM_data * data = (Gfc3d_IPM_data *)options->solverData;
    free(data->reaction_hat);
    free(data->u_hat);
    free(data->reaction_k);
    free(data->u_k);
    free(data->b);
    free(data);
  }

}
static double gfc3d_ipm_select_rho(NumericsMatrix* M, NumericsMatrix* H, int * is_rho_variable, SolverOptions* restrict options)
{
  double rho=0.0;
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] ==
      SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_CONSTANT)
  {
    rho = options->dparam[SICONOS_FRICTION_3D_IPM_RHO];
  }
  else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] ==
           SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_NORM_INF)
  {
    double norm_1_M =   NM_norm_1(M);
    double norm_1_H =   NM_norm_1(H);
    if ((fabs(norm_1_H) > DBL_EPSILON) &&  (fabs(norm_1_M) > DBL_EPSILON))
      rho = norm_1_M/norm_1_H;
    else
      rho = options->dparam[SICONOS_FRICTION_3D_IPM_RHO];
  }
  else if  (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] ==
            SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_RESIDUAL_BALANCING||
            options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] ==
            SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
  {
    rho = options->dparam[SICONOS_FRICTION_3D_IPM_RHO];
    *is_rho_variable = 1 ;
  }
  return rho;
}

void gfc3d_IPM_EXAMPLE(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
                double* restrict velocity, double* restrict globalVelocity,
                int* restrict info, SolverOptions* restrict options)
{
  /* verbose=3; */
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;


  NumericsMatrix* M = NULL;
  NumericsMatrix* H = NULL;

  /* if SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE = SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE,
     we force the copy into a NM_SPARSE storageType */

  if(iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE
     && problem->M->storageType == NM_SPARSE_BLOCK)
  {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    M = NM_create(NM_SPARSE,  problem->M->size0,  problem->M->size1);
    NM_copy_to_sparse(problem->M, M);
  }
  else
  {
    M = problem->M;
  }
  if(iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE
     && problem->H->storageType == NM_SPARSE_BLOCK)
  {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    H = NM_create(NM_SPARSE,  problem->H->size0,  problem->H->size1);
    NM_copy_to_sparse(problem->H, H);
  }
  else
  {
    H = problem->H;
  }

  double* q = problem->q;
  double* mu = problem->mu;

  assert((int)H->size1 == problem->numberOfContacts * problem->dimension);
  assert((int)M->size0 == M->size1);
  assert((int)M->size0 == H->size0); /* size(velocity) ==
                                      * Htrans*globalVelocity */



  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  *info = gfc3d_checkTrivialCaseGlobal(n, q, velocity, reaction, globalVelocity, options);

  if (*info == 0)
    return;

  int contact; /* Number of the current row of blocks in M */

  double norm_q = cblas_dnrm2(n , problem->q , 1);

  double norm_b = cblas_dnrm2(m , problem->b , 1);
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] ==
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES)
  {
    numerics_printf_verbose(1,"---- GFC3D - IPM - Problem information");
    numerics_printf_verbose(1,"---- GFC3D - IPM - 1-norm of M = %g norm of q = %g ", NM_norm_1(problem->M), norm_q);
    numerics_printf_verbose(1,"---- GFC3D - IPM - inf-norm of M = %g ", NM_norm_inf(problem->M));

    numerics_printf_verbose(1,"---- GFC3D - IPM - 1-norm of H = %g norm of b = %g ", NM_norm_1(problem->H), norm_b);
    numerics_printf_verbose(1,"---- GFC3D - IPM - inf-norm of H = %g ", NM_norm_inf(problem->H));
    numerics_printf_verbose(1,"---- GFC3D - IPM -  M is symmetric = %i ", NM_is_symmetric(problem->M));
  }

  int internal_allocation=0;
  if (!options->dWork || options->dWorkSize != 2*m+n)
  {
    gfc3d_IPM_init(problem, options);
    internal_allocation = 1;
  }
  /*****  IPM iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  int is_rho_variable=0;
  double rho = gfc3d_ipm_select_rho(M, H,  &is_rho_variable, options);


  if (rho <= DBL_EPSILON)
    numerics_error("gfc3d_IPM", "dparam[SICONOS_FRICTION_3D_IPM_RHO] must be nonzero");

  /* Compute M + rho H H^T (storage in M)*/
  NumericsMatrix *Htrans =  NM_transpose(H);

  NumericsMatrix *W = NM_new();

  double eta = dparam[SICONOS_FRICTION_3D_IPM_RESTART_ETA];
  double br_tau = dparam[SICONOS_FRICTION_3D_IPM_BALANCING_RESIDUAL_TAU];
  double br_phi = dparam[SICONOS_FRICTION_3D_IPM_BALANCING_RESIDUAL_PHI];

  Gfc3d_IPM_data * data = (Gfc3d_IPM_data *)options->solverData;


  double * v = globalVelocity;

  double * u = data->u;
  double * u_k = data->u_k;
  double * u_hat =  data->u_hat;

  double * reaction_k =  data->reaction_k;
  double * reaction_hat = data->reaction_hat;

  double * b_s = data->b;

  double * tmp_m =  options->dWork;
  double * tmp_n =  &options->dWork[m];

  cblas_dscal(m, 1.0/rho, reaction, 1);

  cblas_dcopy(m , reaction , 1 , reaction_k, 1);
  cblas_dcopy(m , u , 1 , u_k, 1);

  cblas_dcopy(m , reaction , 1 , reaction_hat, 1);
  cblas_dcopy(m , u , 1 , u_hat, 1);

  double rho_k=0.0, rho_ratio=0.0;
  double e_k = INFINITY, e, alpha, r, s, residual, r_scaled, s_scaled;
  double norm_Hr=0.0, norm_HTv=0.0, norm_b_s=0.0, norm_u=0.0;
  double tau , tau_k = 1.0;
  int pos;
  double normUT;

  rho_k=rho;
  int has_rho_changed = 1;


  while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      if (has_rho_changed)
      {
        NM_copy(M, W);
        DEBUG_PRINT("copy of M: "); DEBUG_EXPR(NM_display(W));
        NM_gemm(rho, H, Htrans, 1.0, W);
        DEBUG_PRINT("M + rho H H^T: ");DEBUG_EXPR(NM_display(W));
      }

      /********************/
      /*  0 - Compute b   */
      /********************/

      cblas_dcopy(m,problem->b,1,b_s,1);
      for (contact = 0 ; contact < nc ; ++contact)
      {
        pos = contact * 3;
        normUT = sqrt(u[pos + 1] * u[pos + 1] + u[pos + 2] * u[pos + 2]);
        b_s[pos] +=  problem->mu[contact]*normUT;
      }


      /********************/
      /*  1 - Compute v */
      /********************/

      /* compute the rhs */
      /* q --> v */
      cblas_dcopy(n , q , 1 , v, 1);
      //cblas_dscal(n, -1, v,1);

      /* q +  rho H*( u -b + reaction_k) --> v */

      cblas_dcopy(m , u_hat , 1 , tmp_m, 1);
      cblas_daxpy(m, -1.0, b_s, 1, tmp_m, 1);
      cblas_daxpy(m, 1.0, reaction_hat, 1, tmp_m , 1);
      NM_gemv(rho, H, tmp_m, 1.0, v);


      DEBUG_PRINT("rhs: ");
      DEBUG_EXPR(NV_display(v,n));

      /* Linear system solver */
      /* cblas_dcopy(n , w_k , 1 , v, 1); */
      NM_gesv_expert(W,v,NM_KEEP_FACTORS);
      DEBUG_PRINT("v:");
      DEBUG_EXPR(NV_display(v,n));

      /********************/
      /*  2 - Compute u */
      /********************/

      /* H^T v_k - reaction_k + b */
      cblas_dcopy(m , b_s , 1 , u, 1);
      cblas_daxpy(m, -1.0, reaction_hat, 1, u , 1);
      NM_gemv(1.0, Htrans, v, 1.0, u);

      DEBUG_PRINT("before projection");
      DEBUG_EXPR(NV_display(u,m));

      /* Loop through the contact points */
      for (contact = 0 ; contact < nc ; ++contact)
      {
        pos = contact * 3;
        projectionOnDualCone(&u[pos], mu[contact]);
      }


      DEBUG_EXPR(NV_display(u,m));

      /**********************/
      /*  3 - Compute reaction */
      /**********************/


      /* - H^T v_k + u_k -b_s ->  reaction (residual) */
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] ==
          SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
      {
        cblas_dscal(m , 0.0, reaction, 1);
        NM_gemv(-1.0, Htrans, v, 1.0, reaction);
        norm_HTv = cblas_dnrm2(m , reaction , 1);

        cblas_daxpy(m, -1.0, b_s, 1, reaction , 1);
        cblas_daxpy(m, 1.0, u, 1, reaction , 1);
        norm_b_s =  cblas_dnrm2(m , b_s , 1);
        norm_u =  cblas_dnrm2(m , u , 1);
      }
      else
      {
        cblas_dcopy(m , u, 1 , reaction, 1);
        cblas_daxpy(m, -1.0, b_s, 1, reaction , 1);
        NM_gemv(-1.0, Htrans, v, 1.0, reaction);
      }

      r = cblas_dnrm2(m , reaction , 1);
      DEBUG_EXPR(NV_display(reaction,m));
      /* reaction_hat -  A v_k + u_k -b_s ->  xi */
      cblas_daxpy(m, 1.0, reaction_hat, 1, reaction , 1);


      /*********************************/
      /*  3 - Acceleration and restart */
      /*********************************/

      DEBUG_EXPR(NV_display(u_hat,m));
      DEBUG_EXPR(NV_display(u,m));

      cblas_dcopy(m , u_hat , 1 , tmp_m, 1);
      cblas_daxpy(m, -1.0, u, 1, tmp_m , 1);
      DEBUG_EXPR(NV_display(tmp_m,m));

      cblas_dscal(n, 0.0, tmp_n, 1);

      NM_gemv(1.0*rho, H, tmp_m, 1.0, tmp_n);
      DEBUG_EXPR(NV_display(tmp_n,n));
      s = cblas_dnrm2(n , tmp_n , 1);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] ==
          SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
      {
        cblas_dscal(n, 0.0, tmp_n, 1);
        NM_gemv(1.0*rho, H, reaction, 1.0, tmp_n);
        norm_Hr = cblas_dnrm2(n , tmp_n , 1);
      }

      e =r*r+s*s;

      DEBUG_PRINTF("residual e = %e \n", e);
      DEBUG_PRINTF("residual r = %e \n", r);
      DEBUG_PRINTF("residual s = %e \n", s);
      DEBUG_PRINTF("residual e_k = %e \n", e_k);
      DEBUG_PRINTF("eta  = %e \n", eta);
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_IPM_ACCELERATION ||
          options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_IPM_ACCELERATION_AND_RESTART)
      {
        if (e <  eta * e_k)
        {
          tau  = 0.5 *(1 +sqrt(1.0+4.0*tau_k*tau_k));
          alpha = (tau_k-1.0)/tau;

          cblas_dcopy(m , u , 1 , u_hat, 1);
          cblas_dscal(m, 1+alpha, u_hat,1);
          cblas_daxpy(m, -alpha, u_k, 1, u_hat , 1);
          DEBUG_EXPR(NV_display(u_hat,m));

          cblas_dcopy(m , reaction , 1 , reaction_hat, 1);
          cblas_dscal(m, 1+alpha, reaction_hat,1);
          cblas_daxpy(m, -alpha, reaction_k, 1, reaction_hat , 1);
          DEBUG_EXPR(NV_display(reaction_hat,m));
          DEBUG_PRINTF("tau  = %e, \t tau_k  = %e \t alpha  = %e   \n", tau, tau_k, alpha);
          numerics_printf_verbose(2, "Accelerate :tau  = %e, \t tau_k  = %e, \t alpha  = %e ", tau, tau_k, alpha);
          tau_k=tau;
          e_k=e;
        }
        else
        {
          tau_k=1.0;
          e_k = e_k /eta;
          DEBUG_PRINTF("tau_k  = %e \t alpha  = %e   \n", tau_k);
          numerics_printf_verbose(2," Restart tau_k  = %e", tau_k);
          cblas_dcopy(m , reaction_k , 1 , reaction_hat, 1);
          cblas_dcopy(m , u_k , 1 , u_hat, 1);
        }
      }
      else  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_IPM_NO_ACCELERATION)
      {
        tau_k=1.0;
        e_k = e_k /eta;
        numerics_printf_verbose(2,"Restart tau_k  = %e  \n", tau_k);
        cblas_dcopy(2*m , reaction_k , 1 , reaction_hat, 1);
        cblas_dcopy(2*m , u_k , 1 , u_hat, 1);
      }
      else
      {
        numerics_error("gfc3d_ipm", " options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ACCELERATION] value is not recognize");
      }


      rho_k = rho ;
      numerics_printf_verbose(2, "gfc3d_ipm. residuals : r  = %e, \t  s = %e", r, s);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] ==
          SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
      {
        r_scaled = r / fmax(norm_u,(fmax(norm_HTv, norm_b_s)));
        s_scaled = s / (rho*norm_Hr);
        numerics_printf_verbose(2, "gfc3d_ipm. scaling : norm_u  = %e, \t norm_HTv  = %e, \t norm_b = %e, \t", norm_u,  norm_HTv, norm_b);
        numerics_printf_verbose(2, "gfc3d_ipm. scaled residuals : r_scaled  = %e, \t  s_scaled = %e", r_scaled, s_scaled);
      }
      else
      {
        r_scaled = r;
        s_scaled = s;
      }

      if (is_rho_variable)
      {
        if (r_scaled > br_phi * s_scaled)
        {
          rho = br_tau* rho_k;
          has_rho_changed = 1;
        }
        else if (s_scaled > br_phi * r_scaled)
        {
          rho = rho_k/br_tau;
          has_rho_changed = 1;
        }
        else
        {
          /* keep the value of rho */
          has_rho_changed = 0;
        }
      }
      else
      {
        has_rho_changed = 0;
      }
      numerics_printf_verbose(2, "gfc3d_ipm. rho = %5.2e\t, rho_k = %5.2e\t ", rho, rho_k);
      rho_ratio = rho_k/rho;

      DEBUG_PRINTF("rho =%e\t,rho_k =%e \n", rho, rho_k);

      cblas_dscal(m, rho_ratio, reaction,1);
      cblas_dscal(m, rho_ratio, reaction_hat,1);

      /* Next step */
      cblas_dcopy(m , reaction , 1 , reaction_k, 1);
      cblas_dcopy(m , u , 1 , u_k, 1);

      /*********************************/
      /*  4 - Stopping criterium       */
      /*********************************/


      residual = sqrt(e);
      if (fabs(norm_q) > DBL_EPSILON)
        residual /= norm_q;

      numerics_printf_verbose(1,"---- GFC3D - IPM  - Iteration %i rho = %14.7e, residual = %14.7e, tol = %14.7e", iter, rho, residual, tolerance);

      if (residual < tolerance)
      {
        /* check the full criterion */
        cblas_dscal(m, rho, reaction, 1);
        gfc3d_compute_error(problem,  reaction, velocity, v,  tolerance, options, norm_q, &error);
        if (error < dparam[SICONOS_DPARAM_TOL])
        {
          hasNotConverged = 0;
          numerics_printf_verbose(1,"---- GFC3D - IPM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
        }
        else
        {
          numerics_printf_verbose(1,"---- GFC3D - IPM  - The tolerance on the  residual is not sufficient to reach accuracy (error =  %14.7e)", error);
          tolerance = tolerance * residual/error;
          numerics_printf_verbose(1,"---- GFC3D - IPM  - We reduce the tolerance on the residual to %14.7e", tolerance);
          cblas_dscal(m, 1.0/rho, reaction, 1);
        }
      }
      *info = hasNotConverged;
    }

  if (iter==itermax)
  {
    cblas_dscal(m, rho, reaction, 1);
    gfc3d_compute_error(problem,  reaction, velocity, v,  tolerance, options, norm_q, &error);
    numerics_printf_verbose(1,"---- GFC3D - IPM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
  }

  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;



  /***** Free memory *****/
  NM_free(W);
  NM_free(Htrans);
  if (internal_allocation)
  {
    gfc3d_IPM_free(problem,options);
  }
}



int gfc3d_IPM_setDefaultSolverOptions_EXAMPLE(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the IPM Solver\n");
  }

  options->solverId = SICONOS_GLOBAL_FRICTION_3D_IPM;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;

  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ACCELERATION] =
    SICONOS_FRICTION_3D_IPM_ACCELERATION_AND_RESTART;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] =  SICONOS_FRICTION_3D_IPM_KEEP_STORAGE;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_RHO_STRATEGY] =
    SICONOS_FRICTION_3D_IPM_RHO_STRATEGY_CONSTANT;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
    SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;


  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options->dparam[SICONOS_FRICTION_3D_IPM_RHO] = 1.0;
  options->dparam[SICONOS_FRICTION_3D_IPM_RESTART_ETA] = 0.999;
  options->dparam[SICONOS_FRICTION_3D_IPM_BALANCING_RESIDUAL_TAU]=2.0;
  options->dparam[SICONOS_FRICTION_3D_IPM_BALANCING_RESIDUAL_PHI]=10.0;

  options->internalSolvers = NULL;


  return 0;
}
