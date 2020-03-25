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


#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
#include "SiconosLapack.h"
#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "float.h"
#include "JordanAlgebra.h"
#include "CSparseMatrix.h"
#include "NumericsSparseMatrix.h"


/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"


const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_STR = "GFC3D IPM";


typedef struct
{
  double * globalVelocity;
  double * reaction;
  double * velocity;
}
IPM_starting_point;

typedef struct
{
  double * t_globalVelocity;
  double * t_reaction;
  double * t_velocity;
}
IPM_tmp_point;

typedef struct
{
  NumericsMatrix* mat;
  NumericsMatrix* inv_mat;
}
IPM_change_of_variable;

typedef struct
{
  double alpha_primal; // primal step length
  double alpha_dual;   // dual step length
  double sigma;        // centering parameter
  double barr_param;   // barrier parameter
}
IPM_internal_params;


typedef struct
{

  /* initial interior points */
  IPM_tmp_point* tmp_point;

  /* initial interior points */
  IPM_starting_point* starting_point;

  /* change of variable matrix */
  IPM_change_of_variable* P_mu;

  /* initial internal solver parameters */
  IPM_internal_params* internal_params;

  double **tmp_vault_nd;
  double **tmp_vault_m;
}
Gfc3d_IPM_init_data;


/* ------------------------- Helper functions declaration ------------------------------ */

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
static double getNewtonStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                                  const unsigned int varsCount, const double gamma);

/**
 * Returns the primal constraint vector for global fricprob ( velocity - H @ globalVelocity - w )
 * \param velocity is the vector of relative velocities.
 * \param H is the constraint matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param w is the constraint vector.
 * \param out is the result velocity - H @ globalVelocity - w vector.
 */
static void primalResidualVector(const double * velocity, NumericsMatrix * H,
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
static void dualResidualVector(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f,
                               double * out);

/**
 * Returns the Inf-norm of primal residual vecor ( velocity - H @ globalVelocity - w )
 * \param velocity is the vector of relative velocities.
 * \param H is the constraint matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param w is the constraint vector.
 */
static double primalResidualNorm(const double * velocity, NumericsMatrix * H,
                                 const double * globalVelocity, const double * w);

/**
 * Returns the Inf-norm of the dual residual vector ( M @ globalVelocity + f - H @ reaction )
 * \param M is the mass matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param H is the constraint matrix.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param f is the constraint vector (vector of internal and external forces).
 */
static double dualResidualNorm(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f);

/**
 * Returns the Inf-norm of the cemplementarity residual vector ( velocity o reaction )
 * \param velocity is the vector of relative velocities at each contact point.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param vecSize is the size of reaction and velocity vector.
 * \param varsCount is a count of variables concatenated into vectors reaction and velocity.
 */
static double complemResidualNorm(const double * const velocity, const double * const reaction,
                                  const unsigned int vecSize, const unsigned int varsCount);


/* ------------------------- Helper functions implementation ------------------------------ */


/* Returns the step length for variables update in IPM */
static double getNewtonStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                                  const unsigned int varsCount, const double gamma)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  double * alpha_list = (double*)calloc(varsCount, sizeof(double));

  unsigned int pos;
  double ai, bi, ci, di, alpha, min_alpha;
  double  *xi2, *dxi2, *xi_dxi;

  const double *dxi, *xi;

  dxi2 = (double*)calloc(dimension, sizeof(double));
  xi2 = (double*)calloc(dimension, sizeof(double));
  xi_dxi = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < varsCount; ++i)
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

    if(ai < 0 || (bi < 0 && ai < (bi * bi) / ci))
      alpha = ((-bi - sqrt(di)) / ai);
    else if((fabs(ai) < 1e-12) && (bi < 0))
      alpha = (-ci / (2 * bi));
    else
      alpha = DBL_MAX;

    if(fabs(alpha) < 1e-12)
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
static void primalResidualVector(const double * velocity, NumericsMatrix * H,
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
static void dualResidualVector(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f,
                               double * out)
{
  double m = H->size1;

  /* The memory for the result vectors should be allocated using calloc
   * since H is a sparse matrix. In other case the behaviour will be undefined.*/
  double *Mv = (double*)calloc(m, sizeof(double));
  double *HTr = (double*)calloc(m, sizeof(double));
  double * Mv_minus_f = (double*)calloc(m, sizeof(double));

  // Mv = M @ globalVelocity
  NM_gemv(1.0, M, globalVelocity, 0.0, Mv);
  //printf("Mv norm = %e\n", cblas_dnrm2(m, Mv, 1));


  // Mv_minus_f = M @ globalVelocity - f
  NV_sub(Mv, f, m, Mv_minus_f);

  // HT = H^T
  NumericsMatrix* HT = NM_transpose(H);

  // HTr = H^T @ reaction
  NM_gemv(1.0, HT, reaction, 0.0, HTr);
  //printf("HTr norm = %e\n", cblas_dnrm2(m, HTr, 1));

  // Mv_plus_f_minus_HTr = M @ globalVelocity - f - H^T @ reaction
  NV_sub(Mv_minus_f, HTr, m, out);

  // free allocated memory
  NM_clear(HT);
  free(HT);
  free(Mv);
  free(HTr);
  free(Mv_minus_f);

}

/* Returns the Inf-norm of primal residual vector ( velocity - H @ globalVelocity - w ) */
static double primalResidualNorm(const double * velocity, NumericsMatrix * H,
                                 const double * globalVelocity, const double * w)
{
  double * resid = (double*)calloc(H->size0, sizeof(double));
  primalResidualVector(velocity, H, globalVelocity, w, resid);

  double norm_inf = NV_norm_inf(resid, H->size0);
  free(resid);

  double norm2 = NV_norm_2(w, H->size0);
  return norm_inf / (1 + norm2);
}

/* Returns the Inf-norm of the dual residual vector ( M @ globalVelocity + f - H @ reaction ) */
static double dualResidualNorm(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f)
{
  double * resid = (double*)calloc(H->size1, sizeof(double));
  dualResidualVector(M, globalVelocity, H, reaction, f, resid);
  //printf("resid\n"); NV_display(resid, M->size0);
  double norm_inf = NV_norm_inf(resid, H->size1);
  free(resid);

  return norm_inf / (1 + NV_norm_2(f, H->size1));
}

/* Returns the Inf-norm of the cemplementarity residual vector ( velocity o reaction ) */
static double complemResidualNorm(const double * const velocity, const double * const reaction,
                                  const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  JA_prod(velocity, reaction, vecSize, varsCount, resid);

  double norm2 = NV_norm_2(resid, vecSize);
  free(resid);

  return norm2 / (double)varsCount;
}

static void setErrorArray(double * error, const double pinfeas, const double dinfeas,
                          const double complem, const double barr_param)
{
  error[0] = pinfeas;
  error[1] = dinfeas;
  error[2] = complem;
  error[3] = barr_param;
}

/* static int saveMatrix(NumericsMatrix* m, const char * filename) */
/* { */
/*     NumericsMatrix * md = NM_create(NM_DENSE, m->size0, m->size1); */
/*     NM_to_dense(m, md); */
/*     FILE *f; */
/*     f = fopen(filename, "wb"); */
/*     if (!f) */
/*         return 1; */
/*     for (int i = 0; i < m->size0; ++i) */
/*         for (int j = 0; j < m->size1; ++j) */
/*             fwrite(&(md->matrix0[i+j*md->size0]), sizeof(double), 1, f); */
/*     fclose(f); */
/*     NM_free(md); */
/*     return 0; */
/* } */


/* static int saveVector(double * vec, const unsigned int vecSize, const char * filename) */
/* { */
/*     FILE *f; */
/*     f = fopen(filename, "wb"); */
/*     if (!f) */
/*         return 1; */
/*     for (unsigned int i = 0; i < vecSize; ++i) */
/*         fwrite(&(vec[i]), sizeof(double), 1, f); */
/*     fclose(f); */
/*     return 0; */
/* } */



/* --------------------------- Interior-point method implementation ------------------------------ */
/*
 * Implementation contains the following functions:
 *  - gfc3d_IPM_init - initialize solver (allocate memory)
 *  - gfc3d_IPM_free - deallocate memory
 *  - gfc3d_IPM_setDefaultSolverOptions - setup default solver parameters
 *  - gfc3d_IPM - optimization method
 */
void gfc3d_IPM_init(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;

  if(!options->dWork || options->dWorkSize != (size_t)(m + nd + nd))
  {
    options->dWork = (double*)calloc(m + nd + nd, sizeof(double));
    options->dWorkSize = m + nd + nd;
  }

  /* ------------- initialize starting point ------------- */
  options->solverData=(Gfc3d_IPM_init_data *)malloc(sizeof(Gfc3d_IPM_init_data));
  Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;

  /* --------- allocate memory for tmp point ----------- */
  data->tmp_point = (IPM_tmp_point*)malloc(sizeof(IPM_tmp_point));
  data->tmp_point->t_globalVelocity = (double*)calloc(m, sizeof(double));
  data->tmp_point->t_velocity = (double*)calloc(nd, sizeof(double));
  data->tmp_point->t_reaction = (double*)calloc(nd, sizeof(double));

  /* 1. v */
  data->starting_point = (IPM_starting_point*)malloc(sizeof(IPM_starting_point));
  data->starting_point->globalVelocity = (double*)calloc(m, sizeof(double));
  for(unsigned int i = 0; i < m; ++ i)
    data->starting_point->globalVelocity[i] = 0.01;

  /* 2. u */
  data->starting_point->velocity = (double*)calloc(nd, sizeof(double));
  for(unsigned int i = 0; i < nd; ++ i)
  {
    data->starting_point->velocity[i] = 0.01;
    if(i % d == 0)
      data->starting_point->velocity[i] = 0.1;
  }

  /* 3. r */
  data->starting_point->reaction = (double*)calloc(nd, sizeof(double));
  for(unsigned int i = 0; i < nd; ++ i)
  {
    data->starting_point->reaction[i] = 0.0351;
    if(i % d == 0)
      data->starting_point->reaction[i] = 0.2056;
  }

  /* ------ initialize the change of variable matrix P_mu ------- */
  data->P_mu = (IPM_change_of_variable*)malloc(sizeof(IPM_change_of_variable));
  data->P_mu->mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->mat, nd);
  data->P_mu->mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
    if(i % d == 0)
      NM_zentry(data->P_mu->mat, i, i, 1. / problem->mu[(int)(i/d)]);
    else
      NM_zentry(data->P_mu->mat, i, i, 1.);

  /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
  data->P_mu->inv_mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->inv_mat, nd);
  data->P_mu->inv_mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
    if(i % d == 0)
      NM_zentry(data->P_mu->inv_mat, i, i, problem->mu[(int)(i/d)]);
    else
      NM_zentry(data->P_mu->inv_mat, i, i, 1.);

  /* ------ initial parameters initialization ---------- */
  data->internal_params = (IPM_internal_params*)malloc(sizeof(IPM_internal_params));
  data->internal_params->alpha_primal = 1.0;
  data->internal_params->alpha_dual = 1.0;
  data->internal_params->sigma = 0.0;
  data->internal_params->barr_param = 1.0;


  /* ----- temporary vaults initialization ------- */
  data->tmp_vault_nd = (double**)malloc(17 * sizeof(double*));
  for(unsigned int i = 0; i < 17; ++i)
    data->tmp_vault_nd[i] = (double*)calloc(nd, sizeof(double));

  data->tmp_vault_m = (double**)malloc(2 * sizeof(double*));
  for(unsigned int i = 0; i < 2; ++i)
    data->tmp_vault_m[i] = (double*)calloc(m, sizeof(double));

}

void gfc3d_IPM_free(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  if(options->dWork)
  {
    free(options->dWork);
    options->dWork=NULL;
    options->dWorkSize = 0;
  }
  if(options->solverData)
  {
    Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;

    free(data->starting_point->globalVelocity);
    data->starting_point->globalVelocity = NULL;

    free(data->starting_point->velocity);
    data->starting_point->velocity = NULL;

    free(data->starting_point->reaction);
    data->starting_point->reaction = NULL;

    free(data->starting_point);

    NM_clear(data->P_mu->mat);
    free(data->P_mu->mat);
    data->P_mu->mat = NULL;

    NM_clear(data->P_mu->inv_mat);
    free(data->P_mu->inv_mat);
    data->P_mu->inv_mat = NULL;

    free(data->P_mu);

    for(unsigned int i = 0; i < 17; ++i)
      free(data->tmp_vault_nd[i]);
    free(data->tmp_vault_nd);
    data->tmp_vault_nd = NULL;

    for(unsigned int i = 0; i < 2; ++i)
      free(data->tmp_vault_m[i]);
    free(data->tmp_vault_m);
    data->tmp_vault_m = NULL;

    free(data->tmp_point->t_globalVelocity);
    data->tmp_point->t_globalVelocity = NULL;

    free(data->tmp_point->t_velocity);
    data->tmp_point->t_velocity = NULL;

    free(data->tmp_point->t_reaction);
    data->tmp_point->t_reaction = NULL;

    free(data->tmp_point);

    free(data->internal_params);
  }

}

void gfc3d_IPM(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, double* restrict globalVelocity,
               int* restrict info, SolverOptions* restrict options)
{
  // the size of the problem detection
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;
  unsigned int n = problem->numberOfContacts;

  NumericsMatrix* M = NULL;
  NumericsMatrix* H_tilde = NULL;

  /* if SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE = SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE,
     we force the copy into a NM_SPARSE storageType */
  DEBUG_PRINTF("problem->M->storageType : %i\n",problem->M->storageType);
  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE
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
  DEBUG_PRINTF("problem->M->storageType : %i\n",problem->H->storageType);
  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE
      && problem->H->storageType == NM_SPARSE_BLOCK)
  {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    H_tilde = NM_create(NM_SPARSE,  problem->H->size1,  problem->H->size0);
    NM_copy_to_sparse(NM_transpose(problem->H), H_tilde);
  }
  else
  {
    H_tilde = NM_transpose(problem->H);
  }

  // initialize solver if it is not set
  int internal_allocation=0;
  if(!options->dWork || (options->dWorkSize != (size_t)(m + nd + nd)))
  {
    gfc3d_IPM_init(problem, options);
    internal_allocation = 1;
  }

  Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;
  NumericsMatrix *P_mu = data->P_mu->mat;
  NumericsMatrix *P_mu_inv = data->P_mu->inv_mat;

  double *w_tilde = problem->b;
  double *w = data->tmp_vault_nd[0];
  double *f = problem->q;

  double *iden;

  // change of variable
  // H_tilde --> H
  NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
  // w_tilde --> w
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  // compute -H
  NumericsMatrix *minus_H = NM_create(H->storageType, H->size0, H->size1);
  NM_copy(H, minus_H);
  NM_gemm(-1.0, H, NM_eye(H->size1), 0.0, minus_H);

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;

  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  unsigned int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  bool hasNotConverged = 1;
  unsigned int iteration = 0;
  double pinfeas = -1.;
  double dinfeas = -1.;
  double complem = -1.;
  double error[4];
  error[0] = pinfeas;
  error[1] = dinfeas;
  error[2] = complem;
  error[3] = barr_param;

  double gmm, barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  double *primalConstraint = data->tmp_vault_nd[1];
  double *dualConstraint = data->tmp_vault_m[0];
  double *complemConstraint = data->tmp_vault_nd[2];

  double *d_globalVelocity, *d_velocity, *d_reaction;

  double *rhs = options->dWork;
  double *gv_plus_dgv = data->tmp_vault_m[1];
  double *vr_jprod = data->tmp_vault_nd[3];
  double *v_plus_dv = data->tmp_vault_nd[4];
  double *r_plus_dr = data->tmp_vault_nd[5];
  double *vr_prod_sub_iden = data->tmp_vault_nd[6];
  double *dvdr_jprod = data->tmp_vault_nd[7];


  NumericsMatrix *J;
  long H_nzmax, J_nzmax;
  H_nzmax = NM_triplet(H)->nzmax;
  free(H->matrix2->triplet);
  H->matrix2->triplet = NULL;

  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] ==
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

  double * p = data->tmp_vault_nd[8];
  double * p2 = data->tmp_vault_nd[9];
  double * pinv = data->tmp_vault_nd[10];
  NumericsMatrix* Qp = NULL;
  NumericsMatrix* Qp2 = NULL;
  NumericsMatrix* Qpinv = NULL;
  double * velocity_t = data->tmp_vault_nd[11];
  double * d_velocity_t = data->tmp_vault_nd[12];
  double * d_reaction_t = data->tmp_vault_nd[13];
  double * velocity_t_inv = data->tmp_vault_nd[14];
  double * Qp_velocity_t_inv = data->tmp_vault_nd[15];
  double * tmp1 = data->tmp_vault_nd[16];

  ComputeErrorGlobalPtr computeError = NULL;

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]==
      SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
  {
    computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;
  }
  else
  {
    computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error_convex;
  }
  /* check the full criterion */
  double norm_q = cblas_dnrm2(m, problem->q, 1);
  double norm_b = cblas_dnrm2(nd, problem->b, 1);


  while(iteration < max_iter)
  {
    if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      NesterovToddVector(velocity, reaction, nd, n, p);
      JA_power2(p, nd, n, p2);
      Qp = Quad_repr(p, nd, n);
      Qp2 = Quad_repr(p2, nd, n);
      JA_inv(p, nd, n, pinv);
      Qpinv = Quad_repr(pinv, nd, n);
    }

    pinfeas = primalResidualNorm(velocity, H, globalVelocity, w);
    dinfeas = dualResidualNorm(M, globalVelocity, H, reaction, f);
    complem = complemResidualNorm(velocity, reaction, nd, n);

    setErrorArray(error, pinfeas, dinfeas, complem, barr_param);

    numerics_printf_verbose(-1, "---- GFC3D - IPM - | %3i | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
                            iteration, barr_param, pinfeas, dinfeas, complem, alpha_primal, alpha_dual, sigma);

    // check exit condition
    if(NV_max(error, 4) < tol)
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
    if(!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      NM_insert(J, Arrow_repr(reaction, nd, n), m, m);
      NM_insert(J, Arrow_repr(velocity, nd, n), m, m + nd);
    }
    else
    {
      NM_insert(J, Qp2, m, m);
      NM_insert(J, NM_eye(nd), m, m + nd);
    }
    NM_insert(J, minus_H, m + nd, 0);
    NM_insert(J, NM_eye(nd), m + nd, m);

    /* 2. ---- Predictor step of Mehrotra ---- */

    /*  2.1 Build predictor right-hand side */


    /** Correction of w to take into account the dependence
        on the tangential velocity */
    if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]==
        SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
    {
      for(unsigned int i = 0; i < nd; ++ i)
      {
        if(i % d == 0)
          w[i] = w_tilde[i]/(problem->mu[(int)(i/d)])
                 + sqrt(velocity[i+1]*velocity[i+1]+velocity[i+2]*velocity[i+2]);
      }
    }

    primalResidualVector(velocity, H, globalVelocity, w, primalConstraint);
    dualResidualVector(M, globalVelocity, H, reaction, f, dualConstraint);


    NV_insert(rhs, m + nd + nd, dualConstraint, m, 0);
    if(!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      JA_prod(velocity, reaction, nd, n, complemConstraint);
      NV_insert(rhs, m + nd + nd, complemConstraint, nd, m);
    }
    else
      NV_insert(rhs, m + nd + nd, reaction, nd, m);

    NV_insert(rhs, m + nd + nd, primalConstraint, nd, m + nd);
    cblas_dscal(m + nd + nd, -1.0, rhs, 1);

    /* Newton system solving */
    NM_gesv_expert(J, rhs, NM_KEEP_FACTORS);

    d_globalVelocity = rhs;
    d_velocity = rhs + m;
    d_reaction = rhs + m + nd;


    cblas_dcopy(nd, d_velocity, 1, data->tmp_point->t_velocity, 1);
    cblas_dcopy(nd, d_reaction, 1, data->tmp_point->t_reaction, 1);

    alpha_primal = getNewtonStepLength(velocity, d_velocity, nd, n, 1.);
    alpha_dual = getNewtonStepLength(reaction, d_reaction, nd, n, 1.);
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);
    //gmm = gmmp1 + 0.08 * fmin(alpha_primal, alpha_dual);

    /* ----- Corrector step of Mehrotra ----- */
    cblas_dscal(nd, alpha_primal, data->tmp_point->t_velocity, 1);
    cblas_dscal(nd, alpha_dual, data->tmp_point->t_reaction, 1);

    NV_add(velocity, data->tmp_point->t_velocity, nd, v_plus_dv);
    NV_add(reaction, data->tmp_point->t_reaction, nd, r_plus_dr);

    barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / nd;

    e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * fmin(alpha_primal, alpha_dual)) : sgmp3;
    sigma = fmin(1.0, pow(barr_param_a / barr_param, e));

    /* sigma = pow(barr_param_a / barr_param, 3); */
    /* sigma=0.90;  */
    if(!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      iden = JA_iden(nd, n);
      cblas_dscal(nd, 2 * barr_param * sigma, iden, 1);

      JA_prod(velocity, reaction, nd, n, vr_jprod);
      JA_prod(d_velocity, d_reaction, nd, n, dvdr_jprod);
      NV_sub(vr_jprod, iden, nd, vr_prod_sub_iden);

      free(iden);

      NV_add(vr_prod_sub_iden, dvdr_jprod, nd, complemConstraint);
    }
    else
    {
      NM_gemv(1.0, Qp, velocity, 0.0, velocity_t);
      NM_gemv(1.0, Qp, d_velocity, 0.0, d_velocity_t);
      NM_gemv(1.0, Qpinv, d_reaction, 0.0, d_reaction_t);
      JA_inv(velocity_t, nd, n, velocity_t_inv);
      NM_gemv(1.0, Qp, velocity_t_inv, 0.0, Qp_velocity_t_inv);
      cblas_dscal(nd, 2 * barr_param * sigma, Qp_velocity_t_inv, 1);
      JA_prod(d_velocity_t, d_reaction_t, nd, n, dvdr_jprod);
      NV_sub(reaction, Qp_velocity_t_inv, nd, tmp1);
      NV_add(tmp1, dvdr_jprod, nd, complemConstraint);
    }


    NV_insert(rhs, m + nd + nd, dualConstraint, m, 0);
    NV_insert(rhs, m + nd + nd, complemConstraint, nd, m);
    NV_insert(rhs, m + nd + nd, primalConstraint, nd, m + nd);
    cblas_dscal(m + nd + nd, -1.0, rhs, 1);

    /* Newton system solving */
    NM_gesv_expert(J, rhs, NM_KEEP_FACTORS);

    NM_clear(J);
    free(J);

    d_globalVelocity = rhs;
    d_velocity = rhs + m;
    d_reaction = rhs + m + nd;

    alpha_primal = getNewtonStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getNewtonStepLength(reaction, d_reaction, nd, n, gmm);


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

    /* printf("reaction:");NV_display(reaction,nd); */
    /* printf("velocity:");NV_display(velocity,nd); */

    /* for (int kk =0 ; kk < nd ; kk++) */
    /* { */
    /*   printf("r[kk]*v[kk]  = %42.32e\n", reaction[kk]*velocity[kk] ); */
    /* } */


    /* printf("barr_param  = %42.32e\n", barr_param); */
    /* printf("barr_param  = %42.32e\n", cblas_ddot(nd, reaction, 1, velocity, 1) ); */
    iteration++;
    /* ----- return to original variables ------ */
    NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
    /* cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1); */

    NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
    /* cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1); */



    double err;
    (*computeError)(problem,
                    data->tmp_point->t_reaction, data->tmp_point->t_velocity, globalVelocity,
                    tol, options,
                    norm_q, norm_b,  &err);
    numerics_printf_verbose(-1,"---- GFC3D - IPM  - Iteration %i, full error = %14.7e", iteration, err);
    // check exit condition
    if(err < tol)
    {
      hasNotConverged = 0;
      break;
    }



  }

  /* ----- return to original variables ------ */
  NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
  cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1);

  NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
  cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1);



  double err;
  (*computeError)(problem,
                  reaction, velocity, globalVelocity,
                  tol, options,
                  norm_q, norm_b,  &err);
  numerics_printf_verbose(-1,"---- GFC3D - IPM  - Iteration %i, full error = %14.7e", iteration, err);


  options->dparam[SICONOS_DPARAM_RESIDU] = NV_max(error, 4);
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;

  if(internal_allocation)
  {
    gfc3d_IPM_free(problem,options);
  }

  NM_clear(H_tilde);
  free(H_tilde);
  NM_clear(minus_H);
  free(minus_H);
  NM_clear(H);
  free(H);

  *info = hasNotConverged;
}

void gfc3d_ipm_set_default(SolverOptions* options)
{

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
    SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 1;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]=
    SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO;


  options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-05;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.8;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09;

}
