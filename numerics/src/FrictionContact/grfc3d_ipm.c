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

#include <assert.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <time.h>  // for clock()

#include "Friction_cst.h"
#include "JordanAlgebra.h"  // for JA functions
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "NumericsVector.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "cs.h"
#include "float.h"
#include "gfc3d_ipm.h"  // for primalResidual, dualResidual, ...
#include "grfc3d_Solvers.h"
#include "grfc3d_compute_error.h"  // for grfc3d_compute_error
#include "numerics_verbose.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

#ifdef WITH_MA57
#include "NM_MA57.h"
#include "lbl.h"
#endif

#define EPS 1e-40

typedef struct {
  double *globalVelocity;  // v
  double *velocity;        // u
  double *reaction;        // r
} IPM_point;

typedef struct {
  double *velocity_1;  // velocity_1 = (t, u_bar)
  double *velocity_2;  // velocity_2 = (t_prime, u_tilde)
  double *reaction_1;  // reaction_1 = (r0, r_bar)
  double *reaction_2;  // reaction_2 = (r0, r_tilde)
  double *t;
  double *t_prime;
} IPM_grfc3d_point;

typedef struct {
  /* initial interior points */
  IPM_point *starting_point;  // initial point
  IPM_point *original_point;  // original point which is not changed by the matrix P_mu
  IPM_grfc3d_point *grfc3d_point;

  /* change of variable matrix */
  IPM_change_of_variable *P_mu;

  /* initial internal solver parameters */
  IPM_internal_params *internal_params;

  double **tmp_vault_m;
  double **tmp_vault_nd;
  double **tmp_vault_n_dminus2;
  double **tmp_vault_n;
  double **tmp_vault_n_dplus1;
} Grfc3d_IPM_data;

// Note FP: these functions (getStepLength ...) are also defined in gfc3d_ipm with the same
// names. They are almost the same. These might lead to quite a mess ... To be reviewed by
// whoever is concerned with the ipm implementation.
/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double.
 */
static double grfc3d_getStepLength(const double *const x, const double *const dx,
                                   const unsigned int vecSize, const unsigned int varsCount,
                                   const double gamma) {
  int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  float_type aL, bL, cL, dL, alphaL, nxb;
  double min_alpha;

  min_alpha = 1e20;  // 1.0;

  for (unsigned int i = 0; i < varsCount; ++i) {
    pos = i * dimension;
    aL = dnrm2l(dimension - 1, dx + pos + 1);
    aL = (dx[pos] - aL) * (dx[pos] + aL);
    bL = x[pos] * dx[pos];
    for (int k = 1; k < dimension; bL -= x[pos + k] * dx[pos + k], k++);
    nxb = dnrm2l(dimension - 1, x + pos + 1);
    // cL = (x[pos] - nxb)*(x[pos] + nxb);
    cL = x[pos] - nxb;
    if (cL <= 0.)
      cL = DBL_EPSILON * (x[pos] + nxb);
    else
      cL = (x[pos] - nxb) *
           (x[pos] + nxb);  // to avoid negative number b/c of different data types
    dL = bL * bL - aL * cL;
    if (aL < 0 || (bL < 0 && dL > 0))
      if (bL > 0)
        alphaL = -(bL + sqrtl(dL)) / aL;
      else
        alphaL = cL / (-bL + sqrtl(dL));
    else if ((fabsl(aL) == 0.0) && (bL < 0))
      alphaL = -cL / bL / 2;
    else
      alphaL = DBL_MAX;
    min_alpha = ((alphaL < min_alpha) ? alphaL : min_alpha);
  }
  min_alpha = gamma * min_alpha;
  min_alpha = ((min_alpha < 1.0) ? min_alpha : 1.0);
  return min_alpha;
}

/* Rel gap = gapVal / (1 + abs(primal value) + abs(dual value)) */
static double grfc3d_relGap(NumericsMatrix *M, const double *f, const double *w,
                            const double *globalVelocity, const double *reaction,
                            const unsigned int nd, const unsigned int m, const double gapVal) {
  double *Mv = (double *)calloc(m, sizeof(double));
  double vMv, pval, dval;

  NM_gemv(0.5, M, globalVelocity, 0.0, Mv);
  vMv = cblas_ddot(m, globalVelocity, 1, Mv, 1);
  free(Mv);
  pval = vMv - cblas_ddot(m, f, 1, globalVelocity, 1);
  dval = -vMv - cblas_ddot(nd, w, 1, reaction, 1);
  return gapVal / (1 + fabs(pval) + fabs(dval));
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product
 * velocity o reaction  */
static double grfc3d_complemResidualNorm(const double *const velocity,
                                         const double *const reaction,
                                         const unsigned int vecSize,
                                         const unsigned int varsCount) {
  double *resid = (double *)calloc(vecSize, sizeof(double));
  JA_prod(velocity, reaction, vecSize, varsCount, resid);
  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  free(resid);
  return norm2;
}

/* ------------------------- Helper functions ------------------------------ */
/** Return a speacial sub-vector such that
 * the 1st element is always taken and
 * so do from i-th to j-th elements,
 * starting index is 1
 */
static void extract_vector(const double *const vec, const size_t vecSize, const int varsCount,
                           const size_t i, const size_t j, double *out) {
  assert(vec);
  assert(i >= 1);
  assert(i <= j);
  assert(out);

  size_t vec_dim = (int)(vecSize / varsCount);
  assert(j <= vec_dim);

  size_t out_dim = vec_dim - 2;
  assert(out_dim > 0);

  size_t posX = 0;
  size_t posY = 0;

  for (size_t k = 0; k < varsCount; k++) {
    out[posX++] = vec[posY];
    for (size_t l = i - 1; l < j; l++) {
      out[posX++] = vec[posY + l];
    }
    posY += vec_dim;
  }
}

/** Compute a block matrix J of form
 *      |  1     0     0                 1     0     0                    ...  |
 *      |                                                                      |
 *      |  0     1     0                 0     0     0                    ...  |
 *      |                                                                      |
 *      |  0     0     1                 0     0     0                    ...  |
 *      |                                                                      |
 *      |  0     0     0                 0     1     0                    ...  |
 *      |                                                                      |
 *      |  0     0     0                 0     0     1                    ...  |
 *  J = |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      | ...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ...  |
 */
static NumericsMatrix *compute_J_matrix(const size_t varsCount) {
  assert(varsCount > 0);

  NumericsMatrix *J = NM_create(NM_SPARSE, 5 * varsCount, 3 * varsCount * 2);
  NumericsMatrix *J_1 = NM_create(NM_SPARSE, 5, 3);
  NumericsMatrix *J_2 = NM_create(NM_SPARSE, 5, 3);

  long J_nzmax = 3 * 2 * varsCount;
  long J_1_nzmax = 3;
  long J_2_nzmax = 3;

  NM_triplet_alloc(J, J_nzmax);
  NM_triplet_alloc(J_1, J_1_nzmax);
  NM_triplet_alloc(J_2, J_1_nzmax);

  J->matrix2->origin = NSM_TRIPLET;
  J_1->matrix2->origin = NSM_TRIPLET;
  J_2->matrix2->origin = NSM_TRIPLET;

  NM_insert(J_1, NM_eye(3), 0, 0);
  NM_insert(J_2, NM_eye(1), 0, 0);
  NM_insert(J_2, NM_eye(2), 3, 1);

  for (unsigned i = 0; i < varsCount; i++) {
    NM_insert(J, J_1, i * 5, i * 3);
    NM_insert(J, J_2, i * 5, varsCount * 3 + i * 3);
  }

  return J;
}

/*
 * Return members of NT matrix F family (Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196)
 * f[in,out]: element of NT matrix matrix F
 * wf[in,out]: a coefficient related to f
 * F2[in,out] = F^2        Finv[in,out] = F^-1      Finv2[in,out] = F^-2
 * A member computed must be allocated before.
 * If an input is NULL, then this member will not be computed.
 */
static void family_of_F(const double *const x, const double *const z, const size_t vecSize,
                        const size_t varsCount, double *f, float_type *wf, NumericsMatrix *F,
                        NumericsMatrix *Finv, NumericsMatrix *F2, NumericsMatrix *Finv2) {
  size_t dimension = (size_t)(vecSize / varsCount);
  int f_NULL, wf_NULL;
  f_NULL = 0;
  wf_NULL = 0;
  NumericsMatrix *F_block, *Finv_block, *F2_block, *Finv2_block;

  if (!f)  // f is always allocated to use
  {
    f_NULL = 1;
    f = (double *)calloc(vecSize, sizeof(double));
  }
  if (!wf) {
    wf_NULL = 1;
    wf = (float_type *)calloc(varsCount, sizeof(float_type));
  }
  if (F) {
    F_block = NM_create(NM_DENSE, dimension, dimension);
  }
  if (Finv) {
    Finv_block = NM_create(NM_DENSE, dimension, dimension);
  }
  if (F2) {
    F2_block = NM_create(NM_DENSE, dimension, dimension);
  }
  if (Finv2) {
    Finv2_block = NM_create(NM_DENSE, dimension, dimension);
  }

  float_type gamx, gamz, gamf, w, w2;
  for (size_t i = 0; i < vecSize; i += dimension) {
    gamx = ld_gammal(x + i, dimension);
    gamz = ld_gammal(z + i, dimension);
    w = sqrtl(gamz / gamx);

    wf[(int)(i / dimension)] = w;
    if (F2 || Finv2) w2 = gamz / gamx;

    f[i] = z[i] / w + w * x[i];
    for (size_t j = 1; j < dimension; ++j) {
      f[i + j] = z[i + j] / w - w * x[i + j];
    }

    gamf = ld_gammal(f + i, dimension);
    for (size_t j = 0; j < dimension; ++j) {
      f[i + j] = f[i + j] / gamf;
    }

    // Construct matrix F
    if (F) {
      F_block->matrix0[0] = f[i] * w;
      for (size_t j = 1; j < dimension; ++j) {
        F_block->matrix0[j] = f[i + j] * w;
        F_block->matrix0[j * dimension] = F_block->matrix0[j];
      }
      for (size_t j = 1; j < dimension; ++j) {
        for (size_t k = 1; k < dimension; ++k)
          F_block->matrix0[j + dimension * k] =
              (j == k) ? (1 + f[i + j] * f[i + k] / (1 + f[i])) * w
                       : f[i + j] * f[i + k] / (1 + f[i]) * w;
      }
      NM_insert(F, F_block, i, i);
    }

    // Construct matrix F^-1
    if (Finv) {
      Finv_block->matrix0[0] = f[i] / w;
      for (size_t j = 1; j < dimension; ++j) {
        Finv_block->matrix0[j] = -f[i + j] / w;
        Finv_block->matrix0[j * dimension] = Finv_block->matrix0[j];
      }
      for (size_t j = 1; j < dimension; ++j) {
        for (size_t k = 1; k < dimension; ++k)
          Finv_block->matrix0[j + dimension * k] =
              (j == k) ? (1 + f[i + j] * f[i + k] / (1 + f[i])) / w
                       : f[i + j] * f[i + k] / (1 + f[i]) / w;
      }
      NM_insert(Finv, Finv_block, i, i);
    }

    // Construct matrix F^2
    if (F2) {
      F2_block->matrix0[0] = w2 * cblas_ddot(dimension, f + i, 1, f + i, 1);
      for (size_t j = 1; j < dimension; ++j) {
        F2_block->matrix0[j] = 2 * f[i] * f[i + j] * w2;
        F2_block->matrix0[j * dimension] = F2_block->matrix0[j];
      }
      for (size_t j = 1; j < dimension; ++j) {
        for (size_t k = 1; k < dimension; ++k)
          F2_block->matrix0[j + dimension * k] =
              (j == k) ? (1 + 2 * f[i + j] * f[i + k]) * w2 : 2 * f[i + j] * f[i + k] * w2;
      }
      NM_insert(F2, F2_block, i, i);
    }

    // Construct matrix F^-2
    if (Finv2) {
      Finv2_block->matrix0[0] = cblas_ddot(dimension, f + i, 1, f + i, 1) / w2;
      for (size_t j = 1; j < dimension; ++j) {
        Finv2_block->matrix0[j] = -2 * f[i] * f[i + j] / w2;
        Finv2_block->matrix0[j * dimension] = Finv2_block->matrix0[j];
      }
      for (size_t j = 1; j < dimension; ++j) {
        for (size_t k = 1; k < dimension; ++k)
          Finv2_block->matrix0[j + dimension * k] =
              (j == k) ? (1 + 2 * f[i + j] * f[i + k]) / w2 : 2 * f[i + j] * f[i + k] / w2;
      }
      NM_insert(Finv2, Finv2_block, i, i);
    }
  }

  if (f_NULL > 0) {
    free(f);
    f = NULL;
  }
  if (wf_NULL > 0) {
    free(wf);
    wf = NULL;
  }
  if (F) {
    NM_clear(F_block);
    free(F_block);
  }
  if (Finv) {
    NM_clear(Finv_block);
    free(Finv_block);
  }
  if (F2) {
    NM_clear(F2_block);
    free(F2_block);
  }
  if (Finv2) {
    NM_clear(Finv2_block);
    free(Finv2_block);
  }
}

/* [OLD VERSION] Return the matrix P^-1 where P is the matrix satisfying Jac = P*P'. Using the
 * formula F for the construction  */
static NumericsMatrix *Pinv_F(const double *const f, const double *const g,
                              const float_type *const wf, const float_type *const wg,
                              const size_t vecSize, const size_t varsCount) {
  size_t dim = (size_t)(vecSize / varsCount);  // dim must be 3
  size_t d5 = dim + 2;                         // d5 must be 5

  NumericsMatrix *out = NM_create(NM_SPARSE, 5 * varsCount, 5 * varsCount);
  NM_triplet_alloc(out, (5 + 2 * (2 * 2)) * varsCount);

  NumericsMatrix *out15 = NM_create(NM_DENSE, 1, 5);
  NumericsMatrix *out22 = NM_create(NM_DENSE, 2, 2);
  double *othor = (double *)calloc(2, sizeof(double));

  float_type coef, coef_tmp, tmp1, tmp2, tmp3, tmp4;
  coef = 1.0;
  coef_tmp = 1.0;

  for (size_t i = 0; i < varsCount; i++) {
    tmp1 = dnrm2l(dim, f + i * dim);
    tmp2 = 2 * f[i * dim] * dnrm2l(dim - 1, f + i * dim + 1) / dnrm2l(dim, f + i * dim);
    tmp3 = dnrm2l(dim, g + i * dim);
    tmp4 = 2 * g[i * dim] * dnrm2l(dim - 1, g + i * dim + 1) / dnrm2l(dim, g + i * dim);
    out15->matrix0[0] = 1. / sqrtl((tmp1 - tmp2) * (tmp1 + tmp2) / (wf[i] * wf[i]) +
                                   (tmp3 - tmp4) * (tmp3 + tmp4) / (wg[i] * wg[i]));

    coef = 2 * f[i * dim] * out15->matrix0[0] / dnrm2sqrl(dim, f + i * dim);
    for (size_t k = 1; k < dim; k++) {
      out15->matrix0[k] = coef * f[i * dim + k];
    }

    coef = 2 * g[i * dim] * out15->matrix0[0] / dnrm2sqrl(dim, g + i * dim);
    for (size_t k = 1; k < dim; k++) {
      out15->matrix0[k + 2] = coef * g[i * dim + k];
    }
    NM_insert(out, out15, i * d5, i * d5);

    // For matrix 2x2 at (1,1)
    coef = wf[i] / dnrm2sqrl(dim - 1, f + i * dim + 1);
    coef_tmp = 1. / dnrm2l(dim, f + i * dim);
    othor[0] = -f[i * dim + 2];
    othor[1] = f[i * dim + 1];
    size_t pos = 0;
    for (size_t k = 1; k < dim; k++) {
      for (size_t l = 1; l < dim; l++) {
        out22->matrix0[pos++] =
            coef * (coef_tmp * f[i * dim + k] * f[i * dim + l] + othor[k - 1] * othor[l - 1]);
      }
    }
    NM_insert(out, out22, i * d5 + 1, i * d5 + 1);

    // For matrix 2x2 at (3,3)
    coef = wg[i] / dnrm2sqrl(dim - 1, g + i * dim + 1);
    coef_tmp = 1. / dnrm2l(dim, g + i * dim);
    othor[0] = -g[i * dim + 2];
    othor[1] = g[i * dim + 1];
    pos = 0;
    for (size_t k = 1; k < dim; k++) {
      for (size_t l = 1; l < dim; l++) {
        out22->matrix0[pos++] =
            coef * (coef_tmp * g[i * dim + k] * g[i * dim + l] + othor[k - 1] * othor[l - 1]);
      }
    }
    NM_insert(out, out22, i * d5 + 3, i * d5 + 3);
  }

  NM_clear(out15);
  free(out15);
  NM_clear(out22);
  free(out22);
  free(othor);
  return out;
}

/* Return the matrix P^-1 where P is the matrix satisfying Jac = P*P'. Using the formula Qp for
 * the construction */
static NumericsMatrix *Pinv(const double *u1, const double *r1, const double *u2,
                            const double *r2, const size_t vecSize, const size_t varsCount) {
  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3
  assert(d3 == 3);
  size_t d5 = d3 + 2;  // d5 = 5

  double *x = (double *)calloc(vecSize, sizeof(double));
  double *z = (double *)calloc(vecSize, sizeof(double));
  Nesterov_Todd_vector(3, u1, r1, vecSize, varsCount, x);  // x = pinv2_bar
  Nesterov_Todd_vector(3, u2, r2, vecSize, varsCount, z);  // z = pinv2_tilde

  NumericsMatrix *out = NM_create(NM_SPARSE, 5 * varsCount, 5 * varsCount);
  NM_triplet_alloc(out, (5 + 2 * (2 * 2)) * varsCount);
  CSparseMatrix *out_triplet = out->matrix2->triplet;

  double *othor = (double *)calloc(2, sizeof(double));

  float_type p0inv = 0., data = 0., nub = 0., nrb = 0., det_u = 0., det_r = 0.;
  float_type nxb = 0., nzb = 0., nx = 0., nz = 0., nx2 = 0., nz2 = 0.;

  int id3, id5;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * d3;
    id5 = i * d5;

    nxb = dnrm2l(d3 - 1, x + id3 + 1);
    nzb = dnrm2l(d3 - 1, z + id3 + 1);
    nx = dnrm2l(d3, x + id3);
    nz = dnrm2l(d3, z + id3);
    nx2 = dnrm2sqrl(d3, x + id3);
    nz2 = dnrm2sqrl(d3, z + id3);

    p0inv = nx * nz /
            sqrtl((x[id3] - nxb) * (x[id3] - nxb) * (x[id3] + nxb) * (x[id3] + nxb) * nz2 +
                  (z[id3] - nzb) * (z[id3] - nzb) * (z[id3] + nzb) * (z[id3] + nzb) * nx2);

    // Assign data for out[0,0]
    cs_entry(out_triplet, id5, id5, p0inv);

    // Assign data for out[0,1:2]
    for (size_t k = 1; k < d3; k++) {
      data = -2. * p0inv * x[id3] * x[id3 + k] / nx2;
      cs_entry(out_triplet, id5, id5 + k, data);
    }

    // Assign data for out[0,3:4]
    for (size_t k = 1; k < d3; k++) {
      data = -2. * p0inv * z[id3] * z[id3 + k] / nz2;
      cs_entry(out_triplet, id5, id5 + 2 + k, data);
    }

    // Assign data for matrix A = out[1:2,1:2]
    othor[0] = -x[id3 + 2];
    othor[1] = x[id3 + 1];
    // Compute det(r1), det(u1)
    nrb = dnrm2l(d3 - 1, r1 + id3 + 1);
    nub = dnrm2l(d3 - 1, u1 + id3 + 1);
    det_r = (r1[id3] + nrb) * (r1[id3] - nrb);
    // det_r = r1[id3]-nrb;
    if (det_r < 0.)
      det_r = fabsl(det_r);
    else if (det_r == 0.)
      det_r = (r1[id3] + nrb) * 1e-20;

    det_u = (u1[id3] + nub) * (u1[id3] - nub);
    // det_u = u1[id3]-nub;
    if (det_u < 0.)
      det_u = fabsl(det_u);
    else if (det_u == 0.)
      det_u = (u1[id3] + nub) * 1e-20;

    for (size_t k = 1; k < d3; k++) {
      for (size_t l = 1; l < d3; l++) {
        data = (x[id3 + k] * x[id3 + l] / nx +
                powl(det_r / det_u, 0.25) * othor[k - 1] * othor[l - 1]) /
               dnrm2sqrl(d3 - 1, x + id3 + 1);
        cs_entry(out_triplet, id5 + k, id5 + l, data);
      }
    }

    // Assign data for matrix C = out[3:4,3:4]
    othor[0] = -z[id3 + 2];
    othor[1] = z[id3 + 1];
    // Compute det(r2), det(u2)
    nrb = dnrm2l(d3 - 1, r2 + id3 + 1);
    nub = dnrm2l(d3 - 1, u2 + id3 + 1);
    det_r = (r2[id3] + nrb) * (r2[id3] - nrb);
    // det_r = r2[id3]-nrb;
    if (det_r < 0.)
      det_r = fabsl(det_r);
    else if (det_r == 0.)
      det_r = (r2[id3] + nrb) * 1e-20;

    det_u = (u2[id3] + nub) * (u2[id3] - nub);
    // det_u = u2[id3]-nub;
    if (det_u < 0.)
      det_u = fabsl(det_u);
    else if (det_u == 0.)
      det_u = (u2[id3] + nub) * 1e-20;

    for (size_t k = 1; k < d3; k++) {
      for (size_t l = 1; l < d3; l++) {
        data = (z[id3 + k] * z[id3 + l] / nz +
                powl(det_r / det_u, 0.25) * othor[k - 1] * othor[l - 1]) /
               dnrm2sqrl(d3 - 1, z + id3 + 1);
        cs_entry(out_triplet, id5 + k + 2, id5 + l + 2, data);
      }
    }
  }

  free(x);
  free(z);
  free(othor);

  return out;
}

/* Return the matrix P^-1*y where P is the matrix satisfying Jac = P*P'. Using the formula Qp
 * for the construction */
static void Pinvy(const double *u1, const double *r1, const double *u2, const double *r2,
                  const size_t vecSize, const size_t varsCount, const double *y, double *out) {
  if (!out) {
    printf("\n[ERROR] Pinvy - output has not been yet allocated.\n");
    return;
  }

  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3
  assert(d3 == 3);
  size_t d5 = d3 + 2;  // d5 = 5

  double *x = (double *)calloc(vecSize, sizeof(double));
  double *z = (double *)calloc(vecSize, sizeof(double));
  Nesterov_Todd_vector(3, u1, r1, vecSize, varsCount, x);  // x = pinv2_bar
  Nesterov_Todd_vector(3, u2, r2, vecSize, varsCount, z);  // z = pinv2_tilde

  double *othor = (double *)calloc(2, sizeof(double));

  float_type p0inv = 0., data = 0., nub = 0., nrb = 0., det_u = 0., det_r = 0.;
  float_type nxb = 0., nzb = 0., nx = 0., nz = 0., nx2 = 0., nz2 = 0.;

  int id3, id5;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * d3;
    id5 = i * d5;

    nxb = dnrm2l(d3 - 1, x + id3 + 1);
    nzb = dnrm2l(d3 - 1, z + id3 + 1);
    nx = dnrm2l(d3, x + id3);
    nz = dnrm2l(d3, z + id3);
    nx2 = dnrm2sqrl(d3, x + id3);
    nz2 = dnrm2sqrl(d3, z + id3);

    p0inv = nx * nz /
            sqrtl((x[id3] - nxb) * (x[id3] - nxb) * (x[id3] + nxb) * (x[id3] + nxb) * nz2 +
                  (z[id3] - nzb) * (z[id3] - nzb) * (z[id3] + nzb) * (z[id3] + nzb) * nx2);

    // out[0,0]
    out[id5] =
        p0inv *
        (y[id5] - 2. * x[id3] * cblas_ddot(d3 - 1, x + id3 + 1, 1, y + id5 + 1, 1) / nx2 -
         2. * z[id3] * cblas_ddot(d3 - 1, z + id3 + 1, 1, y + id5 + 3, 1) / nz2);

    // out[1:2]
    othor[0] = -x[id3 + 2];
    othor[1] = x[id3 + 1];
    // Compute det(r1), det(u1)
    nrb = dnrm2l(d3 - 1, r1 + id3 + 1);
    nub = dnrm2l(d3 - 1, u1 + id3 + 1);
    det_r = (r1[id3] + nrb) * (r1[id3] - nrb);
    // det_r = r1[id3]-nrb;
    if (det_r < 0.)
      det_r = fabsl(det_r);
    else if (det_r == 0.)
      det_r = (r1[id3] + nrb) * 1e-20;

    det_u = (u1[id3] + nub) * (u1[id3] - nub);
    // det_u = u1[id3]-nub;
    if (det_u < 0.)
      det_u = fabsl(det_u);
    else if (det_u == 0.)
      det_u = (u1[id3] + nub) * 1e-20;

    for (size_t k = 1; k < d3; k++) {
      out[id5 + k] = (cblas_ddot(d3 - 1, x + id3 + 1, 1, y + id5 + 1, 1) * x[id3 + k] / nx +
                      powl(det_r / det_u, 0.25) *
                          cblas_ddot(d3 - 1, othor, 1, y + id5 + 1, 1) * othor[k - 1]) /
                     dnrm2sqrl(d3 - 1, x + id3 + 1);
    }

    // out[3:4]
    othor[0] = -z[id3 + 2];
    othor[1] = z[id3 + 1];
    // Compute det(r2), det(u2)
    nrb = dnrm2l(d3 - 1, r2 + id3 + 1);
    nub = dnrm2l(d3 - 1, u2 + id3 + 1);
    det_r = (r2[id3] + nrb) * (r2[id3] - nrb);
    // det_r = r2[id3]-nrb;
    if (det_r < 0.)
      det_r = fabsl(det_r);
    else if (det_r == 0.)
      det_r = (r2[id3] + nrb) * 1e-20;

    det_u = (u2[id3] + nub) * (u2[id3] - nub);
    // det_u = u2[id3]-nub;
    if (det_u < 0.)
      det_u = fabsl(det_u);
    else if (det_u == 0.)
      det_u = (u2[id3] + nub) * 1e-20;

    for (size_t k = 1; k < d3; k++) {
      out[id5 + k + 2] =
          (cblas_ddot(d3 - 1, z + id3 + 1, 1, y + id5 + 3, 1) * z[id3 + k] / nz +
           powl(det_r / det_u, 0.25) * cblas_ddot(d3 - 1, othor, 1, y + id5 + 3, 1) *
               othor[k - 1]) /
          dnrm2sqrl(d3 - 1, z + id3 + 1);
    }
  }

  free(x);
  free(z);
  free(othor);
}

/* Return the matrix (P^-1)'*y where P is the matrix satisfying Jac = P*P'. Using the formula
 * Qp for the construction */
static void PinvTy(const double *u1, const double *r1, const double *u2, const double *r2,
                   const size_t vecSize, const size_t varsCount, const double *y,
                   double *out) {
  if (!out) {
    printf("\n[ERROR] PinvTy - output has not been yet allocated.\n");
    return;
  }

  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3
  assert(d3 == 3);
  size_t d5 = d3 + 2;  // d5 = 5

  double *x = (double *)calloc(vecSize, sizeof(double));
  double *z = (double *)calloc(vecSize, sizeof(double));
  Nesterov_Todd_vector(3, u1, r1, vecSize, varsCount, x);  // x = pinv2_bar
  Nesterov_Todd_vector(3, u2, r2, vecSize, varsCount, z);  // z = pinv2_tilde

  double *othor = (double *)calloc(2, sizeof(double));

  float_type p0inv = 0., data = 0., nub = 0., nrb = 0., det_u = 0., det_r = 0.;
  float_type nxb = 0., nzb = 0., nx = 0., nz = 0., nx2 = 0., nz2 = 0.;

  int id3, id5;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * d3;
    id5 = i * d5;

    nxb = dnrm2l(d3 - 1, x + id3 + 1);
    nzb = dnrm2l(d3 - 1, z + id3 + 1);
    nx = dnrm2l(d3, x + id3);
    nz = dnrm2l(d3, z + id3);
    nx2 = dnrm2sqrl(d3, x + id3);
    nz2 = dnrm2sqrl(d3, z + id3);

    p0inv = nx * nz /
            sqrtl((x[id3] - nxb) * (x[id3] - nxb) * (x[id3] + nxb) * (x[id3] + nxb) * nz2 +
                  (z[id3] - nzb) * (z[id3] - nzb) * (z[id3] + nzb) * (z[id3] + nzb) * nx2);

    // out[0,0]
    out[id5] = p0inv * y[id5];

    // out[1:2]
    othor[0] = -x[id3 + 2];
    othor[1] = x[id3 + 1];
    // Compute det(r1), det(u1)
    nrb = dnrm2l(d3 - 1, r1 + id3 + 1);
    nub = dnrm2l(d3 - 1, u1 + id3 + 1);
    det_r = (r1[id3] + nrb) * (r1[id3] - nrb);
    // det_r = r1[id3]-nrb;
    if (det_r < 0.)
      det_r = fabsl(det_r);
    else if (det_r == 0.)
      det_r = (r1[id3] + nrb) * 1e-20;

    det_u = (u1[id3] + nub) * (u1[id3] - nub);
    // det_u = u1[id3]-nub;
    if (det_u < 0.)
      det_u = fabsl(det_u);
    else if (det_u == 0.)
      det_u = (u1[id3] + nub) * 1e-20;

    for (size_t k = 1; k < d3; k++) {
      out[id5 + k] = -2. * x[id3] * y[id5] * p0inv * x[id3 + k] / nx2 +
                     (cblas_ddot(d3 - 1, x + id3 + 1, 1, y + id5 + 1, 1) * x[id3 + k] / nx +
                      powl(det_r / det_u, 0.25) *
                          cblas_ddot(d3 - 1, othor, 1, y + id5 + 1, 1) * othor[k - 1]) /
                         dnrm2sqrl(d3 - 1, x + id3 + 1);
    }

    // out[3:4]
    othor[0] = -z[id3 + 2];
    othor[1] = z[id3 + 1];
    // Compute det(r2), det(u2)
    nrb = dnrm2l(d3 - 1, r2 + id3 + 1);
    nub = dnrm2l(d3 - 1, u2 + id3 + 1);
    det_r = (r2[id3] + nrb) * (r2[id3] - nrb);
    // det_r = r2[id3]-nrb;
    if (det_r < 0.)
      det_r = fabsl(det_r);
    else if (det_r == 0.)
      det_r = (r2[id3] + nrb) * 1e-20;

    det_u = (u2[id3] + nub) * (u2[id3] - nub);
    // det_u = u2[id3]-nub;
    if (det_u < 0.)
      det_u = fabsl(det_u);
    else if (det_u == 0.)
      det_u = (u2[id3] + nub) * 1e-20;

    for (size_t k = 1; k < d3; k++) {
      out[id5 + k + 2] =
          -2. * z[id3] * y[id5] * p0inv * z[id3 + k] / nz2 +
          (cblas_ddot(d3 - 1, z + id3 + 1, 1, y + id5 + 3, 1) * z[id3 + k] / nz +
           powl(det_r / det_u, 0.25) * cblas_ddot(d3 - 1, othor, 1, y + id5 + 3, 1) *
               othor[k - 1]) /
              dnrm2sqrl(d3 - 1, z + id3 + 1);
    }
  }

  free(x);
  free(z);
  free(othor);
}

/* Return the matrix P^-1*H where P is the matrix satisfying Jac = P*P'. Using the formula Qp.
 */
static NumericsMatrix *multiply_PinvH(const double *u1, const double *r1, const double *u2,
                                      const double *r2, const size_t vecSize,
                                      const size_t varsCount, NumericsMatrix *H) {
  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3
  assert(d3 == 3);
  size_t d5 = d3 + 2;  // d5 = 5
  size_t id3 = 0, id5 = 0;

  double *x = (double *)calloc(vecSize, sizeof(double));
  double *z = (double *)calloc(vecSize, sizeof(double));
  Nesterov_Todd_vector(3, u1, r1, vecSize, varsCount, x);  // x = pinv2_bar
  Nesterov_Todd_vector(3, u2, r2, vecSize, varsCount, z);  // z = pinv2_tilde

  double *othor1 = (double *)calloc(2, sizeof(double));
  double *othor2 = (double *)calloc(2, sizeof(double));

  float_type p0inv = 0., data = 0.;
  float_type nub = 0., nrb = 0., det_u1 = 0., det_r1 = 0., det_u2 = 0., det_r2 = 0.;
  float_type nxb = 0., nzb = 0.;

  NumericsMatrix *out = NM_new();

  if (H->storageType != NM_SPARSE) {
    fprintf(stderr,
            "Numerics, GRFC3D IPM, multiply_PinvH failed, only accept for NM_SPARSE of H.\n");
    exit(EXIT_FAILURE);
  }

  CSparseMatrix *H_csc = NM_csc(H);
  CSparseMatrix *out_csc =
      cs_spalloc(H->size0, H->size1, H_csc->nzmax, 1, 0); /* allocate result */

  CS_INT *Hp, *Hi;
  Hp = H_csc->p;
  Hi = H_csc->i;

  CS_INT *outp, *outi;
  outp = out_csc->p;

  CS_ENTRY *Hx, *outx;
  Hx = H_csc->x;

  CS_INT nz = 0;

  bool multiplied = 0;
  for (CS_INT k = 0; k < H->size1; k++)  // traverse all cols of H
  {
    outp[k] = nz; /* column k of out starts here */

    /* reallocate if needed */
    if (nz + H->size1 > out_csc->nzmax &&
        !cs_sprealloc(out_csc, 2 * (out_csc->nzmax) + H->size1)) {
      return NULL; /* out of memory */
    }
    outi = out_csc->i;
    outx = out_csc->x; /* out->i and out->x may be reallocated */

    for (size_t i = 0; i < varsCount; i++)  // traverse all blocks[5x5] of Pinv
    {
      id3 = i * d3;
      id5 = i * d5;

      // p0inv = 1./sqrtl( dnrm2sqrl(d3,x+id3) + dnrm2sqrl(d3,z+id3)
      // - 4.*x[id3]*x[id3]*dnrm2sqrl(d3-1,x+id3+1)/dnrm2sqrl(d3,x+id3)
      // - 4.*z[id3]*z[id3]*dnrm2sqrl(d3-1,z+id3+1)/dnrm2sqrl(d3,z+id3) ); p0inv =
      // dnrm2l(d3,x+id3)*dnrm2l(d3,z+id3)/sqrtl((x[id3]-dnrm2l(d3-1,x+id3+1))*(x[id3]-dnrm2l(d3-1,x+id3+1))*(x[id3]+dnrm2l(d3-1,x+id3+1))*(x[id3]+dnrm2l(d3-1,x+id3+1))*dnrm2sqrl(d3,z+id3)
      // +
      // (z[id3]-dnrm2l(d3-1,z+id3+1))*(z[id3]-dnrm2l(d3-1,z+id3+1))*(z[id3]+dnrm2l(d3-1,z+id3+1))*(z[id3]+dnrm2l(d3-1,z+id3+1))*dnrm2sqrl(d3,x+id3));
      nxb = dnrm2l(d3 - 1, x + id3 + 1);
      nzb = dnrm2l(d3 - 1, z + id3 + 1);
      p0inv = dnrm2l(d3, x + id3) * dnrm2l(d3, z + id3) /
              sqrtl((x[id3] - nxb) * (x[id3] - nxb) * (x[id3] + nxb) * (x[id3] + nxb) *
                        dnrm2sqrl(d3, z + id3) +
                    (z[id3] - nzb) * (z[id3] - nzb) * (z[id3] + nzb) * (z[id3] + nzb) *
                        dnrm2sqrl(d3, x + id3));

      othor1[0] = -x[id3 + 2];
      othor1[1] = x[id3 + 1];
      // Compute det(r1), det(u1)
      nrb = dnrm2l(d3 - 1, r1 + id3 + 1);
      nub = dnrm2l(d3 - 1, u1 + id3 + 1);
      // det_r1 = (r1[id3]+nrb)*(r1[id3]-nrb);
      det_r1 = r1[id3] - nrb;
      if (det_r1 < 0.)
        det_r1 = (r1[id3] + nrb) * DBL_EPSILON;
      else
        det_r1 = (r1[id3] + nrb) * (r1[id3] - nrb);

      // det_u1 = (u1[id3]+nub)*(u1[id3]-nub);
      det_u1 = u1[id3] - nub;
      if (det_u1 <= 0.)
        det_u1 = (u1[id3] + nub) * DBL_EPSILON;
      else
        det_u1 = (u1[id3] + nub) * (u1[id3] - nub);

      othor2[0] = -z[id3 + 2];
      othor2[1] = z[id3 + 1];
      // Compute det(r2), det(u2)
      nrb = dnrm2l(d3 - 1, r2 + id3 + 1);
      nub = dnrm2l(d3 - 1, u2 + id3 + 1);
      // det_r2 = (r2[id3]+nrb)*(r2[id3]-nrb);
      det_r2 = r2[id3] - nrb;
      if (det_r2 < 0.)
        det_r2 = (r2[id3] + nrb) * DBL_EPSILON;
      else
        det_r2 = (r2[id3] + nrb) * (r2[id3] - nrb);

      // det_u2 = (u2[id3]+nub)*(u2[id3]-nub);
      det_u2 = u2[id3] - nub;
      if (det_u2 <= 0.)
        det_u2 = (u2[id3] + nub) * DBL_EPSILON;
      else
        det_u2 = (u2[id3] + nub) * (u2[id3] - nub);

      for (size_t j = 0; j < d5; j++)  // traverse all rows-block of Pinv
      {
        /* multiplication and storage */
        if (j == 0)  // 1st row of P^-1
        {
          outx[nz] = 0.;
          for (CS_INT p = Hp[k]; p < Hp[k + 1]; p++)  // traverse all existential rows of H
          {
            // rows of H such that they belongs to each block of P^-1
            if (Hi[p] == id5) {
              outx[nz] += Hx[p];
              multiplied = 1;
            }

            if (Hi[p] == id5 + 1) {
              outx[nz] += -2. * x[id3] * x[id3 + 1] * Hx[p] / dnrm2sqrl(d3, x + id3);
              multiplied = 1;
            }
            if (Hi[p] == id5 + 2) {
              outx[nz] += -2. * x[id3] * x[id3 + 2] * Hx[p] / dnrm2sqrl(d3, x + id3);
              multiplied = 1;
            }

            if (Hi[p] == id5 + 3) {
              outx[nz] += -2. * z[id3] * z[id3 + 1] * Hx[p] / dnrm2sqrl(d3, z + id3);
              multiplied = 1;
            }
            if (Hi[p] == id5 + 4) {
              outx[nz] += -2. * z[id3] * z[id3 + 2] * Hx[p] / dnrm2sqrl(d3, z + id3);
              multiplied = 1;
            }

          }  // end rows of H
          outx[nz] *= p0inv;
          if (multiplied) {
            outi[nz++] = j + id5;
            multiplied = 0;
          }
        }  // end 1st row of P^-1

        else if (j == 1 || j == 2)  // A^-1/2 of P^-1
        {
          outx[nz] = 0.;
          for (CS_INT p = Hp[k]; p < Hp[k + 1]; p++)  // traverse all existential rows of H
          {
            // get the rows of H that belongs to A^-1/2 for each block
            if (Hi[p] == id5 + 1) {
              multiplied = 1;
              if (j == 1)
                outx[nz] += (x[id3 + 1] * x[id3 + 1] / dnrm2l(d3, x + id3) +
                             powl(det_r1 / det_u1, 0.25) * othor1[0] * othor1[0]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, x + id3 + 1);  // 1st row of A^-1/2
              if (j == 2)
                outx[nz] += (x[id3 + 2] * x[id3 + 1] / dnrm2l(d3, x + id3) +
                             powl(det_r1 / det_u1, 0.25) * othor1[1] * othor1[0]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, x + id3 + 1);  // 2nd row
            }

            if (Hi[p] == id5 + 2) {
              multiplied = 1;
              if (j == 1)
                outx[nz] += (x[id3 + 1] * x[id3 + 2] / dnrm2l(d3, x + id3) +
                             powl(det_r1 / det_u1, 0.25) * othor1[0] * othor1[1]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, x + id3 + 1);  // 1st row of A^-1/2
              if (j == 2)
                outx[nz] += (x[id3 + 2] * x[id3 + 2] / dnrm2l(d3, x + id3) +
                             powl(det_r1 / det_u1, 0.25) * othor1[1] * othor1[1]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, x + id3 + 1);  // 2nd row
            }
          }  // end rows of H
          if (multiplied) {
            outi[nz++] = j + id5;
            multiplied = 0;
          }
        }  // end A^-1/2

        else if (j == 3 || j == 4)  // C^-1/2 of P^-1
        {
          outx[nz] = 0.;
          for (CS_INT p = Hp[k]; p < Hp[k + 1]; p++)  // traverse all existential rows of H
          {
            // get the rows of H that belongs to C^-1/2 for each block
            if (Hi[p] == id5 + 3) {
              multiplied = 1;
              if (j == 3)
                outx[nz] += (z[id3 + 1] * z[id3 + 1] / dnrm2l(d3, z + id3) +
                             powl(det_r2 / det_u2, 0.25) * othor2[0] * othor2[0]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, z + id3 + 1);  // 1st row of C^-1/2
              if (j == 4)
                outx[nz] += (z[id3 + 2] * z[id3 + 1] / dnrm2l(d3, z + id3) +
                             powl(det_r2 / det_u2, 0.25) * othor2[1] * othor2[0]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, z + id3 + 1);  // 2nd row
            }

            if (Hi[p] == id5 + 4) {
              multiplied = 1;
              if (j == 3)
                outx[nz] += (z[id3 + 1] * z[id3 + 2] / dnrm2l(d3, z + id3) +
                             powl(det_r2 / det_u2, 0.25) * othor2[0] * othor2[1]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, z + id3 + 1);  // 1st row of C^-1/2
              if (j == 4)
                outx[nz] += (z[id3 + 2] * z[id3 + 2] / dnrm2l(d3, z + id3) +
                             powl(det_r2 / det_u2, 0.25) * othor2[1] * othor2[1]) *
                            Hx[p] / dnrm2sqrl(d3 - 1, z + id3 + 1);  // 2nd row
            }
          }  // end rows of H
          if (multiplied) {
            outi[nz++] = j + id5;
            multiplied = 0;
          }
        }  // end C^-1/2
      }  // end traverse all rows-block[5x5] of Pinv
    }  // end traverse all blocks[5x5] of Pinv
  }  // end traverse all cols of H

  outp[H->size1] = nz;
  cs_sprealloc(out_csc, 0);

  /*   out->storageType = H->storageType; */
  /*   numericsSparseMatrix(out)->csc = out_csc; */
  /*   out->size0 = (int)out->matrix2->csc->m; */
  /*   out->size1 = (int)out->matrix2->csc->n; */
  /*   numericsSparseMatrix(out)->origin = NSM_CSC; */

  free(x);
  free(z);
  free(othor1);
  free(othor2);
  return out;
}

#include "cs.h"
/* remove duplicate entries and zero entries from A */
static CS_INT cs_dupl_zeros(cs *A) {
  CS_INT i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w;
  CS_ENTRY *Ax;
  if (!CS_CSC(A)) return (0); /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  w = cs_malloc(m, sizeof(CS_INT));  /* get workspace */
  if (!w) return (0);                /* out of memory */
  for (i = 0; i < m; i++) w[i] = -1; /* row i not yet seen */
  for (j = 0; j < n; j++) {
    q = nz; /* column j will start at q */
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p]; /* A(i,j) is nonzero */
      // printf("p=%ld, Ai[p]=%ld, nz=%ld, i=%ld, q=%ld, w[i]=%ld, Ax[nz]=%f,
      // ",p,Ai[p],nz,i,q,w[i],Ax[p]);
      if (w[i] >= q || Ax[p] == 0) {
        Ax[w[i]] += Ax[p]; /* A(i,j) is a duplicate */
      } else {
        w[i] = nz;  /* record where row i occurs */
        Ai[nz] = i; /* keep A(i,j) */
        Ax[nz++] = Ax[p];
      }
      // printf("Ai[nz]=%ld\n", Ai [nz-1]);
    }
    Ap[j] = q; /* record start of column j */
  }
  Ap[n] = nz;                  /* finalize A */
  cs_free(w);                  /* free workspace */
  return (cs_sprealloc(A, 0)); /* remove extra space from A */
}

/* L = chol (A, [pinv parent cp]), pinv is optional */

static csn *cs_chol_2(const cs *A, const css *S, size_t iteration) {
  // CS_ENTRY d, lki, *Lx, *x, *Cx ;
  CS_ENTRY *Lx, *Cx;
  float_type d, lki, *x;
  CS_INT top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci;
  cs *L, *C, *E;
  csn *N;
  if (!CS_CSC(A) || !S || !S->cp || !S->parent) return (NULL);
  n = A->n;
  N = cs_calloc(1, sizeof(csn));        /* allocate result */
  c = cs_malloc(2 * n, sizeof(CS_INT)); /* get CS_INT workspace */
  // x = cs_malloc (n, sizeof (CS_ENTRY)) ;    /* get CS_ENTRY workspace */
  x = cs_malloc(n, sizeof(float_type)); /* get float_type workspace */
  cp = S->cp;
  pinv = S->pinv;
  parent = S->parent;
  C = pinv ? cs_symperm(A, pinv, 1) : ((cs *)A);
  // C = (cs *)A;
  E = pinv ? C : NULL; /* E is alias for A, or a copy E=A(p,p) */
  if (!N || !c || !x || !C) return (cs_ndone(N, E, c, x, 0));
  s = c + n;
  Cp = C->p;
  Ci = C->i;
  Cx = C->x;
  N->L = L = cs_spalloc(n, n, cp[n], 1, 0); /* allocate result */
  if (!L) return (cs_ndone(N, E, c, x, 0));
  Lp = L->p;
  Li = L->i;
  Lx = L->x;
  for (k = 0; k < n; k++) Lp[k] = c[k] = cp[k];
  for (k = 0; k < n; k++) /* compute L(k,:) for L*L' = C */
  {
    /* --- Nonzero pattern of L(k,:) ------------------------------------ */
    top = cs_ereach(C, k, parent, s, c); /* find pattern of L(k,:) */
    if (iteration == 8) printf("top = %lld\n", top);
    x[k] = 0.;                          /* x (0:k) is now zero */
    for (p = Cp[k]; p < Cp[k + 1]; p++) /* x = full(triu(C(:,k))) */
    {
      if (Ci[p] <= k) x[Ci[p]] = Cx[p];
      // if (Ci [p] <= k) x [Ci [p]] = (float_type)Cx [p] ;
    }
    d = x[k]; /* d = C(k,k) */
    if (iteration == 8) printf("d = C(k,k) = %5.40Le\n", d);
    x[k] = 0.; /* clear x for k+1st iteration */
    /* --- Triangular solve --------------------------------------------- */
    for (; top < n; top++) /* solve L(0:k-1,0:k-1) * x = C(:,k) */
    {
      i = s[top];             /* s [top..n-1] is pattern of L(k,:) */
      lki = x[i] / Lx[Lp[i]]; /* L(k,i) = x (i) / L(i,i) */
      x[i] = 0.;              /* clear x for k+1st iteration */
      for (p = Lp[i] + 1; p < c[i]; p++) {
        x[Li[p]] -= Lx[p] * lki;
      }
      // d -= lki * CS_CONJ (lki) ;            /* d = d - L(k,i)*L(k,i) */
      d -= lki * lki; /* d = d - L(k,i)*L(k,i) */
      p = c[i]++;
      Li[p] = k; /* store L(k,i) in column i */
      // Lx [p] = CS_CONJ (lki) ;
      Lx[p] = lki;
      if (iteration == 8) {
        printf(
            "k=%lld, i=%lld, p=%lld, Li[p]=%lld, Lp[p]=%lld, Lx [p]=%5.40e, lki=%5.40Le, "
            "d=%5.40Le\n",
            k, i, p, Li[p], Lp[p], Lx[p], lki, d);
      }
    }
    /* --- Compute L(k,k) ----------------------------------------------- */
    // if (CS_REAL (d) <= 0 || CS_IMAG (d) != 0)
    //   return (cs_ndone (N, E, c, x, 0)) ; /* not pos def */
    // if (CS_REAL (d) <= 0 || CS_IMAG (d) != 0)
    if (creall(d) <= 0 || cimagl(d) != 0) {
      printf("\n\n Factorized matrix has a negative eigenvalue at the cone i = %lld. \n\n",
             k / 5 + 1);
      return (cs_ndone(N, E, c, x, 0)); /* not pos def */
    }
    // if (d < 0.) d = 1e-40;
    p = c[k]++;
    Li[p] = k; /* store L(k,k) = sqrt (d) in column k */
    // Lx [p] = sqrt (d) ;
    Lx[p] = sqrtl(d);
    if (iteration == 8)
      printf("k=%lld, p2=%lld, Li[p]=%lld, Lp[p]=%lld, Lx [p]=%5.40e\n", k, p, Li[p], Lp[p],
             Lx[p]);
  }
  Lp[n] = cp[n]; /* finalize L */
  // if(iteration==16) printf("\n\n cs_chol_2 001  \n\n");
  return (cs_ndone(N, E, c, x, 1)); /* success: free E,s,x; return N */
}

static NumericsMatrix *compute_JQinv2Jt(const double *u1, const double *r1, const double *u2,
                                        const double *r2, const size_t vecSize,
                                        const size_t varsCount) {
  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3
  assert(d3 == 3);
  size_t d5 = d3 + 2;  // d5 = 5

  double *x = (double *)calloc(vecSize, sizeof(double));
  double *z = (double *)calloc(vecSize, sizeof(double));
  Nesterov_Todd_vector(3, u1, r1, vecSize, varsCount, x);  // x = pinv2_bar
  Nesterov_Todd_vector(3, u2, r2, vecSize, varsCount, z);  // z = pinv2_tilde

  NumericsMatrix *out = NM_create(NM_SPARSE, 5 * varsCount, 5 * varsCount);
  NM_triplet_alloc(out, (9 + 2 * (2 * 2)) * varsCount);
  CSparseMatrix *out_triplet = out->matrix2->triplet;

  float_type data = 0., nub = 0., nrb = 0., det_u = 0., det_r = 0.;

  int id3, id5;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * d3;
    id5 = i * d5;

    // Assign data for out[0,0]
    cs_entry(out_triplet, id5, id5, dnrm2sqrl(d3, x + id3) + dnrm2sqrl(d3, z + id3));

    // Assign data for out[0,1:2] & out[1:2, 0]
    for (size_t k = 1; k < d3; k++) {
      data = 2. * x[id3] * x[id3 + k];
      cs_entry(out_triplet, id5, id5 + k, data);
      cs_entry(out_triplet, id5 + k, id5, data);
    }

    // Assign data for out[0,3:4] & out[3:4, 0]
    for (size_t k = 1; k < d3; k++) {
      data = 2. * z[id3] * z[id3 + k];
      cs_entry(out_triplet, id5, id5 + 2 + k, data);
      cs_entry(out_triplet, id5 + 2 + k, id5, data);
    }

    // Assign data for matrix A = out[1:2,1:2]
    // Compute det(r1), det(u1)
    nrb = dnrm2l(d3 - 1, r1 + id3 + 1);
    nub = dnrm2l(d3 - 1, u1 + id3 + 1);
    // det_r = (r1[id3]+nrb)*(r1[id3]-nrb);
    det_r = r1[id3] - nrb;
    if (det_r <= 0.)
      det_r = (r1[id3] + nrb) * DBL_EPSILON;
    else
      det_r = (r1[id3] + nrb) * (r1[id3] - nrb);

    // det_u = (u1[id3]+nub)*(u1[id3]-nub);
    det_u = u1[id3] - nub;
    if (det_u < 0.)
      det_u = (u1[id3] + nub) * DBL_EPSILON;
    else
      det_u = (u1[id3] + nub) * (u1[id3] - nub);

    for (size_t k = 1; k < d3; k++) {
      for (size_t l = 1; l < d3; l++) {
        if (k == l)
          data = sqrtl(det_u / det_r) + 2. * x[id3 + k] * x[id3 + l];
        else
          data = 2. * x[id3 + k] * x[id3 + l];
        cs_entry(out_triplet, id5 + k, id5 + l, data);
      }
    }

    // Assign data for matrix C = out[3:4,3:4]
    // Compute det(r2), det(u2)
    nrb = dnrm2l(d3 - 1, r2 + id3 + 1);
    nub = dnrm2l(d3 - 1, u2 + id3 + 1);
    // det_r = (r2[id3]+nrb)*(r2[id3]-nrb);
    det_r = r2[id3] - nrb;
    if (det_r <= 0.)
      det_r = (r2[id3] + nrb) * DBL_EPSILON;
    else
      det_r = (r2[id3] + nrb) * (r2[id3] - nrb);

    // det_u = (u2[id3]+nub)*(u2[id3]-nub);
    det_u = u2[id3] - nub;
    if (det_u < 0.)
      det_u = (u2[id3] + nub) * DBL_EPSILON;
    else
      det_u = (u2[id3] + nub) * (u2[id3] - nub);

    for (size_t k = 1; k < d3; k++) {
      for (size_t l = 1; l < d3; l++) {
        if (k == l)
          data = sqrtl(det_u / det_r) + 2. * z[id3 + k] * z[id3 + l];
        else
          data = 2. * z[id3 + k] * z[id3 + l];
        cs_entry(out_triplet, id5 + k + 2, id5 + l + 2, data);
      }
    }
  }

  free(x);
  free(z);
  return out;
}

/* Return the matrix L^-1*H where L is Cholesky factor satisfying J*Q^-2*J' = L*L'. Using the
 * formula Qp. */
// static  NumericsMatrix *  multiply_LinvH(const double *u1, const double *r1, const double
// *u2, const double *r2, const size_t vecSize, const size_t varsCount, NumericsMatrix *H,
// CSparseMatrix **chol_L, FILE *file)
static NumericsMatrix *multiply_LinvH(const double *u1, const double *r1, const double *u2,
                                      const double *r2, const size_t vecSize,
                                      const size_t varsCount, NumericsMatrix *H,
                                      CSparseMatrix **chol_L, size_t iteration) {
  NumericsMatrix *LinvH = NM_new();
  NM_types storage = H->storageType;

  switch (storage) {
      /* case NM_DENSE: */
      /*   break; */
    case NM_SPARSE: {
      NumericsMatrix *JQJ =
          compute_JQinv2Jt(u1, r1, u2, r2, vecSize, varsCount);  // JQJ = J*Q^-2*J'
      assert(JQJ);

      CSparseMatrix *JQJ_csc = NM_csc(JQJ);

      CSparseMatrix *B = NM_csc(H);

      // Cholesky factor
      // css* S = cs_schol(1, JQJ_csc);
      CS_INT n = JQJ_csc->n;

      CS_INT cumsum = n;

      css *S = cs_calloc(1, sizeof(css));
      S->pinv = cs_malloc(n, sizeof(CS_INT));
      S->parent = cs_malloc(n, sizeof(CS_INT));
      S->cp = cs_malloc(n + 1, sizeof(CS_INT));

      cumsum = 5;
      for (int i = 0; i < n; i++) {
        if (cumsum <= 0) cumsum = 5;
        S->pinv[i] = i;
        S->parent[i] = i + 1;
        if (i == 0)
          S->cp[i] = 0;
        else {
          S->cp[i] = S->cp[i - 1] + cumsum;
          cumsum--;
        }
      }
      S->parent[n - 1] = -1;
      S->cp[n] = S->cp[n - 1] + 1;
      S->unz = S->lnz = S->cp[n];

      csn *N = cs_chol_2(JQJ_csc, S, iteration);  // L = N->L
      if (!N) {
        if (JQJ) NM_free(JQJ);
        if (S) cs_sfree(S);
        return NULL;
      }

      // X = L\B
      CSparseMatrix *X = cs_spalloc(B->m, B->n, B->nzmax, 1, 0); /* allocate result */

      CS_ENTRY *x, *b, *Xx, *Bx;
      CS_INT *xi, *pinv, top, k, i, p, *Bp, *Bi, *Xp, *Xi;

      x = cs_malloc(n, sizeof(CS_ENTRY));    /* get CS_ENTRY workspace */
      b = cs_malloc(n, sizeof(CS_ENTRY));    /* get CS_ENTRY workspace */
      xi = cs_malloc(2 * n, sizeof(CS_INT)); /* get CS_INT workspace */

      CS_INT xnz = 0;
      Xp = X->p;
      Xi = X->i;
      Xx = X->x;
      Xp[0] = 0;
      pinv = S->pinv;

      /* ---  X = L\B ---------------------------------------------- */
      for (k = 0; k < B->n; k++) {
        // printf("\n\n OK 001 c \n\n");
        /* permutation of the rows  of B(:,k) */
        // for(p = B->p [k] ; p < B->p [k+1] ; p++) x [B->i[p]] = B->x [p] ; /* scatter B  */
        // for(p = B->p [k] ; p < B->p [k+1] ; p++)
        // {
        //   CS_INT i_old= B->i[p];
        //   B->i[p] = pinv[i_old]; /* permute row indices with S->pinv */
        //   B->x[p] = x[i_old];
        // }
        /* call spsolve */
        // if(iteration==16) printf("\n\n multiply_LinvH 005 A \n\n");
        top = cs_spsolve(N->L, B, k, xi, x, NULL, 1); /* x = L\B(:,col) */
        // if(iteration==16) printf("\n\n multiply_LinvH 005 B \n\n");
        // printf("\n\n OK 001 d \n\n");
        // printf("\n x[:,%zu] = [", k);
        // for (int j=0; j<n; j++) printf("%f, ", x[j]);
        // printf("]; \n");
        /* store the result in X */
        if (Xp[k] + n - top > X->nzmax &&
            !cs_sprealloc(X, 2 * (X->nzmax) + n - top)) /* realloc X if need */
        {
          printf("\n\n multiply_LinvH, X = L\\B,  NULL!!! \n\n");
          if (x) free(x);
          if (b) free(b);
          if (xi) free(xi);
          if (JQJ) NM_free(JQJ);
          if (S) cs_sfree(S);
          if (N->U) cs_spfree(N->U);
          if (N->pinv) cs_free(N->pinv);
          if (N->B) cs_free(N->B);
          return NULL; /* (cs_done(X, w, x, 0)) ;   */ /* out of memory */
        }

        Xp = X->p;
        Xi = X->i;
        Xx = X->x;
        for (p = top; p < n; p++) {
          i = xi[p];   /* x(i) is nonzero */
          Xi[xnz] = i; /* store the result in X */
          Xx[xnz++] = x[i];
        }
        Xp[k + 1] = Xp[k] + n - top;
      }

      cs_sprealloc(X, 0);

      LinvH->storageType = H->storageType;
      numericsSparseMatrix(LinvH)->csc = X;
      LinvH->size0 = (int)LinvH->matrix2->csc->m;
      LinvH->size1 = (int)LinvH->matrix2->csc->n;
      numericsSparseMatrix(LinvH)->origin = NSM_CSC;

      *chol_L = N->L;  // for storage Cholesky factor

      if (x) free(x);
      if (b) free(b);
      if (xi) free(xi);
      if (JQJ) NM_free(JQJ);
      if (S) cs_sfree(S);
      if (N->U) cs_spfree(N->U);
      if (N->pinv) cs_free(N->pinv);
      if (N->B) cs_free(N->B);

      break;
    }

    default:
      fprintf(stderr,
              "Numerics, GRFC3D IPM, multiply_LinvH failed, unknown storage type for H.\n");
      exit(EXIT_FAILURE);
  }
  return LinvH;
}

/* Return the matrix L^-1*H where L is Cholesky factor satisfying J*Q^-2*J' = L*L'. Using the
 * formula Qp. */
// static  NumericsMatrix *  multiply_LinvH(const double *u1, const double *r1, const double
// *u2, const double *r2, const size_t vecSize, const size_t varsCount, NumericsMatrix *H,
// CSparseMatrix **chol_L, FILE *file)
static NumericsMatrix *multiply_UinvH(CSparseMatrix *chol_U, NumericsMatrix *H) {
  NumericsMatrix *UinvH = NM_new();

  NM_types storage = H->storageType;

  switch (storage) {
      /* case NM_DENSE: */
      /*   break; */
    case NM_SPARSE: {
      CSparseMatrix *B = NM_csc(H);
      CS_INT n = chol_U->n;

      // X = U\B
      CSparseMatrix *X = cs_spalloc(B->m, B->n, B->nzmax, 1, 0); /* allocate result */

      CS_ENTRY *x, *b, *Xx, *Bx;
      CS_INT *xi, *pinv, top, k, i, p, *Bp, *Bi, *Xp, *Xi;

      x = cs_malloc(n, sizeof(CS_ENTRY));    /* get CS_ENTRY workspace */
      b = cs_malloc(n, sizeof(CS_ENTRY));    /* get CS_ENTRY workspace */
      xi = cs_malloc(2 * n, sizeof(CS_INT)); /* get CS_INT workspace */

      CS_INT xnz = 0;
      Xp = X->p;
      Xi = X->i;
      Xx = X->x;
      Xp[0] = 0;

      /* ---  X = L\B ---------------------------------------------- */
      for (k = 0; k < B->n; k++) {
        top = cs_spsolve(chol_U, B, k, xi, x, NULL, 0); /* x = L\B(:,col) */

        /* store the result in X */
        if (Xp[k] + n - top > X->nzmax &&
            !cs_sprealloc(X, 2 * (X->nzmax) + n - top)) /* realloc X if need */
        {
          printf("\n\n multiply_UinvH, X = L\\B,  NULL!!! \n\n");
          if (x) free(x);
          if (b) free(b);
          if (xi) free(xi);
          return NULL; /* (cs_done(X, w, x, 0)) ;   */ /* out of memory */
        }

        Xp = X->p;
        Xi = X->i;
        Xx = X->x;
        for (p = top; p < n; p++) {
          i = xi[p];   /* x(i) is nonzero */
          Xi[xnz] = i; /* store the result in X */
          Xx[xnz++] = x[i];
        }
        Xp[k + 1] = Xp[k] + n - top;
      }

      cs_sprealloc(X, 0);

      UinvH->storageType = H->storageType;
      numericsSparseMatrix(UinvH)->csc = X;
      UinvH->size0 = (int)UinvH->matrix2->csc->m;
      UinvH->size1 = (int)UinvH->matrix2->csc->n;
      numericsSparseMatrix(UinvH)->origin = NSM_CSC;

      if (x) free(x);
      if (b) free(b);
      if (xi) free(xi);
      break;
    }

    default:
      fprintf(stderr,
              "Numerics, GRFC3D IPM, multiply_UinvH failed, unknown storage type for H.\n");
      exit(EXIT_FAILURE);
  }

  return UinvH;
}

static NumericsMatrix *compute_factor_U(const double *u1, const double *r1, const double *u2,
                                        const double *r2, const size_t vecSize,
                                        const size_t varsCount) {
  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3

  assert(d3 == 3);
  size_t d5 = d3 + 2;  // d5 = 5

  double *x = (double *)calloc(vecSize, sizeof(double));
  double *z = (double *)calloc(vecSize, sizeof(double));
  Nesterov_Todd_vector(3, u1, r1, vecSize, varsCount, x);  // x = pinv2_bar
  Nesterov_Todd_vector(3, u2, r2, vecSize, varsCount, z);  // z = pinv2_tilde

  NumericsMatrix *out = NM_create(NM_SPARSE, 5 * varsCount, 5 * varsCount);
  NM_triplet_alloc(out, 11 * varsCount);
  CSparseMatrix *out_triplet = out->matrix2->triplet;

  float_type p0 = 0., detx = 0., detz = 0., tmpx = 0., tmpz = 0., nx2 = 0., nz2 = 0.;
  float_type nxb = 0., nzb = 0.;

  int id3, id5;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * d3;
    id5 = i * d5;

    nxb = dnrm2l(2, x + id3 + 1);
    nzb = dnrm2l(2, z + id3 + 1);
    detx = (x[id3] - nxb) * (x[id3] + nxb);
    detz = (z[id3] - nzb) * (z[id3] + nzb);

    if (detx <= 0.) detx = fabsl(detx);
    if (detz <= 0.) detz = fabsl(detz);

    tmpx = x[id3] * x[id3] + x[id3 + 2] * x[id3 + 2] -
           x[id3 + 1] * x[id3 + 1];  // tmp = x0^2 - x1^2 + x2^2
    tmpz = z[id3] * z[id3] + z[id3 + 2] * z[id3 + 2] - z[id3 + 1] * z[id3 + 1];

    nx2 = dnrm2sqrl(d3, x + id3);
    nz2 = dnrm2sqrl(d3, z + id3);

    // out[0,0]
    p0 = sqrtl(detx * detx * nz2 + detz * detz * nx2) /
         (dnrm2l(d3, x + id3) * dnrm2l(d3, z + id3));
    cs_entry(out_triplet, id5, id5, p0);

    // out[0, 1:2]
    // if (detx <= 0.) detx = EPS;
    cs_entry(out_triplet, id5, id5 + 1, 2 * x[id3] * x[id3 + 1] * sqrtl(detx / (nx2 * tmpx)));
    cs_entry(out_triplet, id5, id5 + 2, 2 * x[id3] * x[id3 + 2] / sqrtl(tmpx));

    // out[1:2, 1:2]
    cs_entry(out_triplet, id5 + 1, id5 + 1, sqrtl(nx2 * detx / tmpx));
    cs_entry(out_triplet, id5 + 1, id5 + 2, 2 * x[id3 + 1] * x[id3 + 2] / sqrtl(tmpx));
    cs_entry(out_triplet, id5 + 2, id5 + 2, sqrtl(tmpx));

    // out[0, 3:4]
    // if (detz <= 0.) detz = EPS;
    cs_entry(out_triplet, id5, id5 + 3, 2 * z[id3] * z[id3 + 1] * sqrtl(detz / (nz2 * tmpz)));
    cs_entry(out_triplet, id5, id5 + 4, 2 * z[id3] * z[id3 + 2] / sqrtl(tmpz));

    // out[3:4, 3:4]
    cs_entry(out_triplet, id5 + 3, id5 + 3, sqrtl(nz2 * detz / tmpz));
    cs_entry(out_triplet, id5 + 3, id5 + 4, 2 * z[id3 + 1] * z[id3 + 2] / sqrtl(tmpz));
    cs_entry(out_triplet, id5 + 4, id5 + 4, sqrtl(tmpz));
  }
  free(x);
  free(z);
  return out;
}

/* Writing problem data under a Matlab format in a file  */
/* The data are printed under the form of a dense format */
/* problem: min .5 v'*M*v + f'*v, s.t. H*v + w \in K (Lorentz cone)
   d = dimenion
   n = number of contact points
   m = number of degrees of freedom
   M = m x m matrix
   f = m-vector
   H = n*d x m matrix
   w = n*d-vector */
static void printDataProbMatlabFile(NumericsMatrix *M, double *f, NumericsMatrix *H, double *w,
                                    int d, int n, int m, double *mu, double *mu_r,
                                    FILE *file) {
  // printf("\n\n printDataProbMatlabFile OK \n\n");
  fprintf(file, "d = %3i;\n", d);
  fprintf(file, "n = %6i;\n", n);
  fprintf(file, "m = %6i;\n", m);

  fprintf(file, "M = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(M), 0, file);
  fprintf(file, "];\n");
  // fprintf(file,"M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");
  fprintf(file, "M = sparse(M(:,1), M(:,2), M(:,3));\n");

  fprintf(file, "H = [");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, file);
  fprintf(file, "];\n");
  // fprintf(file,"H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");
  fprintf(file, "H = sparse(H(:,1), H(:,2), H(:,3));\n");

  fprintf(file, "f = [");
  for (int i = 0; i < m; i++) {
    fprintf(file, "%22.14e; ", -f[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "w = [");
  for (int i = 0; i < n * d; i++) {
    fprintf(file, "%22.14e; ", w[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "mu = [");
  for (int i = 0; i < n; i++) {
    fprintf(file, "%22.14e; ", mu[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "mu_r = [");
  for (int i = 0; i < n; i++) {
    fprintf(file, "%22.14e; ", mu_r[i]);
  }
  fprintf(file, "];\n");
}

/* print iteres under a Matlab format in a file */
/* iteration = index of the iteration
   v = global velocity
   u = velocity
   r = reaction
   d = dimenion
   n = number of contact points
   m = number of degrees of freedom
*/
static void printInteresProbMatlabFile(int iteration, double *v, double *u_1, double *u_2,
                                       double *r_1, double *r_2, int d, int n, int m,
                                       FILE *file) {
  fprintf(file, "v(%3i,:) = [", iteration + 1);
  for (int i = 0; i < m; i++) {
    fprintf(file, "%20.26e, ", v[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "u1(%3i,:) = [", iteration + 1);
  for (int i = 0; i < n * (d - 2); i++) {
    fprintf(file, "%20.26e, ", u_1[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "u2(%3i,:) = [", iteration + 1);
  for (int i = 0; i < n * (d - 2); i++) {
    fprintf(file, "%20.26e, ", u_2[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "r1(%3i,:) = [", iteration + 1);
  for (int i = 0; i < n * (d - 2); i++) {
    fprintf(file, "%20.26e, ", r_1[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "r2(%3i,:) = [", iteration + 1);
  for (int i = 0; i < n * (d - 2); i++) {
    fprintf(file, "%20.16e, ", r_2[i]);
  }
  fprintf(file, "];\n");

  return;
}

static void printInteresProbPythonFile(int iteration, double *v, double *u, double *u_1,
                                       double *u_2, double *r, double *r_1, double *r_2, int d,
                                       int n, int m, FILE *file) {
  int n_dminus2 = n * (d - 2);
  int nd = n * d;

  // v
  if (iteration == 0)
    fprintf(file, "v = np.array([");
  else
    fprintf(file, "x = np.array([");
  for (int i = 0; i < m; i++) {
    if (i == m - 1)
      fprintf(file, "%20.26e", v[i]);
    else
      fprintf(file, "%20.26e, ", v[i]);
  }
  fprintf(file, "])\n");
  if (iteration > 0) fprintf(file, "v = np.vstack((v, x))\n");

  // u
  if (iteration == 0)
    fprintf(file, "u = np.array([");
  else
    fprintf(file, "x = np.array([");
  for (int i = 0; i < nd; i++) {
    if (i == nd - 1)
      fprintf(file, "%20.26e", u[i]);
    else
      fprintf(file, "%20.26e, ", u[i]);
  }
  fprintf(file, "])\n");
  if (iteration > 0) fprintf(file, "u = np.vstack((u, x))\n");

  // u1
  if (iteration == 0)
    fprintf(file, "u1 = np.array([");
  else
    fprintf(file, "x = np.array([");
  for (int i = 0; i < n_dminus2; i++) {
    if (i == n_dminus2 - 1)
      fprintf(file, "%20.26e", u_1[i]);
    else
      fprintf(file, "%20.26e, ", u_1[i]);
  }
  fprintf(file, "])\n");
  if (iteration > 0) fprintf(file, "u1 = np.vstack((u1, x))\n");

  // u2
  if (iteration == 0)
    fprintf(file, "u2 = np.array([");
  else
    fprintf(file, "x = np.array([");
  for (int i = 0; i < n_dminus2; i++) {
    if (i == n_dminus2 - 1)
      fprintf(file, "%20.26e", u_2[i]);
    else
      fprintf(file, "%20.26e, ", u_2[i]);
  }
  fprintf(file, "])\n");
  if (iteration > 0) fprintf(file, "u2 = np.vstack((u2, x))\n");

  // r
  if (iteration == 0)
    fprintf(file, "r = np.array([");
  else
    fprintf(file, "x = np.array([");
  for (int i = 0; i < nd; i++) {
    if (i == nd - 1)
      fprintf(file, "%20.26e", r[i]);
    else
      fprintf(file, "%20.26e, ", r[i]);
  }
  fprintf(file, "])\n");
  if (iteration > 0) fprintf(file, "r = np.vstack((r, x))\n");

  // r1
  if (iteration == 0)
    fprintf(file, "r1 = np.array([");
  else
    fprintf(file, "x = np.array([");
  for (int i = 0; i < n_dminus2; i++) {
    if (i == n_dminus2 - 1)
      fprintf(file, "%20.26e", r_1[i]);
    else
      fprintf(file, "%20.26e, ", r_1[i]);
  }
  fprintf(file, "])\n");
  if (iteration > 0) fprintf(file, "r1 = np.vstack((r1, x))\n");

  // r2
  if (iteration == 0)
    fprintf(file, "r2 = np.array([");
  else
    fprintf(file, "x = np.array([");
  for (int i = 0; i < n_dminus2; i++) {
    if (i == n_dminus2 - 1)
      fprintf(file, "%20.26e", r_2[i]);
    else
      fprintf(file, "%20.26e, ", r_2[i]);
  }
  fprintf(file, "])\n");
  if (iteration > 0) fprintf(file, "r2 = np.vstack((r2, x))\n");

  return;
}

static void printVectorMatlabFile(int iteration, double *vec, int vecSize, FILE *file) {
  fprintf(file, "vector(%4i,:) = [", iteration + 1);
  for (int i = 0; i < vecSize; i++) {
    fprintf(file, "%24.16e, ", vec[i]);
  }
  fprintf(file, "];\n");
  return;
}

/* This function replaces for grfc3d_compute_error */
/* Compute:
    + Relative dual residual:
      ++ error_dual   = |Mv - H'r - f|/max{|Mv|, |Hr|, |f|}               if max >= tol
                      = |Mv - H'r - f|                                    otherwise
    + Relative primal residual:
      ++ error_primal = |u - Hv - w|/max{|H'v|, |w|, |u|}                 if max >= tol
                      = |u - Hv - w|                                      otherwise
    + Projection error
      For convex case
      ++ error_proj   = |r - projectionOnRollingCone(r-u)|/max{|r|, |u|}  if max >= tol
                      = |r - projectionOnRollingCone(r-u)|                otherwise

      For non-convex case
      ++ error_proj   = |r - projectionOnRollingCone(r-u-mu*|uT|-mur*||wR)|/max{|r|, |u|}  if
   max >= tol = |r - projectionOnRollingCone(r-u-mu*|uT|-mur*||wR)|                otherwise
*/
static void compute_errors(NumericsMatrix *M, NumericsMatrix *H, const double *w,
                           const double *f, double *r, double *u, double *v,
                           double *primalConstraint, double *pinfeas, double *dualConstraint,
                           double *dinfeas, double tolerance, double *proj_error,
                           double *full_error, int problemIsNotConvex) {
  /* Checks inputs */
  if (M == NULL || H == NULL || w == NULL || f == NULL || r == NULL || u == NULL || v == NULL)
    numerics_error("compute_errors", "null input");

  size_t nd = H->size0;
  size_t m = H->size1;
  size_t n = nd / 5;
  // printf("\n\nn = %zu, m = %zu\n\n", n, m);
  double norm_r = cblas_dnrm2(nd, r, 1);
  double norm_u = cblas_dnrm2(nd, u, 1);
  double max_val = 0.;

  double worktmp[5];

  /* --- Relative dual residual = |-Mv + Hr + f|/max{|Mv|, |Hr|, |f|} --- */
  double *HTr = (double *)calloc(m, sizeof(double));

  NM_gemv(1.0, M, v, 0.0, dualConstraint);  // dualConstraint = Mv
  max_val = cblas_dnrm2(m, dualConstraint, 1);
  cblas_daxpy(m, -1.0, f, 1, dualConstraint, 1);    // dualConstraint = Mv - f
  NM_tgemv(1.0, H, r, 0.0, HTr);                    // HTr = H'r
  cblas_daxpy(m, -1.0, HTr, 1, dualConstraint, 1);  // dualConstraint = Mv - f - H'r
  max_val = fmax(max_val, cblas_dnrm2(m, f, 1));
  max_val = fmax(max_val, cblas_dnrm2(m, HTr, 1));  // max_val = max{|Mv|, |f|, |H'r|, typr}

  *dinfeas = cblas_dnrm2(m, dualConstraint, 1);
  if (max_val >= tolerance) *dinfeas /= max_val;
  free(HTr);

  /* --- Relative primal residual = |H'v + w - u|/max{|H'v|, |w|, |u|} --- */
  max_val = 0.;
  NM_gemv(-1.0, H, v, 0.0, primalConstraint);  // primalConstraint = -Hv
  max_val = cblas_dnrm2(nd, primalConstraint, 1);
  cblas_daxpy(nd, -1.0, w, 1, primalConstraint, 1);  // primalConstraint = -Hv - w
  cblas_daxpy(nd, 1.0, u, 1, primalConstraint, 1);   // primalConstraint = u - Hv - w
  max_val = fmax(max_val, norm_u);
  max_val = fmax(max_val, cblas_dnrm2(nd, w, 1));  // max_val = max{|Hv|, |u|, |w|, typu}

  *pinfeas = cblas_dnrm2(nd, primalConstraint, 1);  // max_val = 1;
  if (max_val >= tolerance) *pinfeas /= max_val;

  /* --- Projection error = |r - projectionOnRollingCone(r-u)|/max{|r|, |u|} for convex case
   * --- */
  /* --- Projection error = |r - projectionOnRollingCone(r-u-mu*|uT|-mur*||wR)|/max{|r|, |u|}
   * for non-convex case --- */
  *proj_error = 0.;
  for (size_t i = 0; i < n; i++) {
    grfc3d_unitary_compute_and_add_error(&r[i * 5], &u[i * 5], 1., 1., proj_error, worktmp,
                                         problemIsNotConvex);
  }
  *proj_error = sqrt(*proj_error);

  max_val = fmax(norm_u, norm_r);  // max_val = 1;
  if (max_val >= tolerance) *proj_error /= max_val;

  /* --- Full error = Relative dual residual + Relative primal residual + Projection error ---
   */
  *full_error = *dinfeas + *pinfeas + *proj_error;
}

static void print_NAN_in_matrix(const NumericsMatrix *const m) {
  if (!m) {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  // printf("========== Numerics Matrix\n");
  // printf("========== size0 = %i, size1 = %i\n", m->size0, m->size1);

  switch (m->storageType) {
    case NM_DENSE: {
      // printf("========== storageType = NM_DENSE\n");
      break;
    }
    case NM_SPARSE_BLOCK: {
      assert(m->matrix1);
      // printf("========== storageType =  NM_SPARSE_BLOCK\n");
      break;
    }
    case NM_SPARSE: {
      assert(m->matrix2);
      // printf("========== storageType = NM_SPARSE\n");
      switch (m->matrix2->origin) {
        case NSM_TRIPLET: {
          // printf("========== origin =  NSM_TRIPLET\n");
          break;
        }
        case NSM_HALF_TRIPLET: {
          // printf("========== origin =  NSM_HALF_TRIPLET\n");
          break;
        }
        case NSM_CSC: {
          // printf("========== origin =  NSM_CSC\n");
          break;
        }
        case NSM_CSR: {
          // printf("========== origin =  NSM_CSR\n");
          break;
        }
        default: {
          // fprintf(stderr, "NM_display ::  unknown origin %d for sparse matrix\n",
          // m->matrix2->origin);
        }
      }

      // printf("========== size0 = %i, size1 = %i\n", m->size0, m->size1);
      CSparseMatrix *A;
      if (m->matrix2->triplet) {
        // printf("========== a matrix in format triplet is stored\n");
        A = m->matrix2->triplet;
      } else if (m->matrix2->csc) {
        // printf("========== a matrix in format csc is stored\n");
        A = m->matrix2->csc;
      } else if (m->matrix2->trans_csc) {
        // printf("========== a matrix in format trans_csc is stored\n");
        A = m->matrix2->trans_csc;
      }

      CS_INT p, nz, *Ap, *Ai;
      CS_ENTRY *Ax;

      Ap = A->p;
      Ai = A->i;
      Ax = A->x;
      nz = A->nz;

      for (p = 0; p < nz; p++) {
        if (Ax)
          if (isnan(Ax[p])) {
            printf("    %g %g : ", (double)(Ai[p]), (double)(Ap[p]));
            printf("%g\n", Ax[p]);
          }
      }
      break;
    }
    default: {
      fprintf(stderr, "display for NumericsMatrix: matrix type %d not supported!\n",
              m->storageType);
    }
  }
}

static void is_in_int_of_Lcone(const double *const x, const size_t vecSize,
                               const size_t varsCount) {
  size_t dim = vecSize / varsCount, id3 = 0;
  assert(dim == 3);
  float_type diffL = 0.;
  double diff = 0.;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * dim;
    // diff = (float_type)x[id3]-dnrm2l(dim-1, x+id3+1);
    diff = x[id3] - cblas_dnrm2(dim - 1, x + id3 + 1, 1);
    if (diff <= 0.) {
      printf("\n(double)     x0 = %9.65e", x[id3]);
      printf("\n(l doub)     x0 = %9.65Le", (float_type)x[id3]);
      printf("\n(double) x0+eps = %9.65e", x[id3] + DBL_EPSILON);
      printf("\n(double)|x_bar| = %9.65e", cblas_dnrm2(dim - 1, x + id3 + 1, 1));
      printf("\n(l doub)|x_bar| = %9.65Le\n", dnrm2l(dim - 1, x + id3 + 1));

      /*       printf("\ni = %zu: (double)x0 - (double)|x_bar| = (double) %9.65e\n",
       * i, */
      /*              diff); */

      diffL = (float_type)(x[id3] - cblas_dnrm2(dim - 1, x + id3 + 1, 1));
      printf("\ni = %zu: (double)x0 - (double)|x_bar| = (l doub) %9.65Le\n", i, diffL);

      diff = (double)(x[id3] - dnrm2l(dim - 1, x + id3 + 1));
      printf("\ni = %zu: (double)x0 - (l doub)|x_bar| = (double) %9.65e\n", i, diff);

      diffL = x[id3] - dnrm2l(dim - 1, x + id3 + 1);
      printf("\ni = %zu: (double)x0 - (l doub)|x_bar| = (l doub) %9.65Le\n", i, diffL);
    }
  }
}

static void update_w(double *w, double *w_origin, const double *velocity, const size_t vecSize,
                     const size_t varsCount, int update) {
  if (update == 0) return;

  size_t dim = vecSize / varsCount;
  assert(dim == 5);
  for (size_t i = 0; i < vecSize; ++i) {
    if (i % dim == 0) /* w[i] = w_origin[i]/(problem->mu[(int)(i/d)]) */
      w[i] = w_origin[i] +
             sqrt(velocity[i + 1] * velocity[i + 1] + velocity[i + 2] * velocity[i + 2]) +
             sqrt(velocity[i + 3] * velocity[i + 3] + velocity[i + 4] * velocity[i + 4]);
  }
}

static double compute_min_steplenght_of4(const double *x, const double *dx, const double *y,
                                         const double *dy, const double *z, const double *dz,
                                         const double *t, const double *dt,
                                         const size_t vecSize, const size_t varsCount,
                                         double gamma) {
  double alpha_primal_1 = grfc3d_getStepLength(x, dx, vecSize, varsCount, gamma);
  double alpha_primal_2 = grfc3d_getStepLength(y, dy, vecSize, varsCount, gamma);
  double alpha_dual_1 = grfc3d_getStepLength(z, dz, vecSize, varsCount, gamma);
  double alpha_dual_2 = grfc3d_getStepLength(t, dt, vecSize, varsCount, gamma);

  return fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
}

/** Compute a block matrix J of form
 *      |  a0    b1    b2                c0    d1    d2                   ...  |
 *      |                                                                      |
 *      |  b1    A11   A12               0     0     0                    ...  |
 *      |                                                                      |
 *      |  b2    A21   A22               0     0     0                    ...  |
 *      |                                                                      |
 *      |  0     0     0                 d1    C11   C12                  ...  |
 *      |                                                                      |
 *      |  0     0     0                 d2    C21   C22                  ...  |
 *  J = |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      | ...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ...  |
 */

static NumericsMatrix *compute_JQinv(const double *u1, const double *r1, const double *u2,
                                     const double *r2, const size_t vecSize,
                                     const size_t varsCount) {
  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3

  assert(d3 == 3);
  size_t d5 = d3 + 2;  // d5 = 5

  double *x = (double *)calloc(vecSize, sizeof(double));
  double *z = (double *)calloc(vecSize, sizeof(double));
  Nesterov_Todd_vector(1, u1, r1, vecSize, varsCount, x);  // x = pinv_bar
  Nesterov_Todd_vector(1, u2, r2, vecSize, varsCount, z);  // z = pinv_tilde

  NumericsMatrix *out = NM_create(NM_SPARSE, 5 * varsCount, 6 * varsCount);
  NM_triplet_alloc(out, (10 + 2 * (2 * 2)) * varsCount);
  CSparseMatrix *out_triplet = out->matrix2->triplet;

  float_type data = 0., nub = 0., nrb = 0., det_u = 0., det_r = 0.;

  size_t id3, id5;
  size_t nx3 = 3 * varsCount;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * d3;
    id5 = i * d5;

    // a0
    cs_entry(out_triplet, id5, id3, dnrm2sqrl(d3, x + id3));

    // b
    for (size_t k = 1; k < d3; k++) {
      data = 2. * x[id3] * x[id3 + k];
      cs_entry(out_triplet, id5, id3 + k, data);
      cs_entry(out_triplet, id5 + k, id3, data);
    }

    // A
    // det(r1), det(u1)
    nrb = dnrm2l(d3 - 1, r1 + id3 + 1);
    nub = dnrm2l(d3 - 1, u1 + id3 + 1);
    // det_r = (r1[id3]+nrb)*(r1[id3]-nrb);
    det_r = r1[id3] - nrb;
    if (det_r <= 0.)
      det_r = (r1[id3] + nrb) * DBL_EPSILON;
    else
      det_r = (r1[id3] + nrb) * (r1[id3] - nrb);

    // det_u = (u1[id3]+nub)*(u1[id3]-nub);
    det_u = u1[id3] - nub;
    if (det_u < 0.)
      det_u = (u1[id3] + nub) * DBL_EPSILON;
    else
      det_u = (u1[id3] + nub) * (u1[id3] - nub);

    for (size_t k = 1; k < d3; k++) {
      for (size_t l = 1; l < d3; l++) {
        if (k == l)
          data = powl(det_u / det_r, 0.25) + 2. * x[id3 + k] * x[id3 + l];
        else
          data = 2. * x[id3 + k] * x[id3 + l];
        cs_entry(out_triplet, id5 + k, id3 + l, data);
      }
    }

    // c0
    cs_entry(out_triplet, id5, id3 + nx3, dnrm2sqrl(d3, z + id3));

    // d
    for (size_t k = 1; k < d3; k++) {
      data = 2. * z[id3] * z[id3 + k];
      cs_entry(out_triplet, id5, id3 + k + nx3, data);
      cs_entry(out_triplet, id5 + k + 2, id3 + nx3, data);
    }

    // C
    // det(r2), det(u2)
    nrb = dnrm2l(d3 - 1, r2 + id3 + 1);
    nub = dnrm2l(d3 - 1, u2 + id3 + 1);
    // det_r = (r2[id3]+nrb)*(r2[id3]-nrb);
    det_r = r2[id3] - nrb;
    if (det_r <= 0.)
      det_r = (r2[id3] + nrb) * DBL_EPSILON;
    else
      det_r = (r2[id3] + nrb) * (r2[id3] - nrb);

    // det_u = (u2[id3]+nub)*(u2[id3]-nub);
    det_u = u2[id3] - nub;
    if (det_u < 0.)
      det_u = (u2[id3] + nub) * DBL_EPSILON;
    else
      det_u = (u2[id3] + nub) * (u2[id3] - nub);

    for (size_t k = 1; k < d3; k++) {
      for (size_t l = 1; l < d3; l++) {
        if (k == l)
          data = powl(det_u / det_r, 0.25) + 2. * z[id3 + k] * z[id3 + l];
        else
          data = 2. * z[id3 + k] * z[id3 + l];
        cs_entry(out_triplet, id5 + k + 2, id3 + l + nx3, data);
      }
    }
  }

  free(x);
  free(z);
  return out;
}

static void print_neg_eigval(const double *x, const size_t vecSize, const size_t varsCount) {
  size_t d3 = (size_t)(vecSize / varsCount);  // d3 = 3
  assert(d3 == 3);

  size_t id3 = 0;
  double min_eigval = 0.;
  for (size_t i = 0; i < varsCount; i++) {
    id3 = i * d3;
    min_eigval = x[id3] - cblas_dnrm2(d3 - 1, x + id3 + 1, 1);
    if (min_eigval < 0.) printf("Cone %zu: min_eigval = %e\n", i, min_eigval);
  }
}

static double detMat(NumericsMatrix *A) {
  double det = 1.;

  CSparseMatrix *A_csc = NM_csc(A);

  css *S;
  csn *N;
  CS_INT n;
  if (!CS_CSC(A_csc)) return (0); /* check inputs */
  n = A_csc->n;
  S = cs_sqr(1, A_csc, 0);          /* ordering and symbolic analysis */
  N = cs_lu(A_csc, S, DBL_EPSILON); /* numeric LU factorization */

  cs *L = N->L;
  CS_INT p, j, *Lp, *Li;
  CS_ENTRY *Lx;
  if (!CS_CSC(L)) return (0); /* check inputs */
  n = L->n;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;
  for (j = 0; j < n; j++) {
    det *= Lx[Lp[j]];
  }

  cs *U = N->U;
  CS_INT *Up, *Ui;
  CS_ENTRY *Ux;
  if (!CS_CSC(U)) return (0); /* check inputs */
  n = U->n;
  Up = U->p;
  Ui = U->i;
  Ux = U->x;
  for (j = n - 1; j >= 0; j--) {
    det *= Ux[Up[j + 1] - 1];
  }

  cs_sfree(S);
  cs_nfree(N);

  return det;
}

/* --------------------------- Interior-point method implementation
 * ------------------------------ */
/* General rolling friction contact problem */
/* Convex case: */
/* problem: min .5 v'*M*v + f'*v, s.t. H*v + w \in F (rolling friction cone)
   d = dimenion
   n = number of contact points
   m = number of degrees of freedom
   M = m x m matrix
   f = m-vector
   H = n*d x m matrix
   w = n*d-vector */

void grfc3d_IPM(GlobalRollingFrictionContactProblem *restrict problem,
                double *restrict reaction, double *restrict velocity,
                double *restrict globalVelocity, int *restrict info,
                SolverOptions *restrict options) {
  clock_t t1 = clock();

  printf("\n\n#################### grfc3d_IPM is starting ####################\n\n");

  /* -------------------------- Variable declaration -------------------------- */
  // the size of the problem detection
  size_t m = problem->M->size0;
  size_t nd = problem->H->size1;
  size_t d = problem->dimension;

  size_t n = problem->numberOfContacts;  // n = 1; nd = 5;
  size_t m_plus_nd = m + nd;
  size_t d_minus_2 = d - 2;
  size_t d_plus_1 = d + 1;
  size_t n_dminus2 = n * d_minus_2;
  size_t n_dplus1 = n * d_plus_1;
  size_t pos_t = 0;

  size_t id3 = 0;  // id3 = i*d_minus_2 used for the loop of cones
  size_t id5 = 0;  // id5 = i*d         used for the loop of cones

  NumericsMatrix *M = NULL, *minus_M = NULL, *Minv = NULL, *HMinv = NULL, *HMinvHt = NULL;
  NumericsMatrix *H_origin = NULL, *minus_H = NULL, *minus_Ht = NULL, *Ht = NULL;

  /* symmetrization of the matrix M */
  if (!(NM_is_symmetric(problem->M))) {
    printf("#################### SYMMETRIZATION ####################\n");
    NumericsMatrix *MT = NM_transpose(problem->M);
    problem->M = NM_add(1 / 2., problem->M, 1 / 2., MT);
    NM_free(MT);
  }

  /* if SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE =
     SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE, we force the copy into a NM_SPARSE
     storageType */
  DEBUG_PRINTF("problem->M->storageType : %i\n", problem->M->storageType);
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] ==
          SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE &&
      problem->M->storageType == NM_SPARSE_BLOCK) {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    printf(
        "\n\n\n######################### FORCE SPARSE STORAGE "
        "#########################\n\n\n");
    M = NM_create(NM_SPARSE, problem->M->size0, problem->M->size1);
    NM_copy_to_sparse(problem->M, M, DBL_EPSILON);
  } else {
    M = problem->M;
  }
  size_t M_nzmax = NM_nnz(M);
  int block_number_of_M = M->size0 / 3;
  unsigned int *blocksizes_of_M = NULL;

  DEBUG_PRINTF("problem->M->storageType : %i\n", problem->H->storageType);
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] ==
          SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE &&
      problem->H->storageType == NM_SPARSE_BLOCK) {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    H_origin = NM_create(NM_SPARSE, problem->H->size1, problem->H->size0);
    NM_copy_to_sparse(NM_transpose(problem->H), H_origin, DBL_EPSILON);
  } else {
    H_origin = NM_transpose(
        problem->H);  // H <== H' because some different storages in fclib and paper
  }

  printf("nnz H = %8zu density = %9.4f, \t mu = %.2e\n", NM_nnz(problem->H),
         NM_nnz(problem->H) / 1.0 / nd / m, problem->mu[0]);

  // initialize solver if it is not set
  int internal_allocation = 0;
  if (!options->dWork) {
    grfc3d_IPM_init(problem, options);
    internal_allocation = 1;
  }

  Grfc3d_IPM_data *data = (Grfc3d_IPM_data *)options->solverData;
  NumericsMatrix *P_mu = data->P_mu->mat;

  double *w_origin = problem->b;
  double *f = problem->q;

  double alpha_primal_1 = 0.0;
  double alpha_primal_2 = 0.0;
  double alpha_dual_1 = 0.0;
  double alpha_dual_2 = 0.0;
  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;
  double sigma_mu = 0.0;

  // Copy the initial data into executive vars
  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  size_t no_n = 0, no_m = 0, no_nd = 0, no_ndm2 = 0;

  double *t = data->tmp_vault_n[no_n++];
  double *t_prime = data->tmp_vault_n[no_n++];
  for (size_t i = 0; i < n; ++i) {
    t[i] = 2.0;
    t_prime[i] = 1.0;
  }
  double *velocity_1 = data->tmp_vault_n_dminus2[no_ndm2++];  // = (t, u_bar)
  double *velocity_2 = data->tmp_vault_n_dminus2[no_ndm2++];  // = (t', u_tilde)
  double *reaction_1 = data->tmp_vault_n_dminus2[no_ndm2++];  // = (r0, r_bar)
  double *reaction_2 = data->tmp_vault_n_dminus2[no_ndm2++];  // = (r0, r_tilde)

  // For Newton directions
  double *d_globalVelocity = data->tmp_vault_m[no_m++];

  double *d_velocity = data->tmp_vault_nd[no_nd++];
  double *d_reaction = data->tmp_vault_nd[no_nd++];

  double *d_velocity_1 = data->tmp_vault_n_dminus2[no_ndm2++];
  double *d_velocity_2 = data->tmp_vault_n_dminus2[no_ndm2++];
  double *d_reaction_1 = data->tmp_vault_n_dminus2[no_ndm2++];
  double *d_reaction_2 = data->tmp_vault_n_dminus2[no_ndm2++];

  double *d_t = data->tmp_vault_n[no_n++];
  double *d_t_prime = data->tmp_vault_n[no_n++];

  double *rhs = options->dWork;

  double *rhs_tmp = NULL, *rhs_tmp2 = NULL;

  double *velocity_tmp1 = NULL, *velocity_tmp2 = NULL, *reaction_tmp1 = NULL,
         *reaction_tmp2 = NULL;
  double *d_velocity_tmp1 = NULL, *d_velocity_tmp2 = NULL, *d_reaction_tmp1 = NULL,
         *d_reaction_tmp2 = NULL;

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  size_t max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp0 = 0.999;
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  size_t iteration = 0;
  double pinfeas = 1e300;
  double pinfeas_new = 1e300;
  double dinfeas = 1e300;
  double complem_1 = 1e300;
  double complem_2 = 1e300;
  double gapVal = 1e300;
  double relgap = 1e300;

  double u1dotr1 = 1e300;  // u1 = velecity_1, r1 = reaction_1
  double u2dotr2 = 1e300;  // u2 = velecity_2, r2 = reaction_2

  double udotr = 1e300;
  double proj_error = 1e300;  // projection error

  double residu_LS1_m = 0.0, residu_LS2_m = 0.0;
  double residu_LS1_nd = 0.0, residu_LS2_nd = 0.0;
  double residu_LS1_ndplus1 = 0.0, residu_LS2_ndplus1 = 0.0;

  long blocks_nzmax = 3 * 2 * n;  // for 3x3 no scaling

  NumericsMatrix *Jac = NULL, *Jactmp = NULL; /* Jacobian matrix */

  long Jac_nzmax;
  int jacobian_is_nan = 0;

  NumericsMatrix *J = compute_J_matrix(n); /* use for Jac */
  NumericsMatrix *Jt = NULL;

  double full_error = 1e300;
  int nRefine = 0;            // = info of refinement solver, otherwise refine is empty
  double residu_refine = 0.;  // = residu of refinement solver

  double *w = data->tmp_vault_nd[no_nd++];
  NM_gemv(1.0, P_mu, w_origin, 0.0, w);  // w_origin --> w

  double gmm = gmmp0;
  double barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  // compute -f
  // cblas_dscal(m, -1.0, f, 1); // f <== -f because some different storages in fclib and paper
  // double *minus_f = (double*)calloc(m, sizeof(double));
  // cblas_dcopy(m, f, 1, minus_f, 1);
  // cblas_dscal(m, -1.0, minus_f, 1);

  // change of variable
  // H_origin --> H
  NumericsMatrix *H = NM_multiply(P_mu, H_origin);
  size_t H_nzmax = NM_nnz(H);

  /* -------------------------- Declaration -------------------------- */
  // For 3x3 no scaling
  NumericsMatrix *block_1 = NULL, *block_2 = NULL;
  NumericsMatrix *arrowMat_u1 = NULL, *arrowMat_u2 = NULL, *arrowMat_r1 = NULL,
                 *arrowMat_r2 = NULL;
  NumericsMatrix *Z = NULL, *ZJT = NULL;

  // For NT scaling
  NumericsMatrix *Qp_bar = NULL, *Qp_tilde = NULL;
  NumericsMatrix *Qp2_bar = NULL, *Qp2_tilde = NULL;
  NumericsMatrix *Qpinv_bar = NULL, *Qpinv_tilde = NULL;
  NumericsMatrix *Qpinv2_bar = NULL, *Qpinv2_tilde = NULL;
  NumericsMatrix *Qinv = NULL, *Qinv2 = NULL;
  NumericsMatrix *JQinv = NULL, *JQinvT = NULL;
  NumericsMatrix *JQinv2 = NULL;
  NumericsMatrix *JQJ = NULL;
  NumericsMatrix *P_inv = NULL, *P_invT = NULL, *HMHP = NULL, *PHMHP = NULL;
  NumericsMatrix *P_inv_F = NULL;
  NumericsMatrix *PinvH = NULL, *PinvH_T = NULL;
  NumericsMatrix *identity = NULL;

  NumericsMatrix *chol_U = NULL;
  CSparseMatrix *chol_L = NULL, *chol_U_csc = NULL, *chol_UT_csc = NULL;

  double *p_bar = NULL, *p_tilde = NULL, *pinv_tilde = NULL, *pinv2_bar = NULL;
  double *p2_bar = NULL, *p2_tilde = NULL, *pinv_bar = NULL, *pinv2_tilde = NULL;

  double *d_velocity_1_hat = NULL, *d_velocity_2_hat = NULL;
  double *d_reaction_1_check = NULL, *d_reaction_2_check = NULL;
  double *velocity_1_inv = NULL, *velocity_2_inv = NULL;
  double *velocity_1_hat = NULL, *velocity_2_hat = NULL;
  double *velocity_1_hat_inv = NULL, *velocity_2_hat_inv = NULL;
  double *velocity_1_hat_inv_dvhat_drcheck_1 = NULL,
         *velocity_2_hat_inv_dvhat_drcheck_2 = NULL;
  double *velocity_hat_inv_dvhat_drcheck = NULL;

  double *tmp_n_dplus1 = NULL, *tmp_nd = NULL, *tmp_nd2 = NULL;

  double *Qinv2x_bar = NULL, *Qinv2x_tilde = NULL;

  double *Hvw = NULL;
  double *rhs_save = NULL, *rhs_dx = NULL;
  double *iden = NULL;

  double *Hrf = NULL, *HMHrfw = NULL, *rdr = NULL, *MfHrdr = NULL;

  /* -------------------------- Allocation -------------------------- */
  // For residuals
  double *dualConstraint = data->tmp_vault_m[no_m++];

  double *primalConstraint = data->tmp_vault_nd[no_nd++];

  // For predictor step
  double *v_plus_dv = data->tmp_vault_nd[no_nd++];  // v_plus_dv = velocity +
                                                    // alpha_primal * d_velocity
  double *r_plus_dr = data->tmp_vault_nd[no_nd++];  // r_plus_dr = reaction +
                                                    // alpha_primal * d_reaction

  double *complemConstraint_1 = data->tmp_vault_n_dminus2[no_ndm2++];
  double *complemConstraint_2 = data->tmp_vault_n_dminus2[no_ndm2++];

  // For RHS
  double *dvdr_jprod_1 = data->tmp_vault_n_dminus2[no_ndm2++];  // dvdr_jprod_1 = dv_1 o dr_1
  double *dvdr_jprod_2 = data->tmp_vault_n_dminus2[no_ndm2++];  // dvdr_jprod_2 = dv_2 o dr_2

  double *tmp_n_dminus2_1 = data->tmp_vault_n_dminus2[no_ndm2++];
  double *tmp_n_dminus2_2 = data->tmp_vault_n_dminus2[no_ndm2++];

  switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
      rhs_save = (double *)calloc(m + nd + n_dplus1, sizeof(double));
      rhs_dx = (double *)calloc(m + nd + n_dplus1, sizeof(double));

      minus_H = NM_create(H->storageType, H->size0, H->size1);  // -H
      NM_copy(H, minus_H);
      NM_scal(-1.0, minus_H);

      minus_Ht = NM_transpose(minus_H);  // -H'
      Jt = NM_transpose(J);              // J'
      break;

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
      p_bar = data->tmp_vault_n_dminus2[no_ndm2++];
      p_tilde = data->tmp_vault_n_dminus2[no_ndm2++];
      p2_bar = data->tmp_vault_n_dminus2[no_ndm2++];
      p2_tilde = data->tmp_vault_n_dminus2[no_ndm2++];

      d_velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_1_check = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_2_check = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_inv = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_inv = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_1_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat_inv_dvhat_drcheck_1 = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv_dvhat_drcheck_2 = data->tmp_vault_n_dminus2[no_ndm2++];

      rhs_save = (double *)calloc(m + nd + n_dplus1, sizeof(double));  // for printing

      minus_H = NM_create(H->storageType, H->size0, H->size1);  // -H
      NM_copy(H, minus_H);
      NM_scal(-1.0, minus_H);

      minus_Ht = NM_transpose(minus_H);  // -H'

      Jt = NM_transpose(J);  // J'
      break;

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv:
      velocity_hat_inv_dvhat_drcheck = data->tmp_vault_n_dplus1[0];

      Hvw = data->tmp_vault_nd[no_nd++];

      d_velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_1_check = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_2_check = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_1_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat_inv_dvhat_drcheck_1 = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv_dvhat_drcheck_2 = data->tmp_vault_n_dminus2[no_ndm2++];

      rhs_save = (double *)calloc(m + nd + n_dplus1, sizeof(double));  // for printing

      minus_H = NM_create(H->storageType, H->size0, H->size1);  // -H
      NM_copy(H, minus_H);
      NM_scal(-1.0, minus_H);

      minus_Ht = NM_transpose(minus_H);  // -H'
      identity = NM_eye(n_dplus1);
      break;

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ:
      velocity_hat_inv_dvhat_drcheck = data->tmp_vault_n_dplus1[0];

      Hvw = data->tmp_vault_nd[no_nd++];

      d_velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_1_check = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_2_check = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_1_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat_inv_dvhat_drcheck_1 = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv_dvhat_drcheck_2 = data->tmp_vault_n_dminus2[no_ndm2++];

      Qinv2x_bar = data->tmp_vault_n_dminus2[no_ndm2++];
      Qinv2x_tilde = data->tmp_vault_n_dminus2[no_ndm2++];

      rhs_save = (double *)calloc(m + nd, sizeof(double));  // for printing

      Ht = NM_transpose(H);  // H'
      minus_M = NM_create(M->storageType, M->size0, M->size1);
      NM_copy(M, minus_M);
      NM_scal(-1.0, minus_M);
      break;

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH:
      velocity_hat_inv_dvhat_drcheck = data->tmp_vault_n_dplus1[0];

      Hvw = data->tmp_vault_nd[no_nd++];
      tmp_nd = data->tmp_vault_nd[no_nd++];

      d_velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_1_check = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_2_check = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_1_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat_inv_dvhat_drcheck_1 = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv_dvhat_drcheck_2 = data->tmp_vault_n_dminus2[no_ndm2++];

      Qinv2x_bar = data->tmp_vault_n_dminus2[no_ndm2++];
      Qinv2x_tilde = data->tmp_vault_n_dminus2[no_ndm2++];

      rhs_save = (double *)calloc(m + nd, sizeof(double));  // for printing

      minus_M = NM_create(M->storageType, M->size0, M->size1);
      NM_copy(M, minus_M);
      NM_scal(-1.0, minus_M);
      identity = NM_eye(nd);
      break;

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ:
      Hrf = data->tmp_vault_m[no_m++];
      HMHrfw = data->tmp_vault_nd[no_nd++];
      rdr = data->tmp_vault_nd[no_nd++];
      MfHrdr = data->tmp_vault_m[no_m++];

      velocity_hat_inv_dvhat_drcheck = data->tmp_vault_n_dplus1[0];

      d_velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_1_check = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_2_check = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_1_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat_inv_dvhat_drcheck_1 = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv_dvhat_drcheck_2 = data->tmp_vault_n_dminus2[no_ndm2++];

      Qinv2x_bar = data->tmp_vault_n_dminus2[no_ndm2++];
      Qinv2x_tilde = data->tmp_vault_n_dminus2[no_ndm2++];

      blocksizes_of_M = (unsigned int *)malloc(block_number_of_M * sizeof(unsigned int));
      for (int i = 0; i < block_number_of_M; i++) *(blocksizes_of_M + i) = 3;
      Minv = NM_inverse_diagonal_block_matrix(M, block_number_of_M, blocksizes_of_M);
      free(blocksizes_of_M);
      blocksizes_of_M = NULL;

      rhs_save = (double *)calloc(nd, sizeof(double));  // for printing
      Ht = NM_transpose(H);                             // H'
      HMinv = NM_multiply(H, Minv);
      HMinvHt = NM_multiply(HMinv, Ht);

      break;

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
      Hrf = data->tmp_vault_m[no_m++];
      HMHrfw = data->tmp_vault_nd[no_nd++];
      // PHMHrfw = data->tmp_vault_nd[no_nd++];
      tmp_nd = data->tmp_vault_nd[no_nd++];
      rdr = data->tmp_vault_nd[no_nd++];
      MfHrdr = data->tmp_vault_m[no_m++];

      velocity_hat_inv_dvhat_drcheck = data->tmp_vault_n_dplus1[0];

      d_velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_1_check = data->tmp_vault_n_dminus2[no_ndm2++];
      d_reaction_2_check = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_1_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv = data->tmp_vault_n_dminus2[no_ndm2++];

      velocity_1_hat_inv_dvhat_drcheck_1 = data->tmp_vault_n_dminus2[no_ndm2++];
      velocity_2_hat_inv_dvhat_drcheck_2 = data->tmp_vault_n_dminus2[no_ndm2++];

      Qinv2x_bar = data->tmp_vault_n_dminus2[no_ndm2++];
      Qinv2x_tilde = data->tmp_vault_n_dminus2[no_ndm2++];

      blocksizes_of_M = (unsigned int *)malloc(block_number_of_M * sizeof(unsigned int));
      for (int i = 0; i < block_number_of_M; i++) *(blocksizes_of_M + i) = 3;
      Minv = NM_inverse_diagonal_block_matrix(M, block_number_of_M, blocksizes_of_M);
      free(blocksizes_of_M);
      blocksizes_of_M = NULL;

      rhs_save = (double *)calloc(nd, sizeof(double));  // for printing
      Ht = NM_transpose(H);                             // H'
      HMinv = NM_multiply(H, Minv);
      HMinvHt = NM_multiply(HMinv, Ht);

      identity = NM_eye(nd);
      break;

    default:
      printf("Some vars are not allocated.\n");
  }

  /* -------------------------- Print into matlab file -------------------------- */
  FILE *iterates, *iterates_python;
  FILE *matrixH;
  char matlab_name[100];
  sprintf(matlab_name, "iteratesNC%zu.m", n);

  char python_name[100];
  sprintf(python_name, "PrimitiveSoup-nc-%zu.py", n);

  int reinit = 0;
  while (1) {
    if (reinit == 1)  // reset solver only once
    {
      reinit = 0;
      iteration = 0;
      gmm = gmmp0;
      hasNotConverged = 1;

      /* 1. v */
      for (size_t i = 0; i < m; ++i) {
        globalVelocity[i] = 0.01;
        d_globalVelocity[i] = 0.0;  // reset d_globalVelocity
      }

      /* 2. & 3. u, r */
      for (size_t i = 0; i < nd; ++i) {
        velocity[i] = 0.001;
        if (i % d == 0) velocity[i] = 3.0;
        d_velocity[i] = 0.0;  // reset d_velocity

        reaction[i] = 0.04;
        if (i % d == 0) reaction[i] = 0.5;
        d_reaction[i] = 0.0;  // reset d_reaction

        // reset Hvw
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] ==
                SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ ||
            options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] ==
                SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH)
          Hvw[i] = 0.;
      }

      /* 4. t,t' */
      for (size_t i = 0; i < n; ++i) {
        t[i] = 2.0;
        t_prime[i] = 1.0;
        d_t[i] = 0.0;
        d_t_prime[i] = 0.0;  // reset d_t, d_t_prime
      }
    }

    /* writing data in a Matlab file */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
      // iterates = fopen("iterates.m", "w");
      iterates = fopen(matlab_name, "w");

      printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, problem->mu_r, iterates);
    }

    // ComputeErrorGlobalRollingPtr computeError = NULL;
    // computeError = (ComputeErrorGlobalRollingPtr)&grfc3d_compute_error;

    /* -------------------------- Display problem info -------------------------- */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] ==
        SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES) {
      numerics_printf_verbose(1, "---- GRFC3D - IPM - Problem information");
      numerics_printf_verbose(1, "---- GRFC3D - IPM - 1-norm of M = %g norm of f = %g ",
                              NM_norm_1(M), norm_f);
      numerics_printf_verbose(1, "---- GRFC3D - IPM - inf-norm of M = %g ", NM_norm_inf(M));

      numerics_printf_verbose(1, "---- GRFC3D - IPM - 1-norm of H = %g norm of w = %g ",
                              NM_norm_1(problem->H), norm_w);
      numerics_printf_verbose(1, "---- GRFC3D - IPM - inf-norm of H = %g ",
                              NM_norm_inf(problem->H));
      numerics_printf_verbose(1, "---- GRFC3D - IPM - M is symmetric = %i ",
                              NM_is_symmetric(M));

      numerics_printf_verbose(1, "---- GRFC3D - IPM - M size = (%i, %i) ", M->size0, M->size1);
      numerics_printf_verbose(1, "---- GRFC3D - IPM - H size = (%i, %i) ", problem->H->size0,
                              problem->H->size1);
    }

    /* -------------------------- Display IPM iterations -------------------------- */
    // numerics_printf_verbose(-1, "problem dimensions d, n, m: %1i, %6i, %6i\n",d, n, m);
    numerics_printf_verbose(-1, "problem dimensions n, nd x m: %1i, %6i x %-6i", n, nd, m);
    switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL: {
        numerics_printf_verbose(-1, "LS solution: 3x3 no scaling\n");
        numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = %d\n",
                                options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT]);
        break;
      }
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2: {
        numerics_printf_verbose(-1, "LS solution: 3x3 NT scaling with Qp2");
        numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = %d\n",
                                options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT]);
        break;
      }
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv: {
        numerics_printf_verbose(-1, "LS solution: 3x3 NT scaling with J*Qp^-1");
        numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = %d\n",
                                options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT]);
        break;
      }
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ: {
        numerics_printf_verbose(-1, "LS solution: 2x2 NT scaling with J*Qp^-2*J'");
        numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = %d\n",
                                options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT]);
        break;
      }
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH: {
        numerics_printf_verbose(-1, "LS solution: 2x2 NT scaling with P^-1*H");
        numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = %d\n",
                                options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT]);
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES)
          numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY = %d\n",
                                  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY]);
        break;
      }
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ: {
        numerics_printf_verbose(-1, "LS solution: 1x1 JQJ NT scaling");
        numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = %d\n",
                                options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT]);
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES)
          numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY = %d\n",
                                  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY]);
        break;
      }
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH: {
        numerics_printf_verbose(-1, "LS solution: 1x1 QPH NT scaling");
        numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = %d\n",
                                options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT]);
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES)
          numerics_printf_verbose(-1, "SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY = %d\n",
                                  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY]);
        break;
      }
      default: {
        printf("ERROR\n");
      }
    }
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
        SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
      numerics_printf_verbose(
          -1,
          "| it  No| rel gap | pinfeas | dinfeas |  <u, r> | proj err| complem1| complem2| "
          "full err| barparam| alpha_p | alpha_d | sigma | |dv|/|v| | |du|/|u| | |dr|/|r| | "
          "1st residu m nd n(d+1) | 2nd residu m nd n(d+1) |resi refine|");
      numerics_printf_verbose(
          -1,
          "-----------------------------------------------------------------------------------"
          "-----------------------------------------------------------------------------------"
          "------------------------------------------------------");
    } else {
      // numerics_printf_verbose(-1, "| it| rel gap | pinfeas | dinfeas |  <u, r> | proj err|
      // complem1| complem2| full err| barparam| alpha_p | alpha_d | sigma | |dv|/|v| |
      // |du|/|u| | |dr|/|r| |residu m|residu nd|  n(d+1) |"); numerics_printf_verbose(-1, "|
      // it| rel gap | pinfeas | dinfeas |  <u, r> | proj err| complem1| complem2| full err|
      // barparam| alpha_p | alpha_d | sigma | |dv|/|v| | |du|/|u| | |dr|/|r| | 1st residu m nd
      // n(d+1) | 2nd residu m nd n(d+1) |"); numerics_printf_verbose(-1,
      // "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
      numerics_printf_verbose(-1,
                              "| it| rel gap | pinfeas | dinfeas |  <u, r> | proj err| "
                              "complem1| complem2| barparam|  alpha  | sigma | |dv|/|v| | "
                              "|du|/|u| | |dr|/|r| | residu   m  nd  n(d+1) |");
      numerics_printf_verbose(
          -1,
          "-----------------------------------------------------------------------------------"
          "-----------------------------------------------------------------------------");
    }

    /* -------------------------- Check the full criterion -------------------------- */
    while (iteration < max_iter) {
      /* -------------------------- Extract vectors -------------------------- */
      /* 2. velocity_1 = (t, u_bar), velocity_2 = (t_prime, u_tilde) */
      extract_vector(velocity, nd, n, 2, 3, velocity_1);
      extract_vector(velocity, nd, n, 4, 5, velocity_2);
      for (size_t i = 0; i < n; i++) {
        id3 = i * d_minus_2;
        velocity_1[id3] = t[i];
        velocity_2[id3] = t_prime[i];
      }

      /* 3. reaction_1 = (r0, r_bar), reaction_2 = (r0, r_tilde) */
      extract_vector(reaction, nd, n, 2, 3, reaction_1);
      extract_vector(reaction, nd, n, 4, 5, reaction_2);

      /* writing data in a Matlab file */
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
        printInteresProbMatlabFile(iteration, globalVelocity, velocity_1, velocity_2,
                                   reaction_1, reaction_2, d, n, m, iterates);
      }

      /* writing solution data in a Python file */
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE]) {
        if (iteration == 0) iterates_python = fopen(python_name, "w");
        printInteresProbPythonFile(iteration, globalVelocity, velocity, velocity_1, velocity_2,
                                   reaction, reaction_1, reaction_2, d, n, m, iterates_python);
      }

      // Update w
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1 &&
          iteration > 0)  // & (totalresidual <= 100*tol))
      {
        // printf("update w\n");
        for (unsigned int i = 0; i < nd; ++i) {
          if (i % d == 0) {
            w[i] =
                w_origin[i] +
                sqrt(velocity[i + 1] * velocity[i + 1] + velocity[i + 2] * velocity[i + 2]) +
                sqrt(velocity[i + 3] * velocity[i + 3] + velocity[i + 4] * velocity[i + 4]);
          }
        }
        /* for(unsigned int i = 0; i < nd; ++ i) printf("%20.14e\n",w[i]); */
      }
      // primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas_new, tol);
      // // primalConstraint = u - H*v -w

      // // Compute again primalConstraint with new w
      // NM_gemv(-1.0, H, globalVelocity, 0.0, primalConstraint);  // primalConstraint = -Hv
      // double max_val = cblas_dnrm2(nd, primalConstraint, 1);
      // cblas_daxpy(nd, -1.0, w, 1, primalConstraint, 1);         // primalConstraint = -Hv -
      // w cblas_daxpy(nd, 1.0, velocity, 1, primalConstraint, 1); // primalConstraint = u - Hv
      // - w

      // /* Dual gap = (primal value - dual value)/ (1 + abs(primal value) + abs(dual value))
      // */
      // // Note: primal objectif func = 1/2 * v' * M *v + f' * v
      // dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);

      /* Gap value = u'.v */
      gapVal = cblas_ddot(nd, reaction, 1, velocity, 1);
      // gapVal = cblas_ddot(n_dminus2, reaction_1, 1, velocity_1, 1) + cblas_ddot(n_dminus2,
      // reaction_1, 1, velocity_1, 1);

      // Note: primal objectif func = 1/2 * v' * M *v + f' * v
      relgap = grfc3d_relGap(M, f, w, globalVelocity, reaction, nd, m, gapVal);

      barr_param = gapVal / n;
      // barr_param = grfc3d_complemResidualNorm(velocity, reaction, nd, n);
      // barr_param = fabs(barr_param);

      complem_1 = grfc3d_complemResidualNorm(velocity_1, reaction_1, n_dminus2,
                                             n);  // (t, u_bar) o (r0, r_bar)
      complem_2 = grfc3d_complemResidualNorm(velocity_2, reaction_2, n_dminus2,
                                             n);  // (t', u_tilde) o (r0, r_tilde)

      udotr = gapVal;

      if (udotr < 0.)
        udotr =
            1e300;  // To avoid the negative value. Normally, udotr must always be positive.

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)  // non-smooth case
      {
        compute_errors(M, H, w, f, reaction, velocity, globalVelocity, primalConstraint,
                       &pinfeas, dualConstraint, &dinfeas, tol, &proj_error, &full_error, 1);
      }

      else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 0)  // convex case
      {
        compute_errors(M, H, w, f, reaction, velocity, globalVelocity, primalConstraint,
                       &pinfeas, dualConstraint, &dinfeas, tol, &proj_error, &full_error, 0);
      }

      // if (fmax(pinfeas, fmax(dinfeas, fmin(udotr, proj_error))) <= tol)
      if (fmax(pinfeas, fmax(dinfeas, udotr)) <= tol) {
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
          numerics_printf_verbose(
              -1, "| %2i %3d| %7.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
              iteration, nRefine, relgap, pinfeas, dinfeas, udotr, proj_error, complem_1,
              complem_2, full_error, barr_param);
        } else {
          // numerics_printf_verbose(-1, "| %2i| %7.1e | %.1e | %.1e | %.1e | %.1e | %.1e |
          // %.1e | %.1e | %.1e |",
          //                         iteration, relgap, pinfeas, dinfeas, udotr, proj_error,
          //                         complem_1, complem_2, full_error, barr_param);
          numerics_printf_verbose(
              -1, "| %2i| %7.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |", iteration,
              relgap, pinfeas, dinfeas, udotr, proj_error, complem_1, complem_2, barr_param);
        }

        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
          printInteresProbMatlabFile(iteration, globalVelocity, velocity_1, velocity_2,
                                     reaction_1, reaction_2, d, n, m, iterates);

        hasNotConverged = 0;
        break;
      }

      switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
        case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL: {
          /* -------------------------- FIRST linear system -------------------------- */
          /* 1. Build the Jacobian matrix
           *
           *               m     nd      n(d-2)    n(d-2)
           *            |  M     -H'        0         0    | m
           *            |                                  |
           *    Jac   = | -H      0      [...... J ......] | nd
           *            |                                  |
           *            |  0   block_1   arw(r1)      0    | n(d-2)
           *            |                                  |
           *            |  0   block_2      0      arw(r2) | n(d-2)
           *
           *  where
           *
           *                1       2        2
           *            |   t     u_bar'     0  | 1
           *  block_1 = |                       |         x n blocks
           *            | u_bar     tI       0  | 2
           *
           *  and
           *
           *                1       2        2
           *            |   t'      0     u_tilde' | 1
           *  block_2 = |                          |      x n blocks
           *            | u_tilde   0       t'I    | 2
           *
           *
           * ANOTHER VERSION
           *               m     nd     n(d+1)
           *            |  M     -H'     0  | m
           *            |                   |
           *    Jac   = | -H      0      J  | nd
           *            |                   |
           *            |  0     ZJ'     R  | n(d+1)
           *
           *  where
           *  Z = diag(arw(t ,u_bar), arw(t',u_tilde)),
           *  R = diag(arw(r0,r_bar), arw(r0,r_tilde))
           *
           */
          Jac = NM_create(NM_SPARSE, m + nd + n_dplus1, m + nd + n_dplus1);
          Jac_nzmax = M_nzmax + 2 * H_nzmax + 2 * 9 * n * n + 6 * n;
          NM_triplet_alloc(Jac, Jac_nzmax);
          Jac->matrix2->origin = NSM_TRIPLET;

          // if(!identity) identity = NM_eye(nd);
          // NM_scal(-barr_param, identity);

          NM_insert(Jac, M, 0, 0);
          NM_insert(Jac, minus_Ht, 0, m);
          NM_insert(Jac, minus_H, m, 0);
          // NM_insert(Jac, identity, m, m);
          NM_insert(Jac, J, m, m_plus_nd);

          arrowMat_u1 = Arrow_repr(velocity_1, n_dminus2, n);
          arrowMat_u2 = Arrow_repr(velocity_2, n_dminus2, n);
          Z = NM_create(NM_SPARSE, n_dplus1, n_dplus1);
          NM_triplet_alloc(Z, 7 * 2 * n);
          Z->matrix2->origin = NSM_TRIPLET;
          NM_insert(Z, arrowMat_u1, 0, 0);
          NM_insert(Z, arrowMat_u2, n_dminus2, n_dminus2);
          ZJT = NM_multiply(Z, Jt);

          arrowMat_r1 = Arrow_repr(reaction_1, n_dminus2, n);
          arrowMat_r2 = Arrow_repr(reaction_2, n_dminus2, n);

          /* Add the 3rd row into Jac matrix */
          NM_insert(Jac, ZJT, m_plus_nd, m);
          NM_insert(Jac, arrowMat_r1, m_plus_nd, m_plus_nd);
          NM_insert(Jac, arrowMat_r2, m_plus_nd + n_dminus2, m_plus_nd + n_dminus2);

          if (block_1) block_1 = NM_free(block_1);
          if (block_2) block_2 = NM_free(block_2);
          if (arrowMat_u1) arrowMat_u1 = NM_free(arrowMat_u1);
          if (arrowMat_u2) arrowMat_u2 = NM_free(arrowMat_u2);
          if (arrowMat_r1) arrowMat_r1 = NM_free(arrowMat_r1);
          if (arrowMat_r2) arrowMat_r2 = NM_free(arrowMat_r2);
          if (identity) identity = NM_free(identity);

          /* Correction of w to take into account the dependence on the tangential velocity */
          update_w(w, w_origin, velocity, nd, d,
                   options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]);

          /*  2. Build the right-hand-side
           *
           *  rhs = -     <==== ATTENTION to negative sign
           *        [ M*v - H'*r + f ]  m         dualConstraint
           *        [  u - H*v - w   ]  nd        primalConstraint
           *        [   u_1 o r_1    ]  n(d-2)    complemConstraint 1
           *        [   u_2 o r_2    ]  n(d-2)    complemConstraint 2
           */
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, primalConstraint, 1, rhs + m, 1);

          JA_prod(velocity_1, reaction_1, n_dminus2, n,
                  complemConstraint_1);  // complemConstraint_1 = u1 o r1
          JA_prod(velocity_2, reaction_2, n_dminus2, n,
                  complemConstraint_2);  // complemConstraint_2 = u2 o r2
          cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs + m_plus_nd, 1);
          assert((m_plus_nd + 2 * n_dminus2) ==
                 (m_plus_nd + n_dplus1));  // m + nd + n(d+1): size of rhs
          cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs + m_plus_nd + n_dminus2, 1);

          cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
          cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_save, 1);

          /* 3. Solve non-symmetric Newton system without NT scaling via LU factorization */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          NM_LU_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS1_m = dnrm2l(m, rhs_save);
          residu_LS1_nd = dnrm2l(nd, rhs_save + m);
          residu_LS1_ndplus1 = dnrm2l(n_dplus1, rhs_save + m + nd);

          /* 4. Retrieve the directions for PREDICTOR step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          cblas_dcopy(n_dminus2, rhs + m_plus_nd, 1, d_velocity_1, 1);
          assert((m_plus_nd + 2 * n_dminus2) ==
                 (m_plus_nd + n_dplus1));  // m + nd + n(d+1): size of rhs
          cblas_dcopy(n_dminus2, rhs + m_plus_nd + n_dminus2, 1, d_velocity_2, 1);

          for (size_t i = 0; i < n; i++) {
            id5 = i * d;
            id3 = i * d_minus_2;
            d_velocity[id5] = d_velocity_1[id3] + d_velocity_2[id3];
            d_velocity[id5 + 1] = d_velocity_1[id3 + 1];
            d_velocity[id5 + 2] = d_velocity_1[id3 + 2];
            d_velocity[id5 + 3] = d_velocity_2[id3 + 1];
            d_velocity[id5 + 4] = d_velocity_2[id3 + 2];
          }

          /* 5. Computing the affine step-length */

          alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
          alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
          alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
          alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

          alpha_primal =
              fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
          alpha_dual = alpha_primal;

          /* updating the gamma parameter used to compute the step-length */
          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* -------------------------- Predictor step of Mehrotra --------------------------
           */
          cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
          cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
          cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          /* affine barrier parameter */
          barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n;

          /* computing the centralization parameter */
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * alpha_primal * alpha_primal) : sgmp3;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e)) / d;
          // sigma = fmin(1.0, pow(barr_param_a / barr_param, e));

          /* -------------------------- SECOND linear system -------------------------- */
          // 6. Update the RHS
          // 1st term: u_1 o r_1
          // already assigned in complemConstraint_1 & complemConstraint_2

          // 2nd term: 2 * barr_param * sigma * e
          iden = JA_iden(n_dminus2, n);
          cblas_dscal(n_dminus2, 2 * barr_param * sigma, iden, 1);  // iden = 2nd term

          // 3rd term: du_1 o dr_1
          JA_prod(d_velocity_1, d_reaction_1, n_dminus2, n, dvdr_jprod_1);
          JA_prod(d_velocity_2, d_reaction_2, n_dminus2, n, dvdr_jprod_2);

          // Update complemConstraint = 1st term - 2nd term + 3rd term
          NV_sub(complemConstraint_1, iden, n_dminus2, tmp_n_dminus2_1);
          NV_sub(complemConstraint_2, iden, n_dminus2, tmp_n_dminus2_2);
          if (iden) {
            free(iden);
            iden = NULL;
          }
          NV_add(tmp_n_dminus2_1, dvdr_jprod_1, n_dminus2,
                 complemConstraint_1);  // complemConstraint_1 = updated rhs_1
          NV_add(tmp_n_dminus2_2, dvdr_jprod_2, n_dminus2,
                 complemConstraint_2);  // complemConstraint_1 = updated rhs_2

          // Update rhs
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, primalConstraint, 1, rhs + m, 1);
          cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs + m_plus_nd, 1);
          cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs + m_plus_nd + n_dminus2, 1);

          cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
          cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_save, 1);

          /* 7. Solve the 2nd linear system */
          // print_NAN_in_matrix(Jac);
          // if (NV_isnan(rhs, m + nd + n_dplus1)) printf("(2nd sys) NaN in RHS before solving,
          // i = %zu\n", iteration);

          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
              SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
            nRefine = NM_LU_refine(Jac, rhs, 1e-10, max_iter, &residu_refine);
          } else {
            NM_LU_solve(Jac, rhs, 1);
          }

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS2_m = dnrm2l(m, rhs_save);
          residu_LS2_nd = dnrm2l(nd, rhs_save + m);
          residu_LS2_ndplus1 = dnrm2l(n_dplus1, rhs_save + m + nd);

          /* 8. Retrieve the directions for CORRECTOR step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          cblas_dcopy(n_dminus2, rhs + m_plus_nd, 1, d_velocity_1, 1);
          cblas_dcopy(n_dminus2, rhs + m_plus_nd + n_dminus2, 1, d_velocity_2, 1);

          for (size_t i = 0; i < n; i++) {
            id5 = i * d;
            id3 = i * d_minus_2;
            d_velocity[id5] = d_velocity_1[id3] + d_velocity_2[id3];
            d_velocity[id5 + 1] = d_velocity_1[id3 + 1];
            d_velocity[id5 + 2] = d_velocity_1[id3 + 2];
            d_velocity[id5 + 3] = d_velocity_2[id3 + 1];
            d_velocity[id5 + 4] = d_velocity_2[id3 + 2];
            d_t[i] = d_velocity_1[id3];
            d_t_prime[i] = d_velocity_2[id3];
          }

          break;
        }  // end of SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL

        case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2: {
          /* -------------------------- Compute NT directions -------------------------- */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            Qp_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp2_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp2_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            NM_triplet_alloc(Qp_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp2_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp2_tilde, d_minus_2 * d_minus_2 * n);

            family_of_F(velocity_1, reaction_1, n_dminus2, n, NULL, NULL, Qp_bar, Qpinv_bar,
                        Qp2_bar, NULL);
            family_of_F(velocity_2, reaction_2, n_dminus2, n, NULL, NULL, Qp_tilde,
                        Qpinv_tilde, Qp2_tilde, NULL);

          } else if (options->iparam
                         [SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                     SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            Nesterov_Todd_vector(0, velocity_1, reaction_1, n_dminus2, n, p_bar);
            Qp_bar = QRmat(p_bar, n_dminus2, n);
            Nesterov_Todd_vector(2, velocity_1, reaction_1, n_dminus2, n, p2_bar);
            Qp2_bar = QRmat(p2_bar, n_dminus2, n);

            Nesterov_Todd_vector(0, velocity_2, reaction_2, n_dminus2, n, p_tilde);
            Qp_tilde = QRmat(p_tilde, n_dminus2, n);
            Nesterov_Todd_vector(2, velocity_2, reaction_2, n_dminus2, n, p2_tilde);
            Qp2_tilde = QRmat(p2_tilde, n_dminus2, n);
          }

          /* -------------------------- FIRST linear system -------------------------- */
          /*  1. Build the Jacobian matrix
           *
           *           m    nd   n(d+1)
           *        |  M    -H'    0   | m
           *        |                  |
           *  Jac = | -H     0     J   | nd
           *        |                  |
           *        |  0     J'    Q^2 | n(d+1)
           *
           *          n(d-2)    n(d-2)
           *        | Qp_bar      0    | n(d-2)
           *   Q  = |                  |
           *        |  0      Qp_tilde | n(d-2)
           */

          Jac = NM_create(NM_SPARSE, m + nd + n_dplus1, m + nd + n_dplus1);
          Jac_nzmax = M_nzmax + 2 * H_nzmax + 2 * 9 * n * n + 6 * n;
          NM_triplet_alloc(Jac, Jac_nzmax);
          Jac->matrix2->origin = NSM_TRIPLET;
          NM_insert(Jac, M, 0, 0);
          NM_insert(Jac, minus_Ht, 0, m);
          NM_insert(Jac, minus_H, m, 0);
          NM_insert(Jac, J, m, m_plus_nd);
          NM_insert(Jac, Jt, m_plus_nd, m);
          NM_insert(Jac, Qp2_bar, m_plus_nd, m_plus_nd);
          NM_insert(Jac, Qp2_tilde, m_plus_nd + n_dminus2, m_plus_nd + n_dminus2);

          /* Correction of w to take into account the dependence on the tangential velocity */
          update_w(w, w_origin, velocity, nd, d,
                   options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]);

          /*  2. Build the right-hand-side
           *
           *  rhs = -     <==== ATTENTION to negative sign
           *        [ M*v - H'*r + f ]  m         dualConstraint        = a (No have negative
           * sign at the beginning) [  u - H*v - w   ]  nd        primalConstraint      = b [
           * r_1       ]  n(d-2)    complemConstraint 1 | [      r_2       ]  n(d-2)
           * complemConstraint 2 | = c
           */
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, primalConstraint, 1, rhs + m, 1);

          cblas_dcopy(n_dminus2, reaction_1, 1, rhs + m_plus_nd, 1);
          assert((m_plus_nd + 2 * n_dminus2) ==
                 (m_plus_nd + n_dplus1));  // m + nd + n(d+1): size of rhs
          cblas_dcopy(n_dminus2, reaction_2, 1, rhs + m_plus_nd + n_dminus2, 1);

          cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
          cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_save, 1);

          /* 3. Solve full symmetric Newton system with NT scaling via LDLT factorization */

          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          NSM_linearSolverParams(Jac)->solver = NSM_HSL;
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS1_m = dnrm2l(m, rhs_save);
          residu_LS1_nd = dnrm2l(nd, rhs_save + m);
          residu_LS1_ndplus1 = dnrm2l(n_dplus1, rhs_save + m + nd);

          /* 4. Retrieve the directions for PREDICTOR step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          cblas_dcopy(n_dminus2, rhs + m_plus_nd, 1, d_velocity_1, 1);
          cblas_dcopy(n_dminus2, rhs + m_plus_nd + n_dminus2, 1, d_velocity_2, 1);

          for (size_t i = 0; i < n; i++) {
            id5 = i * d;
            id3 = i * d_minus_2;
            d_velocity[id5] = d_velocity_1[id3] + d_velocity_2[id3];
            d_velocity[id5 + 1] = d_velocity_1[id3 + 1];
            d_velocity[id5 + 2] = d_velocity_1[id3 + 2];
            d_velocity[id5 + 3] = d_velocity_2[id3 + 1];
            d_velocity[id5 + 4] = d_velocity_2[id3 + 2];
          }

          /* 5. Computing the affine step-length */
          alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
          alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
          alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
          alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

          alpha_primal =
              fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
          alpha_dual = alpha_primal;

          /* updating the gamma parameter used to compute the step-length */
          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* -------------------------- Predictor step of Mehrotra --------------------------
           */
          cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
          cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
          cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          /* affine barrier parameter */
          barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / (n);

          /* computing the centralization parameter */
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * alpha_primal * alpha_primal) : sgmp3;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e)) / d;

          /* -------------------------- SECOND linear system -------------------------- */
          // 6. Update the RHS
          // 1st term: r_1
          // already assigned in reaction_1 & reaction_2

          /* 2nd term: Qp_bar * { (Qp_bar * u_1)^-1 o [(Qp_bar * du_1) o (Qpinv_1 * dr_1)] } */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            NM_gemv(1.0, Qp_bar, velocity_1, 0.0,
                    velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            NM_gemv(1.0, Qp_tilde, velocity_2, 0.0,
                    velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            NM_gemv(1.0, Qp_bar, d_velocity_1, 0.0,
                    d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            NM_gemv(1.0, Qp_tilde, d_velocity_2, 0.0,
                    d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2

            NM_gemv(1.0, Qpinv_bar, d_reaction_1, 0.0,
                    d_reaction_1_check);  // d_reaction_1_check     = Qpinv_bar * dr_1
            NM_gemv(1.0, Qpinv_tilde, d_reaction_2, 0.0,
                    d_reaction_2_check);  // d_reaction_2_check     = Qpinv_tilde * dr_2

            JA_prod(d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                    dvdr_jprod_1);  // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(
                d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                dvdr_jprod_2);  // dvdr_jprod_2      = (Qp_tilde * du_1) o (Qpinv_tilde * dr_1)

            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // = (Qp_bar * u_1)^-1 o (Qp_bar *
                                                          // du_1) o (Qpinv_bar * dr_1)
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o (Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_1)

            NM_gemv(1.0, Qp_bar, velocity_1_hat_inv_dvhat_drcheck_1, 0.0,
                    complemConstraint_1);  // complemConstraint_1 = 2nd term_1
            NM_gemv(1.0, Qp_tilde, velocity_2_hat_inv_dvhat_drcheck_2, 0.0,
                    complemConstraint_2);  // complemConstraint_2 = 2nd term_2

          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpz(velocity_1, reaction_1, velocity_1, n_dminus2, n,
                  velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            QNTpz(velocity_2, reaction_2, velocity_2, n_dminus2, n,
                  velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            QNTpz(velocity_1, reaction_1, d_velocity_1, n_dminus2, n,
                  d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            QNTpz(velocity_2, reaction_2, d_velocity_2, n_dminus2, n,
                  d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2

            QNTpinvz(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                     d_reaction_1_check);  // d_reaction_1_check     = Qpinv_bar * dr_1
            QNTpinvz(velocity_2, reaction_2, d_reaction_2, n_dminus2, n,
                     d_reaction_2_check);  // d_reaction_2_check     = Qpinv_tilde * dr_2

            JA_prod(d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                    dvdr_jprod_1);  // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(
                d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                dvdr_jprod_2);  // dvdr_jprod_2      = (Qp_tilde * du_1) o (Qpinv_tilde * dr_1)

            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // = (Qp_bar * u_1)^-1 o (Qp_bar *
                                                          // du_1) o (Qpinv_bar * dr_1)
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o (Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_1)

            QNTpz(velocity_1, reaction_1, velocity_1_hat_inv_dvhat_drcheck_1, n_dminus2, n,
                  complemConstraint_1);  // complemConstraint_1 = 2nd term_1
            QNTpz(velocity_2, reaction_2, velocity_2_hat_inv_dvhat_drcheck_2, n_dminus2, n,
                  complemConstraint_2);  // complemConstraint_2 = 2nd term_2
          }

          /* 3rd term: 2 * barr_param * sigma * (u_1)^-1  */
          Jinv(velocity_1, n_dminus2, n, velocity_1_inv);  // velocity_1_inv    = u_1^-1
          Jinv(velocity_2, n_dminus2, n, velocity_2_inv);  // velocity_2_inv    = u_2^-1
          cblas_dscal(n_dminus2, 2 * barr_param * sigma, velocity_1_inv,
                      1);  // velocity_1_inv = 3rd term
          cblas_dscal(n_dminus2, 2 * barr_param * sigma, velocity_2_inv,
                      1);  // velocity_2_inv = 3rd term

          /* Update complemConstraint = r_1 + 2nd term - 3rd term */
          cblas_daxpy(n_dminus2, 1.0, reaction_1, 1, complemConstraint_1,
                      1);  // complemConstraint_1 = r_1 + 2nd term_1
          cblas_daxpy(n_dminus2, 1.0, reaction_2, 1, complemConstraint_2, 1);
          cblas_daxpy(n_dminus2, -1.0, velocity_1_inv, 1, complemConstraint_1,
                      1);  // complemConstraint_1 = r_1 + 2nd term_1 - - 3rd term
          cblas_daxpy(n_dminus2, -1.0, velocity_2_inv, 1, complemConstraint_2, 1);

          // Update rhs
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, primalConstraint, 1, rhs + m, 1);
          cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs + m_plus_nd, 1);
          cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs + m_plus_nd + n_dminus2, 1);

          cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
          cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_save, 1);

          /* 7. Solve the 2nd linear system */
          // print_NAN_in_matrix(Jac);
          // if (NV_isnan(rhs, m + nd + n_dplus1)) printf("(2nd sys) NaN in RHS before solving,
          // i = %zu\n", iteration);

          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
              SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
            cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_save, 1);
            NM_LDLT_refine(J, rhs, rhs_save, 1, tol, 10, 0);
          } else
            NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS2_m = dnrm2l(m, rhs_save);
          residu_LS2_nd = dnrm2l(nd, rhs_save + m);
          residu_LS2_ndplus1 = dnrm2l(n_dplus1, rhs_save + m + nd);

          /* 8. Retrieve the directions for CORRECTOR step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          cblas_dcopy(n_dminus2, rhs + m_plus_nd, 1, d_velocity_1, 1);
          cblas_dcopy(n_dminus2, rhs + m_plus_nd + n_dminus2, 1, d_velocity_2, 1);

          for (size_t i = 0; i < n; i++) {
            id5 = i * d;
            id3 = i * d_minus_2;

            d_velocity[id5] = d_velocity_1[id3] + d_velocity_2[id3];
            d_velocity[id5 + 1] = d_velocity_1[id3 + 1];
            d_velocity[id5 + 2] = d_velocity_1[id3 + 2];
            d_velocity[id5 + 3] = d_velocity_2[id3 + 1];
            d_velocity[id5 + 4] = d_velocity_2[id3 + 2];

            d_t[i] = d_velocity_1[id3];
            d_t_prime[i] = d_velocity_2[id3];
          }

          break;
        }  // end of SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2

        case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv: {
          /* -------------------------- Compute NT directions -------------------------- */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            Qp_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp2_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp2_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            NM_triplet_alloc(Qp_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp2_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp2_tilde, d_minus_2 * d_minus_2 * n);

            family_of_F(velocity_1, reaction_1, n_dminus2, n, NULL, NULL, Qp_bar, Qpinv_bar,
                        Qp2_bar, NULL);
            family_of_F(velocity_2, reaction_2, n_dminus2, n, NULL, NULL, Qp_tilde,
                        Qpinv_tilde, Qp2_tilde, NULL);
          } else if (options->iparam
                         [SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                     SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            JQinv =
                compute_JQinv(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n);
            JQinvT = NM_transpose(JQinv);
          }

          /* -------------------------- FIRST linear system -------------------------- */
          /*  1. Build the Jacobian matrix
           *
           *           m     nd      n(d+1)
           *        |  M     -H'       0    | m
           *        |                       |
           *  Jac = | -H      0      J*Q^-1 | nd
           *        |                       |
           *        |  0    Q^-1*J'    I    | n(d+1)
           *
           *          n(d-2)    n(d-2)
           *        | Qp_bar      0    | n(d-2)
           *   Q  = |                  |
           *        |  0      Qp_tilde | n(d-2)
           */

          Jac = NM_create(NM_SPARSE, m + nd + n_dplus1, m + nd + n_dplus1);
          Jac_nzmax = M_nzmax + 2 * H_nzmax + nd * n_dplus1 + n_dplus1;
          NM_triplet_alloc(Jac, Jac_nzmax);
          Jac->matrix2->origin = NSM_TRIPLET;
          NM_insert(Jac, M, 0, 0);
          NM_insert(Jac, minus_Ht, 0, m);
          NM_insert(Jac, minus_H, m, 0);
          NM_insert(Jac, JQinv, m, m_plus_nd);
          NM_insert(Jac, JQinvT, m_plus_nd, m);
          NM_insert(Jac, identity, m_plus_nd, m_plus_nd);

          /* Correction of w to take into account the dependence on the tangential velocity */
          // update_w(w, w_origin, velocity, nd, d,
          // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]);

          /*  2. Build the right-hand-side
           *
           *  rhs = -     <==== ATTENTION to negative sign
           *        [ M*v - H'*r + f ]  m         dualConstraint        = a (No have negative
           * sign at the beginning) [  u - H*v - w   ]  nd        primalConstraint      = b [
           * Qp_bar*r_1   ]  n(d-2)    complemConstraint 1 | [  Qp_tilde*r_2  ]  n(d-2)
           * complemConstraint 2 | = c
           */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv with formula F: not "
                "yet!\n\n");
          } else if (options->iparam
                         [SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                     SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinvz(velocity_1, reaction_1, reaction_1, n_dminus2, n,
                     complemConstraint_1);  // complemConstraint_1 = Qpinv_bar*r_1
            QNTpinvz(velocity_2, reaction_2, reaction_2, n_dminus2, n, complemConstraint_2);
          }

          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, primalConstraint, 1, rhs + m, 1);
          cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs + m_plus_nd, 1);
          cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs + m_plus_nd + n_dminus2, 1);

          cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
          cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_save, 1);

          /* 3. Solve full symmetric Newton system with NT scaling via LDLT factorization */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          NSM_linearSolverParams(Jac)->solver = NSM_HSL;
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS1_m = dnrm2l(m, rhs_save);
          residu_LS1_nd = dnrm2l(nd, rhs_save + m);
          residu_LS1_ndplus1 = dnrm2l(n_dplus1, rhs_save + m + nd);

          /* 4. Retrieve the directions for PREDICTOR step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);

          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          cblas_dcopy(n_dminus2, rhs + m_plus_nd, 1, tmp_n_dminus2_1,
                      1);  // tmp_n_dminus2_1 <-- d^_velocity_1 = Qp_bar * d_velocity_1
          cblas_dcopy(n_dminus2, rhs + m_plus_nd + n_dminus2, 1, tmp_n_dminus2_2, 1);

          // Recover d_velocity
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv with formula F: not "
                "yet!\n\n");

          } else if (options->iparam
                         [SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                     SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinvz(velocity_1, reaction_1, tmp_n_dminus2_1, n_dminus2, n,
                     d_velocity_1);  // d_velocity_1 = Qpinv_bar * d^_velocity_1
            QNTpinvz(velocity_2, reaction_2, tmp_n_dminus2_2, n_dminus2, n, d_velocity_2);
          }

          for (size_t i = 0; i < n; i++) {
            id5 = i * d;
            id3 = i * d_minus_2;
            d_velocity[id5] = d_velocity_1[id3] + d_velocity_2[id3];
            d_velocity[id5 + 1] = d_velocity_1[id3 + 1];
            d_velocity[id5 + 2] = d_velocity_1[id3 + 2];
            d_velocity[id5 + 3] = d_velocity_2[id3 + 1];
            d_velocity[id5 + 4] = d_velocity_2[id3 + 2];
          }

          /* 5. Computing the affine step-length */
          alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
          alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
          alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
          alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

          alpha_primal =
              fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
          alpha_dual = alpha_primal;

          /* updating the gamma parameter used to compute the step-length */
          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* -------------------------- Predictor step of Mehrotra --------------------------
           */
          cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
          cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
          cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          /* affine barrier parameter */
          barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / (n);

          /* computing the centralization parameter */
          // e = barr_param > sgmp1 ? fmax(1.0, 2. * alpha_primal) : sgmp3;
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * alpha_primal * alpha_primal) : sgmp3;
          // e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(alpha_primal, sgmp4)) : sgmp3;
          // sigma = fmin(1.0, pow(barr_param_a / barr_param, e))/sgmp5;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e)) / 5;

          /* -------------------------- SECOND linear system -------------------------- */
          // 6. Update the RHS
          // 1st term: y_check
          // already assigned in r1_check & r2_check

          /* 2nd terms:          (z_hat)^-1 o [dz_hat_a            o dy_check_a           - 2 *
           * barr_param * sigma * e] } =   (Qp_bar * u_1)^-1 o [(Qp_bar   * du_1)   o
           * (Qpinv_bar   * dr_1) - 2 * barr_param * sigma * e] } (Qp_tilde * u_2)^-1 o
           * [(Qp_tilde * du_2)   o (Qpinv_tilde * dr_2) - 2 * barr_param * sigma * e] }
           */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // Compute (z_hat)^-1
            NM_gemv(1.0, Qp_bar, velocity_1, 0.0,
                    velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            NM_gemv(1.0, Qp_tilde, velocity_2, 0.0,
                    velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            // Compute Jordan product dz_hat_a o dy_check_a
            NM_gemv(1.0, Qp_bar, d_velocity_1, 0.0,
                    d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            NM_gemv(1.0, Qp_tilde, d_velocity_2, 0.0,
                    d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2
            NM_gemv(1.0, Qpinv_bar, d_reaction_1, 0.0,
                    d_reaction_1_check);  // d_reaction_1_check     = Qpinv_bar * dr_1
            NM_gemv(1.0, Qpinv_tilde, d_reaction_2, 0.0,
                    d_reaction_2_check);  // d_reaction_2_check     = Qpinv_tilde * dr_2
            JA_prod(d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                    dvdr_jprod_1);  // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(
                d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                dvdr_jprod_2);  // dvdr_jprod_2      = (Qp_tilde * du_2) o (Qpinv_tilde * dr_2)

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o [(Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_2) - 2
                                                          // * mu * sigma * e]

            // Compute 2nd term
            NM_gemv(1.0, Qpinv_bar, velocity_1_hat_inv_dvhat_drcheck_1, 0.0,
                    tmp_n_dminus2_1);  // tmp_n_dminus2_1 = 2nd term_1
            NM_gemv(1.0, Qpinv_tilde, velocity_2_hat_inv_dvhat_drcheck_2, 0.0,
                    tmp_n_dminus2_2);  // tmp_n_dminus2_2 = 2nd term_2

            cblas_dcopy(n_dminus2, tmp_n_dminus2_1, 1, velocity_hat_inv_dvhat_drcheck, 1);
            cblas_dcopy(n_dminus2, tmp_n_dminus2_2, 1,
                        velocity_hat_inv_dvhat_drcheck + n_dminus2,
                        1);  // velocity_hat_inv_dvhat_drcheck = 2nd term

            /* Update 2nd row = P^-1 * [ (-H*v-w) - J*2nd term ] */
            NM_gemv(-1.0, J, velocity_hat_inv_dvhat_drcheck, 1.0,
                    Hvw);  // Hvw = (-H*v-w) - J*2nd term
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            // Compute (z_hat)^-1
            QNTpz(velocity_1, reaction_1, velocity_1, n_dminus2, n,
                  velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            QNTpz(velocity_2, reaction_2, velocity_2, n_dminus2, n, velocity_2_hat);
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n, velocity_2_hat_inv);

            // Compute Jordan product dz_hat_a o dy_check_a
            QNTpz(velocity_1, reaction_1, d_velocity_1, n_dminus2, n,
                  d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            QNTpz(velocity_2, reaction_2, d_velocity_2, n_dminus2, n, d_velocity_2_hat);
            QNTpinvz(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                     d_reaction_1_check);  // d_reaction_1_check   = Qpinv_bar * dr_1
            QNTpinvz(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, d_reaction_2_check);
            JA_prod(
                d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                dvdr_jprod_1);  // dvdr_jprod_1         = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(d_velocity_2_hat, d_reaction_2_check, n_dminus2, n, dvdr_jprod_2);

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute 2nd term = (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);

            // Update Complementarity constraints
            cblas_daxpy(n_dminus2, 1.0, velocity_1_hat_inv_dvhat_drcheck_1, 1,
                        complemConstraint_1, 1);
            cblas_daxpy(n_dminus2, 1.0, velocity_2_hat_inv_dvhat_drcheck_2, 1,
                        complemConstraint_2, 1);
          }

          // Update rhs
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, primalConstraint, 1, rhs + m, 1);
          cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs + m_plus_nd, 1);
          cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs + m_plus_nd + n_dminus2, 1);

          cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
          cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_save, 1);

          /* 7. Solve the 2nd linear system */
          // print_NAN_in_matrix(Jac);
          // if (NV_isnan(rhs, m + nd + n_dplus1)) printf("(2nd sys) NaN in RHS before solving,
          // i = %zu\n", iteration);

          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
              SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
            double *rhs_TMP = (double *)calloc(m + nd + n_dplus1, sizeof(double));
            cblas_dcopy(m + nd + n_dplus1, rhs, 1, rhs_TMP, 1);
            NM_LDLT_refine(Jac, rhs, rhs_TMP, 1, 1e-12, 10, 0);
            free(rhs_TMP);
          } else
            NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS2_m = dnrm2l(m, rhs_save);
          residu_LS2_nd = dnrm2l(nd, rhs_save + m);
          residu_LS2_ndplus1 = dnrm2l(n_dplus1, rhs_save + m + nd);

          /* 8. Retrieve the directions for CORRECTOR step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);

          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          cblas_dcopy(n_dminus2, rhs + m_plus_nd, 1, tmp_n_dminus2_1,
                      1);  // tmp_n_dminus2_1 <-- d^_velocity_1 = Qp_bar * d_velocity_1
          cblas_dcopy(n_dminus2, rhs + m_plus_nd + n_dminus2, 1, tmp_n_dminus2_2, 1);

          // Recover d_velocity
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv with formula F: not "
                "yet!\n\n");

          } else if (options->iparam
                         [SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                     SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinvz(velocity_1, reaction_1, tmp_n_dminus2_1, n_dminus2, n,
                     d_velocity_1);  // d_velocity_1 = Qpinv_bar * d^_velocity_1
            QNTpinvz(velocity_2, reaction_2, tmp_n_dminus2_2, n_dminus2, n, d_velocity_2);
          }

          for (size_t i = 0; i < n; i++) {
            id5 = i * d;
            id3 = i * d_minus_2;
            d_velocity[id5] = d_velocity_1[id3] + d_velocity_2[id3];
            d_velocity[id5 + 1] = d_velocity_1[id3 + 1];
            d_velocity[id5 + 2] = d_velocity_1[id3 + 2];
            d_velocity[id5 + 3] = d_velocity_2[id3 + 1];
            d_velocity[id5 + 4] = d_velocity_2[id3 + 2];

            d_t[i] = d_velocity_1[id3];
            d_t_prime[i] = d_velocity_2[id3];
          }

          break;
        }  // end of SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv

        case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ: {
          /* -------------------------- Compute NT directions -------------------------- */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            Qp_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv2_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv2_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            NM_triplet_alloc(Qp_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv2_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv2_tilde, d_minus_2 * d_minus_2 * n);

            double *f_NT = (double *)calloc(n_dminus2, sizeof(double));
            double *g_NT = (double *)calloc(n_dminus2, sizeof(double));
            float_type *wf_NT = (float_type *)calloc(n, sizeof(float_type));
            float_type *wg_NT = (float_type *)calloc(n, sizeof(float_type));

            family_of_F(velocity_1, reaction_1, n_dminus2, n, f_NT, wf_NT, Qp_bar, Qpinv_bar,
                        NULL, Qpinv2_bar);
            family_of_F(velocity_2, reaction_2, n_dminus2, n, g_NT, wg_NT, Qp_tilde,
                        Qpinv_tilde, NULL, Qpinv2_tilde);

            P_inv_F = Pinv_F(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n);

            free(f_NT);
            free(g_NT);
            free(wf_NT);
            free(wg_NT);
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            JQJ =
                compute_JQinv2Jt(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n);
          }

          /* -------------------------- FIRST linear system -------------------------- */
          /*  1. Build the reduced Jacobian matrix
           *
           *            m       nd
           *        |  -M       H'     | m
           *  Jac = |                  |
           *        |   H   J*Q^-2*J'  | nd
           */
          Jac = NM_create(NM_SPARSE, m + nd, m + nd);
          Jac_nzmax = M_nzmax + 2 * H_nzmax + nd * nd;
          NM_triplet_alloc(Jac, Jac_nzmax);
          Jac->matrix2->origin = NSM_TRIPLET;
          NM_insert(Jac, minus_M, 0, 0);
          NM_insert(Jac, Ht, 0, m);
          NM_insert(Jac, H, m, 0);
          NM_insert(Jac, JQJ, m, m);

          if (JQJ) JQJ = NM_free(JQJ);

          /* Correction of w to take into account the dependence on the tangential velocity */
          update_w(w, w_origin, velocity, nd, d,
                   options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]);

          /*  2. Build the right-hand-side
           *
           *  rhs = +     <==== positive sign
           *        [ M*v - H'*r + f ]  m
           *        [     -H*v-w     ]  nd
           */
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);

          NM_gemv(-1.0, H, globalVelocity, 0.0, Hvw);  // Hvw = -H*v
          cblas_daxpy(nd, -1.0, w, 1, Hvw, 1);         // Hvw = -H*v - w
          cblas_dcopy(nd, Hvw, 1, rhs + m, 1);
          cblas_dcopy(m + nd, rhs, 1, rhs_save, 1);

          /* 3. Solving full symmetric Newton system with NT scaling via LDLT factorization */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          NSM_linearSolverParams(Jac)->solver = NSM_HSL;
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS1_m = dnrm2l(m, rhs_save);
          residu_LS1_nd = dnrm2l(nd, rhs_save + m);

          /* 4. Retrieve the solutions for predictor step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);

          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];
          }

          /* 5. Computing the affine step-length */
          alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
          alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
          alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
          alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

          alpha_primal =
              fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
          alpha_dual = alpha_primal;

          /* updating the gamma parameter used to compute the step-length */
          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* -------------------------- Predictor step of Mehrotra --------------------------
           */
          cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
          cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
          cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          /* affine barrier parameter */
          barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / (n);

          /* computing the centralization parameter */
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * alpha_primal * alpha_primal) : sgmp3;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e)) / d;

          /* -------------------------- SECOND linear system -------------------------- */
          // 6. Update the RHS
          /* 1st terms: [-H*v-w] */
          // already in var Hvw = -Hv - w

          /* 2nd terms: Qpinv * {         (z_hat)^-1 o [dz_hat_a            o dy_check_a - 2 *
           * barr_param * sigma * e] } = Qpinv * {  (Qp_bar * u_1)^-1 o [(Qp_bar   * du_1)   o
           * (Qpinv_bar   * dr_1) - 2 * barr_param * sigma * e] } Qpinv * {(Qp_tilde * u_2)^-1
           * o [(Qp_tilde * du_2)   o (Qpinv_tilde * dr_2) - 2 * barr_param * sigma * e] }
           */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // Compute (z_hat)^-1
            NM_gemv(1.0, Qp_bar, velocity_1, 0.0,
                    velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            NM_gemv(1.0, Qp_tilde, velocity_2, 0.0,
                    velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            // Compute Jordan product dz_hat_a o dy_check_a
            NM_gemv(1.0, Qp_bar, d_velocity_1, 0.0,
                    d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            NM_gemv(1.0, Qp_tilde, d_velocity_2, 0.0,
                    d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2
            NM_gemv(1.0, Qpinv_bar, d_reaction_1, 0.0,
                    d_reaction_1_check);  // d_reaction_1_check     = Qpinv_bar * dr_1
            NM_gemv(1.0, Qpinv_tilde, d_reaction_2, 0.0,
                    d_reaction_2_check);  // d_reaction_2_check     = Qpinv_tilde * dr_2
            JA_prod(d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                    dvdr_jprod_1);  // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(
                d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                dvdr_jprod_2);  // dvdr_jprod_2      = (Qp_tilde * du_2) o (Qpinv_tilde * dr_2)

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o [(Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_2) - 2
                                                          // * mu * sigma * e]

            // Compute 2nd term
            NM_gemv(1.0, Qpinv_bar, velocity_1_hat_inv_dvhat_drcheck_1, 0.0,
                    tmp_n_dminus2_1);  // tmp_n_dminus2_1 = 2nd term_1
            NM_gemv(1.0, Qpinv_tilde, velocity_2_hat_inv_dvhat_drcheck_2, 0.0,
                    tmp_n_dminus2_2);  // tmp_n_dminus2_2 = 2nd term_2

            cblas_dcopy(n_dminus2, tmp_n_dminus2_1, 1, velocity_hat_inv_dvhat_drcheck, 1);
            cblas_dcopy(n_dminus2, tmp_n_dminus2_2, 1,
                        velocity_hat_inv_dvhat_drcheck + n_dminus2,
                        1);  // velocity_hat_inv_dvhat_drcheck = 2nd term

            /* Update 2nd row = P^-1 * [ (-H*v-w) - J*2nd term ] */
            NM_gemv(-1.0, J, velocity_hat_inv_dvhat_drcheck, 1.0,
                    Hvw);  // Hvw = (-H*v-w) - J*2nd term
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            // Compute (z_hat)^-1
            QNTpz(velocity_1, reaction_1, velocity_1, n_dminus2, n,
                  velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            QNTpz(velocity_2, reaction_2, velocity_2, n_dminus2, n,
                  velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            // Compute Jordan product dz_hat_a o dy_check_a
            QNTpz(velocity_1, reaction_1, d_velocity_1, n_dminus2, n,
                  d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            QNTpz(velocity_2, reaction_2, d_velocity_2, n_dminus2, n,
                  d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2
            QNTpinvz(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                     d_reaction_1_check);  // d_reaction_1_check   = Qpinv_bar * dr_1
            QNTpinvz(velocity_2, reaction_2, d_reaction_2, n_dminus2, n,
                     d_reaction_2_check);  // d_reaction_2_check   = Qpinv_tilde * dr_2
            JA_prod(
                d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                dvdr_jprod_1);  // dvdr_jprod_1         = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                    dvdr_jprod_2);  // dvdr_jprod_2         = (Qp_tilde * du_2) o (Qpinv_tilde
                                    // * dr_2)

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o [(Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_2) - 2
                                                          // * mu * sigma * e]

            // Compute 2nd term
            QNTpinvz(velocity_1, reaction_1, velocity_1_hat_inv_dvhat_drcheck_1, n_dminus2, n,
                     tmp_n_dminus2_1);  // tmp_n_dminus2_1 = 2nd term_1
            QNTpinvz(velocity_2, reaction_2, velocity_2_hat_inv_dvhat_drcheck_2, n_dminus2, n,
                     tmp_n_dminus2_2);  // tmp_n_dminus2_2 = 2nd term_2

            cblas_dcopy(n_dminus2, tmp_n_dminus2_1, 1, velocity_hat_inv_dvhat_drcheck, 1);
            cblas_dcopy(n_dminus2, tmp_n_dminus2_2, 1,
                        velocity_hat_inv_dvhat_drcheck + n_dminus2,
                        1);  // velocity_hat_inv_dvhat_drcheck = 2nd term

            /* Update 2nd row = (-H*v-w) - J*2nd term */
            NM_gemv(-1.0, J, velocity_hat_inv_dvhat_drcheck, 1.0,
                    Hvw);  // Hvw = (-H*v-w) - J*2nd term
          }

          // Update rhs
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, Hvw, 1, rhs + m, 1);

          cblas_dcopy(m + nd, rhs, 1, rhs_save, 1);

          /* 7. Solve the 2nd linear system */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
              SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
            cblas_dcopy(m + nd, rhs, 1, rhs_save, 1);
            NM_LDLT_refine(J, rhs, rhs_save, 1, tol, 10, 0);
          }

          else
            NM_LDLT_solve(Jac, rhs, 1);

          // if (NV_isnan(rhs, m + nd)) printf("(2nd sys) NaN in RHS after solving, i = %zu\n",
          // iteration);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS2_m = dnrm2l(m, rhs_save);
          residu_LS2_nd = dnrm2l(nd, rhs_save + m);

          /* 8. Retrieve the solutions for predictor step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);

          cblas_dcopy(nd, rhs + m, 1, d_reaction, 1);
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1 - 2nd term
            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            // id5 = i*d;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];
            // d_velocity[id5]   = d_velocity_1[id3] + d_velocity_2[id3];

            d_t[i] = d_velocity_1[id3];
            d_t_prime[i] = d_velocity_2[id3];
          }

          break;
        }  // end of SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ

        case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH: {
          /* -------------------------- Compute NT directions -------------------------- */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            Qp_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv2_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qp_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            Qpinv2_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
            NM_triplet_alloc(Qp_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv2_bar, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qp_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv_tilde, d_minus_2 * d_minus_2 * n);
            NM_triplet_alloc(Qpinv2_tilde, d_minus_2 * d_minus_2 * n);

            double *f_NT = (double *)calloc(n_dminus2, sizeof(double));
            double *g_NT = (double *)calloc(n_dminus2, sizeof(double));
            float_type *wf_NT = (float_type *)calloc(n, sizeof(float_type));
            float_type *wg_NT = (float_type *)calloc(n, sizeof(float_type));

            family_of_F(velocity_1, reaction_1, n_dminus2, n, f_NT, wf_NT, Qp_bar, Qpinv_bar,
                        NULL, Qpinv2_bar);
            family_of_F(velocity_2, reaction_2, n_dminus2, n, g_NT, wg_NT, Qp_tilde,
                        Qpinv_tilde, NULL, Qpinv2_tilde);

            P_inv_F = Pinv_F(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n);

            free(f_NT);
            free(g_NT);
            free(wf_NT);
            free(wg_NT);
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES) {
              chol_U = compute_factor_U(velocity_1, reaction_1, velocity_2, reaction_2,
                                        n_dminus2, n);
              chol_U_csc = NM_csc(chol_U);

              PinvH = multiply_UinvH(chol_U_csc, H);
              PinvH_T = NM_transpose(PinvH);
            }

            else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                     SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO) {
              P_inv = Pinv(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n);
              PinvH = NM_multiply(P_inv, H);

              // PinvH = multiply_PinvH(velocity_1, reaction_1, velocity_2, reaction_2,
              // n_dminus2, n, H);     // Same results with NM_multiply but executive time is
              // so dramatic !
              PinvH_T = NM_transpose(PinvH);
            }
          }

          /* -------------------------- FIRST linear system -------------------------- */
          /*  1. Build the reduced Jacobian matrix
           *
           *            m        nd
           *        |  -M     H'P^-1' | m
           *  Jac = |                 |
           *        | P^-1H      I    | nd
           */
          Jac = NM_create(NM_SPARSE, m + nd, m + nd);
          Jac_nzmax = M_nzmax + 2 * m * nd + nd;
          NM_triplet_alloc(Jac, Jac_nzmax);
          Jac->matrix2->origin = NSM_TRIPLET;
          NM_insert(Jac, minus_M, 0, 0);
          NM_insert(Jac, PinvH_T, 0, m);
          NM_insert(Jac, PinvH, m, 0);
          NM_insert(Jac, identity, m, m);

          /* Correction of w to take into account the dependence on the tangential velocity */
          update_w(w, w_origin, velocity, nd, d,
                   options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]);

          /*  2. Build the right-hand-side
           *
           *  rhs = +     <==== positive sign
           *        [ M*v - H'*r + f ]  m
           *        [  P^-1*(-H*v-w) ]  nd
           */
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          NM_gemv(-1.0, H, globalVelocity, 0.0, Hvw);  // Hvw = -H*v
          cblas_daxpy(nd, -1.0, w, 1, Hvw, 1);         // Hvw = -H*v - w

          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH with formula F: not "
                "yet!\n\n");

          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES) {
              if (chol_L || chol_U) {
                cblas_dcopy(nd, Hvw, 1, rhs + m, 1);
                cs_usolve(chol_U_csc, rhs + m);
              } else {
                printf("\n Cholesky factor is NULL! 2. Build the right-hand-side \n");
                break;
              }

            } else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                       SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO) {
              Pinvy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, Hvw,
                    rhs + m);  // rhs+m = P^-1*(-H*v-w)
            }
          }

          cblas_dcopy(m + nd, rhs, 1, rhs_save, 1);

          /* 3. Solving full symmetric Newton system with NT scaling via LDLT factorization */

          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          NSM_linearSolverParams(Jac)->solver = NSM_HSL;
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS1_m = dnrm2l(m, rhs_save);
          residu_LS1_nd = dnrm2l(nd, rhs_save + m);

          /* 4. Retrieve the solutions for predictor step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
          cblas_dcopy(nd, rhs + m, 1, tmp_nd,
                      1);  // tmp_nd <-- d_reaction_reduced = P'*d_reaction

          // Recover d_reaction
          // NM_tgemv(1.0, P_inv, tmp_nd, 0.0, d_reaction);      // d_reaction =
          // P^-1'*d_reaction_reduced
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH with formula F: not "
                "yet!\n\n");

          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES) {
              if (chol_L || chol_U) {
                cblas_dcopy(nd, tmp_nd, 1, d_reaction, 1);
                cs_utsolve(chol_U_csc, d_reaction);
              } else {
                printf("\n Cholesky factor is NULL! \n");
                break;
              }

            } else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                       SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO)
              PinvTy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, tmp_nd,
                     d_reaction);  // d_reaction = P^-1'*d_reaction_reduced
          }
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH with formula F: not "
                "yet!\n\n");

          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];
          }

          /* 5. Computing the affine step-length */
          alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
          alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
          alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
          alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

          alpha_primal =
              fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
          alpha_dual = alpha_primal;

          /* updating the gamma parameter used to compute the step-length */
          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* -------------------------- Predictor step of Mehrotra --------------------------
           */
          cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
          cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
          cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          /* affine barrier parameter */
          barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / (n);

          /* computing the centralization parameter */
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * alpha_primal * alpha_primal) : sgmp3;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e)) / d;

          /* -------------------------- SECOND linear system -------------------------- */
          // 6. Update the RHS
          /* 1st terms: [-H*v-w] */
          // already in var Hvw = -Hv - w

          /* 2nd terms: Qpinv * {         (z_hat)^-1 o [dz_hat_a            o dy_check_a - 2 *
           * barr_param * sigma * e] } = Qpinv * {  (Qp_bar * u_1)^-1 o [(Qp_bar   * du_1)   o
           * (Qpinv_bar   * dr_1) - 2 * barr_param * sigma * e] } Qpinv * {(Qp_tilde * u_2)^-1
           * o [(Qp_tilde * du_2)   o (Qpinv_tilde * dr_2) - 2 * barr_param * sigma * e] }
           */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // Compute (z_hat)^-1
            NM_gemv(1.0, Qp_bar, velocity_1, 0.0,
                    velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            NM_gemv(1.0, Qp_tilde, velocity_2, 0.0,
                    velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            // Compute Jordan product dz_hat_a o dy_check_a
            NM_gemv(1.0, Qp_bar, d_velocity_1, 0.0,
                    d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            NM_gemv(1.0, Qp_tilde, d_velocity_2, 0.0,
                    d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2
            NM_gemv(1.0, Qpinv_bar, d_reaction_1, 0.0,
                    d_reaction_1_check);  // d_reaction_1_check     = Qpinv_bar * dr_1
            NM_gemv(1.0, Qpinv_tilde, d_reaction_2, 0.0,
                    d_reaction_2_check);  // d_reaction_2_check     = Qpinv_tilde * dr_2
            JA_prod(d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                    dvdr_jprod_1);  // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(
                d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                dvdr_jprod_2);  // dvdr_jprod_2      = (Qp_tilde * du_2) o (Qpinv_tilde * dr_2)

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o [(Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_2) - 2
                                                          // * mu * sigma * e]

            // Compute 2nd term
            NM_gemv(1.0, Qpinv_bar, velocity_1_hat_inv_dvhat_drcheck_1, 0.0,
                    tmp_n_dminus2_1);  // tmp_n_dminus2_1 = 2nd term_1
            NM_gemv(1.0, Qpinv_tilde, velocity_2_hat_inv_dvhat_drcheck_2, 0.0,
                    tmp_n_dminus2_2);  // tmp_n_dminus2_2 = 2nd term_2

            cblas_dcopy(n_dminus2, tmp_n_dminus2_1, 1, velocity_hat_inv_dvhat_drcheck, 1);
            cblas_dcopy(n_dminus2, tmp_n_dminus2_2, 1,
                        velocity_hat_inv_dvhat_drcheck + n_dminus2,
                        1);  // velocity_hat_inv_dvhat_drcheck = 2nd term

            /* Update 2nd row = P^-1 * [ (-H*v-w) - J*2nd term ] */
            NM_gemv(-1.0, J, velocity_hat_inv_dvhat_drcheck, 1.0,
                    Hvw);  // Hvw = (-H*v-w) - J*2nd term
            NM_gemv(1.0, P_inv, Hvw, 0.0,
                    tmp_nd);  // tmp_nd = P^-1 * [ (-H*v-w) - J*2nd term ]
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            // Compute (z_hat)^-1
            QNTpz(velocity_1, reaction_1, velocity_1, n_dminus2, n,
                  velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            QNTpz(velocity_2, reaction_2, velocity_2, n_dminus2, n,
                  velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            // Compute Jordan product dz_hat_a o dy_check_a
            QNTpz(velocity_1, reaction_1, d_velocity_1, n_dminus2, n,
                  d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            QNTpz(velocity_2, reaction_2, d_velocity_2, n_dminus2, n,
                  d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2
            QNTpinvz(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                     d_reaction_1_check);  // d_reaction_1_check   = Qpinv_bar * dr_1
            QNTpinvz(velocity_2, reaction_2, d_reaction_2, n_dminus2, n,
                     d_reaction_2_check);  // d_reaction_2_check   = Qpinv_tilde * dr_2
            JA_prod(
                d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                dvdr_jprod_1);  // dvdr_jprod_1         = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                    dvdr_jprod_2);  // dvdr_jprod_2         = (Qp_tilde * du_2) o (Qpinv_tilde
                                    // * dr_2)

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o [(Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_2) - 2
                                                          // * mu * sigma * e]

            // Compute 2nd term
            QNTpinvz(velocity_1, reaction_1, velocity_1_hat_inv_dvhat_drcheck_1, n_dminus2, n,
                     tmp_n_dminus2_1);  // tmp_n_dminus2_1 = 2nd term_1
            QNTpinvz(velocity_2, reaction_2, velocity_2_hat_inv_dvhat_drcheck_2, n_dminus2, n,
                     tmp_n_dminus2_2);  // tmp_n_dminus2_2 = 2nd term_2

            cblas_dcopy(n_dminus2, tmp_n_dminus2_1, 1, velocity_hat_inv_dvhat_drcheck, 1);
            cblas_dcopy(n_dminus2, tmp_n_dminus2_2, 1,
                        velocity_hat_inv_dvhat_drcheck + n_dminus2,
                        1);  // velocity_hat_inv_dvhat_drcheck = 2nd term

            /* Update 2nd row = P^-1 * [ (-H*v-w) - J*2nd term ] */
            NM_gemv(-1.0, J, velocity_hat_inv_dvhat_drcheck, 1.0,
                    Hvw);  // Hvw = (-H*v-w) - J*2nd term
            // NM_gemv(1.0, P_inv, Hvw, 0.0, tmp_nd);                     // tmp_nd = P^-1 * [
            // (-H*v-w) - J*2nd term ]

            if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES) {
              if (chol_L || chol_U) {
                cblas_dcopy(nd, Hvw, 1, tmp_nd, 1);
                cs_usolve(chol_U_csc, tmp_nd);

                // cs_lsolve (chol_L, tmp_nd) ;
              } else {
                printf("\n Cholesky factor is NULL! \n");
                break;
              }
            } else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                       SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO)
              Pinvy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, Hvw,
                    tmp_nd);  // tmp_nd = P^-1 * [ (-H*v-w) - J*2nd term ]
          }

          // Update rhs
          cblas_dcopy(m, dualConstraint, 1, rhs, 1);
          cblas_dcopy(nd, tmp_nd, 1, rhs + m, 1);

          cblas_dcopy(m + nd, rhs, 1, rhs_save, 1);

          /* 7. Solve the 2nd linear system */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
              SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
            cblas_dcopy(m + nd, rhs, 1, rhs_save, 1);
            NM_LDLT_refine(J, rhs, rhs_save, 1, tol, 10, 0);
          }

          else
            NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS2_m = dnrm2l(m, rhs_save);
          residu_LS2_nd = dnrm2l(nd, rhs_save + m);

          /* 8. Retrieve the solutions for predictor step */
          cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
          cblas_dcopy(nd, rhs + m, 1, tmp_nd,
                      1);  // tmp_nd <-- d_reaction_reduced = P'*d_reaction

          // Recover d_reaction
          // NM_tgemv(1.0, P_inv, tmp_nd, 0.0, d_reaction); // d_reaction =
          // P^-1'*d_reaction_reduced
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH with formula F: not "
                "yet!\n\n");

          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES) {
              if (chol_L || chol_U) {
                cblas_dcopy(nd, tmp_nd, 1, d_reaction, 1);
                cs_utsolve(chol_U_csc, d_reaction);

                // cs_ltsolve (chol_L, d_reaction) ;
              } else {
                printf("\n Cholesky factor is NULL! \n");
                break;
              }

            } else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] ==
                       SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO)
              PinvTy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, tmp_nd,
                     d_reaction);  // d_reaction = P^-1'*d_reaction_reduced
          }

          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH with formula F: not "
                "yet!\n\n");

          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1 - 2nd term
            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];

            d_t[i] = d_velocity_1[id3];
            d_t_prime[i] = d_velocity_2[id3];
          }

          break;
        }  // end of SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH

        case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ: {
          /* -------------------------- Compute NT directions -------------------------- */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            printf(
                "\n\nSICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ with formula F: NOT YET\n\n");
            break;
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            JQJ =
                compute_JQinv2Jt(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n);
          }

          /* -------------------------- FIRST linear system -------------------------- */
          /*  1. Build the reduced Jacobian matrix
           *
           *  Jac = H*M^-1*H' + J*Qp^-2*J'
           *
           */
          Jac = NM_create(NM_SPARSE, nd, nd);
          Jac = NM_add(1.0, HMinvHt, 1.0, JQJ);

          if (JQJ) JQJ = NM_free(JQJ);

          /* Correction of w to take into account the dependence on the tangential velocity */
          update_w(w, w_origin, velocity, nd, d,
                   options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]);

          /*  2. Build the right-hand-side
           *
           *  rhs = -H*M^-1*(H'r+f) - w
           *
           */
          cblas_dcopy(m, f, 1, Hrf, 1);
          NM_gemv(1.0, Ht, reaction, 1.0, Hrf);    // Hrf = H'r + f
          NM_gemv(-1.0, HMinv, Hrf, 0.0, HMHrfw);  // HMHrfw = -H*M^-1*(H'r+f)
          cblas_daxpy(nd, -1.0, w, 1, HMHrfw, 1);  // HMHrfw = -H*M^-1*(H'r+f) - w

          cblas_dcopy(nd, HMHrfw, 1, rhs, 1);
          cblas_dcopy(nd, HMHrfw, 1, rhs_save, 1);

          /* 3. Solving full symmetric Newton system with NT scaling via LDLT factorization */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          // NM_Cholesky_solve(Jac, rhs, 1);
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS1_nd = dnrm2l(nd, rhs_save);

          /* 4. Retrieve the solutions for predictor step */
          cblas_dcopy(nd, rhs, 1, d_reaction, 1);  // d_reaction = rhs
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          NV_add(reaction, d_reaction, nd, rdr);  // rdr = r + dr
          cblas_dcopy(m, f, 1, MfHrdr, 1);        // MfHrdr = f
          NM_gemv(1.0, Ht, rdr, 1.0, MfHrdr);     // MfHrdr = f + H'(r+dr)
          NM_gemv(1.0, Minv, MfHrdr, 0.0,
                  d_globalVelocity);  // d_globalVelocity = M^-1*(f + H'(r+dr))
          cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity,
                      1);  // d_globalVelocity = -v + M^-1*(f + H'(r+dr))

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];
          }

          /* 5. Computing the affine step-length */
          alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
          alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
          alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
          alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

          alpha_primal =
              fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
          alpha_dual = alpha_primal;

          /* updating the gamma parameter used to compute the step-length */
          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* -------------------------- Predictor step of Mehrotra --------------------------
           */
          cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
          cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
          cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          /* affine barrier parameter */
          barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / (n);

          /* computing the centralization parameter */
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * alpha_primal * alpha_primal) : sgmp3;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e)) / d;

          /* -------------------------- SECOND linear system -------------------------- */
          /* 6. Update the RHS
           *
           * rhs = -H*M^-1*(H'r+f) - w - J*2nd terms
           *
           */

          /* 2nd terms: Qpinv * {         (z_hat)^-1 o [dz_hat_a            o dy_check_a - 2 *
           * barr_param * sigma * e] } = Qpinv * {  (Qp_bar * u_1)^-1 o [(Qp_bar   * du_1)   o
           * (Qpinv_bar   * dr_1) - 2 * barr_param * sigma * e] } Qpinv * {(Qp_tilde * u_2)^-1
           * o [(Qp_tilde * du_2)   o (Qpinv_tilde * dr_2) - 2 * barr_param * sigma * e] }
           */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            printf(
                "\n\nSICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ with formula F: NOT YET\n\n");
            break;
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            // Compute (z_hat)^-1
            QNTpz(velocity_1, reaction_1, velocity_1, n_dminus2, n,
                  velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            QNTpz(velocity_2, reaction_2, velocity_2, n_dminus2, n,
                  velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            // Compute Jordan product dz_hat_a o dy_check_a
            QNTpz(velocity_1, reaction_1, d_velocity_1, n_dminus2, n,
                  d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            QNTpz(velocity_2, reaction_2, d_velocity_2, n_dminus2, n,
                  d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2
            QNTpinvz(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                     d_reaction_1_check);  // d_reaction_1_check   = Qpinv_bar * dr_1
            QNTpinvz(velocity_2, reaction_2, d_reaction_2, n_dminus2, n,
                     d_reaction_2_check);  // d_reaction_2_check   = Qpinv_tilde * dr_2
            JA_prod(
                d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                dvdr_jprod_1);  // dvdr_jprod_1         = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                    dvdr_jprod_2);  // dvdr_jprod_2         = (Qp_tilde * du_2) o (Qpinv_tilde
                                    // * dr_2)

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o [(Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_2) - 2
                                                          // * mu * sigma * e]

            // Compute 2nd term
            QNTpinvz(velocity_1, reaction_1, velocity_1_hat_inv_dvhat_drcheck_1, n_dminus2, n,
                     tmp_n_dminus2_1);  // tmp_n_dminus2_1 = 2nd term_1
            QNTpinvz(velocity_2, reaction_2, velocity_2_hat_inv_dvhat_drcheck_2, n_dminus2, n,
                     tmp_n_dminus2_2);  // tmp_n_dminus2_2 = 2nd term_2

            cblas_dcopy(n_dminus2, tmp_n_dminus2_1, 1, velocity_hat_inv_dvhat_drcheck, 1);
            cblas_dcopy(n_dminus2, tmp_n_dminus2_2, 1,
                        velocity_hat_inv_dvhat_drcheck + n_dminus2,
                        1);  // velocity_hat_inv_dvhat_drcheck = 2nd term
          }

          // Update rhs = -H*M^-1*(H'r+f) - w - J*2nd term
          cblas_dcopy(nd, HMHrfw, 1, rhs, 1);
          NM_gemv(-1.0, J, velocity_hat_inv_dvhat_drcheck, 1.0, rhs);

          cblas_dcopy(nd, rhs, 1, rhs_save, 1);

          /* 7. Solve the 2nd linear system */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          // NM_Cholesky_solve(Jac, rhs, 1);
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS2_nd = dnrm2l(nd, rhs_save);

          /* 8. Retrieve the solutions for predictor step */
          cblas_dcopy(nd, rhs, 1, d_reaction, 1);  // d_reaction = rhs
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          NV_add(reaction, d_reaction, nd, rdr);  // rdr = r + dr
          cblas_dcopy(m, f, 1, MfHrdr, 1);        // MfHrdr = f
          NM_gemv(1.0, Ht, rdr, 1.0, MfHrdr);     // MfHrdr = f + H'(r+dr)
          NM_gemv(1.0, Minv, MfHrdr, 0.0,
                  d_globalVelocity);  // d_globalVelocity = M^-1*(f + H'(r+dr))
          cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity,
                      1);  // d_globalVelocity = -v + M^-1*(f + H'(r+dr))

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1 - 2nd term
            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];

            d_t[i] = d_velocity_1[id3];
            d_t_prime[i] = d_velocity_2[id3];
          }

          break;
        }  // end of SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ

        case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH: {
          /* -------------------------- Compute NT directions -------------------------- */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            printf(
                "\n\nSICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH with formula F: NOT YET\n\n");
            break;
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            P_inv = Pinv(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n);
          }

          /* -------------------------- FIRST linear system -------------------------- */
          /*  1. Build the reduced Jacobian matrix
           *
           *  Jac = P^-1H*M^-1*H'P^-1' + I
           *
           */
          Jac = NM_create(NM_SPARSE, nd, nd);

          P_invT = NM_transpose(P_inv);
          HMHP = NM_multiply(HMinvHt, P_invT);
          PHMHP = NM_multiply(P_inv, HMHP);

          Jac = NM_add(1.0, PHMHP, 1.0, identity);

          /* Correction of w to take into account the dependence on the tangential velocity */
          update_w(w, w_origin, velocity, nd, d,
                   options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]);

          /*  2. Build the right-hand-side
           *
           *  rhs = P^-1 * [ -H*M^-1*(H'r+f) - w ]
           *
           */
          cblas_dcopy(m, f, 1, Hrf, 1);
          NM_gemv(1.0, Ht, reaction, 1.0, Hrf);    // Hrf = H'r + f
          NM_gemv(-1.0, HMinv, Hrf, 0.0, HMHrfw);  // HMHrfw = -H*M^-1*(H'r+f)
          cblas_daxpy(nd, -1.0, w, 1, HMHrfw, 1);  // HMHrfw = -H*M^-1*(H'r+f) - w

          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            Pinvy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, HMHrfw,
                  rhs);  // rhs = P^-1 * [ -H*M^-1*(H'r+f) - w ]
          }

          cblas_dcopy(nd, rhs, 1, rhs_save, 1);

          /* 3. Solving full symmetric Newton system with NT scaling via LDLT factorization */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          // NM_Cholesky_solve(Jac, rhs, 1);
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS1_nd = dnrm2l(nd, rhs_save);

          /* 4. Retrieve the solutions for predictor step */
          cblas_dcopy(nd, rhs, 1, tmp_nd, 1);  // tmp_nd <-- d_reaction_reduced = P'*d_reaction
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            PinvTy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, tmp_nd,
                   d_reaction);  // d_reaction = P^-1'*d_reaction_reduced
          }
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          NV_add(reaction, d_reaction, nd, rdr);  // rdr = r + dr
          cblas_dcopy(m, f, 1, MfHrdr, 1);        // MfHrdr = f
          NM_gemv(1.0, Ht, rdr, 1.0, MfHrdr);     // MfHrdr = f + H'(r+dr)
          NM_gemv(1.0, Minv, MfHrdr, 0.0,
                  d_globalVelocity);  // d_globalVelocity = M^-1*(f + H'(r+dr))
          cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity,
                      1);  // d_globalVelocity = -v + M^-1*(f + H'(r+dr))

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];
          }

          /* 5. Computing the affine step-length */
          alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
          alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
          alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
          alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

          alpha_primal =
              fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
          alpha_dual = alpha_primal;

          /* updating the gamma parameter used to compute the step-length */
          gmm = gmmp1 + gmmp2 * alpha_primal;

          /* -------------------------- Predictor step of Mehrotra --------------------------
           */
          cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
          cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
          cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
          cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

          /* affine barrier parameter */
          barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / (n);

          /* computing the centralization parameter */
          e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * alpha_primal * alpha_primal) : sgmp3;
          sigma = fmin(1.0, pow(barr_param_a / barr_param, e)) / d;

          /* -------------------------- SECOND linear system -------------------------- */
          /* 6. Update the RHS
           *
           * rhs = -H*M^-1*(H'r+f) - w - J*2nd terms
           *
           */

          /* 2nd terms: Qpinv * {         (z_hat)^-1 o [dz_hat_a            o dy_check_a - 2 *
           * barr_param * sigma * e] } = Qpinv * {  (Qp_bar * u_1)^-1 o [(Qp_bar   * du_1)   o
           * (Qpinv_bar   * dr_1) - 2 * barr_param * sigma * e] } Qpinv * {(Qp_tilde * u_2)^-1
           * o [(Qp_tilde * du_2)   o (Qpinv_tilde * dr_2) - 2 * barr_param * sigma * e] }
           */
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            printf(
                "\n\nSICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH with formula F: NOT YET\n\n");
            break;
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            // Compute (z_hat)^-1
            QNTpz(velocity_1, reaction_1, velocity_1, n_dminus2, n,
                  velocity_1_hat);  // velocity_1_hat = Qp_bar * u_1
            QNTpz(velocity_2, reaction_2, velocity_2, n_dminus2, n,
                  velocity_2_hat);  // velocity_2_hat = Qp_tilde * u_2
            Jinv(velocity_1_hat, n_dminus2, n,
                 velocity_1_hat_inv);  // velocity_1_hat_inv = (Qp_bar * u_1)^-1
            Jinv(velocity_2_hat, n_dminus2, n,
                 velocity_2_hat_inv);  // velocity_2_hat_inv = (Qp_tilde * u_2)^-1

            // Compute Jordan product dz_hat_a o dy_check_a
            QNTpz(velocity_1, reaction_1, d_velocity_1, n_dminus2, n,
                  d_velocity_1_hat);  // d_velocity_1_hat     = Qp_bar * du_1
            QNTpz(velocity_2, reaction_2, d_velocity_2, n_dminus2, n,
                  d_velocity_2_hat);  // d_velocity_2_hat     = Qp_tilde * du_2
            QNTpinvz(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                     d_reaction_1_check);  // d_reaction_1_check   = Qpinv_bar * dr_1
            QNTpinvz(velocity_2, reaction_2, d_reaction_2, n_dminus2, n,
                     d_reaction_2_check);  // d_reaction_2_check   = Qpinv_tilde * dr_2
            JA_prod(
                d_velocity_1_hat, d_reaction_1_check, n_dminus2, n,
                dvdr_jprod_1);  // dvdr_jprod_1         = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
            JA_prod(d_velocity_2_hat, d_reaction_2_check, n_dminus2, n,
                    dvdr_jprod_2);  // dvdr_jprod_2         = (Qp_tilde * du_2) o (Qpinv_tilde
                                    // * dr_2)

            // Jordan product - 2 * mu * sigma * e
            sigma_mu = 2. * barr_param * sigma;
            for (size_t k = 0; k < n_dminus2; k += d_minus_2) {
              dvdr_jprod_1[k] -= sigma_mu;
              dvdr_jprod_2[k] -= sigma_mu;
            }

            // Compute (z_hat)^-1 o [...]
            JA_prod(velocity_1_hat_inv, dvdr_jprod_1, n_dminus2, n,
                    velocity_1_hat_inv_dvhat_drcheck_1);  // velocity_1_hat_inv_dvhat_drcheck_1
                                                          // =   (Qp_bar * u_1)^-1 o [(Qp_bar *
                                                          // du_1)   o (Qpinv_bar * dr_1)   - 2
                                                          // * mu * sigma * e]
            JA_prod(velocity_2_hat_inv, dvdr_jprod_2, n_dminus2, n,
                    velocity_2_hat_inv_dvhat_drcheck_2);  // velocity_2_hat_inv_dvhat_drcheck_2
                                                          // = (Qp_tilde * u_2)^-1 o [(Qp_tilde
                                                          // * du_2) o (Qpinv_tilde * dr_2) - 2
                                                          // * mu * sigma * e]

            // Compute 2nd term
            QNTpinvz(velocity_1, reaction_1, velocity_1_hat_inv_dvhat_drcheck_1, n_dminus2, n,
                     tmp_n_dminus2_1);  // tmp_n_dminus2_1 = 2nd term_1
            QNTpinvz(velocity_2, reaction_2, velocity_2_hat_inv_dvhat_drcheck_2, n_dminus2, n,
                     tmp_n_dminus2_2);  // tmp_n_dminus2_2 = 2nd term_2

            cblas_dcopy(n_dminus2, tmp_n_dminus2_1, 1, velocity_hat_inv_dvhat_drcheck, 1);
            cblas_dcopy(n_dminus2, tmp_n_dminus2_2, 1,
                        velocity_hat_inv_dvhat_drcheck + n_dminus2,
                        1);  // velocity_hat_inv_dvhat_drcheck = 2nd term
          }

          // Update rhs = P^-1 * [ -H*M^-1*(H'r+f) - w - J*2nd term ]
          cblas_dcopy(nd, HMHrfw, 1, tmp_nd, 1);
          NM_gemv(-1.0, J, velocity_hat_inv_dvhat_drcheck, 1.0, tmp_nd);
          Pinvy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, tmp_nd, rhs);

          cblas_dcopy(nd, rhs, 1, rhs_save, 1);

          /* 7. Solve the 2nd linear system */
          jacobian_is_nan = NM_isnan(Jac);
          if (jacobian_is_nan) {
            numerics_printf_verbose(0, "The Jacobian matrix contains NaN");
            break;
          }

          // NM_Cholesky_solve(Jac, rhs, 1);
          NM_LDLT_solve(Jac, rhs, 1);

          NM_gemv(1.0, Jac, rhs, -1.0, rhs_save);
          residu_LS2_nd = dnrm2l(nd, rhs_save);

          /* 8. Retrieve the solutions for predictor step */
          cblas_dcopy(nd, rhs, 1, tmp_nd, 1);  // tmp_nd <-- d_reaction_reduced = P'*d_reaction
          PinvTy(velocity_1, reaction_1, velocity_2, reaction_2, n_dminus2, n, tmp_nd,
                 d_reaction);  // d_reaction = P^-1'*d_reaction_reduced
          extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
          extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

          NV_add(reaction, d_reaction, nd, rdr);  // rdr = r + dr
          cblas_dcopy(m, f, 1, MfHrdr, 1);        // MfHrdr = f
          NM_gemv(1.0, Ht, rdr, 1.0, MfHrdr);     // MfHrdr = f + H'(r+dr)
          NM_gemv(1.0, Minv, MfHrdr, 0.0,
                  d_globalVelocity);  // d_globalVelocity = M^-1*(f + H'(r+dr))
          cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity,
                      1);  // d_globalVelocity = -v + M^-1*(f + H'(r+dr))

          // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
          // Recover du_bar & du_tilde
          NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity = H*dv
          cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity,
                      1);  // d_velocity = H*dv - (u-Hv-w)
          extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
          extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

          // Recover d_t & d_t'
          if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
              SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F) {
            // TO DO
            printf(
                "\n\n SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ with formula F: not yet!\n\n");
          }

          else if (options
                       ->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] ==
                   SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP) {
            QNTpinv2z(velocity_1, reaction_1, d_reaction_1, n_dminus2, n,
                      Qinv2x_bar);  // Qinv2x_bar = Qpbar^-2*d_reaction_1
            QNTpinv2z(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, Qinv2x_tilde);

            cblas_dscal(n_dminus2, -1.0, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1
            cblas_dscal(n_dminus2, -1.0, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, velocity_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1
            cblas_daxpy(n_dminus2, -1.0, velocity_2, 1, Qinv2x_tilde, 1);

            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_1, 1, Qinv2x_bar,
                        1);  // Qinv2x_bar = - Qpbar^-2*d_reaction_1 - velocity_1 - 2nd term
            cblas_daxpy(n_dminus2, -1.0, tmp_n_dminus2_2, 1, Qinv2x_tilde, 1);
          }

          for (size_t i = 0; i < n; i++) {
            id3 = i * d_minus_2;
            d_velocity_1[id3] = Qinv2x_bar[id3];
            d_velocity_2[id3] = Qinv2x_tilde[id3];

            d_t[i] = d_velocity_1[id3];
            d_t_prime[i] = d_velocity_2[id3];
          }

          break;
        }  // end of SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH

        default: {
          printf("ERROR\n");
        }

      }  // end switch

      // 9. Compute again the affine step-length
      alpha_primal_1 = grfc3d_getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
      alpha_primal_2 = grfc3d_getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
      alpha_dual_1 = grfc3d_getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
      alpha_dual_2 = grfc3d_getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

      alpha_primal =
          fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
      alpha_dual = alpha_primal;

      /* updating the gamma parameter used to compute the step-length */
      gmm = gmmp1 + gmmp2 * alpha_primal;

      // Print out all useful parameters
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        numerics_printf_verbose(
            -1,
            "| %2i %3d| %7.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e "
            "| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e %.1e %.1e | %.1e %.1e %.1e | %.1e |",
            iteration, nRefine, relgap, pinfeas, dinfeas, udotr, proj_error, complem_1,
            complem_2, full_error, barr_param, alpha_primal, alpha_dual, sigma,
            cblas_dnrm2(m, d_globalVelocity, 1) / cblas_dnrm2(m, globalVelocity, 1),
            cblas_dnrm2(nd, d_velocity, 1) / cblas_dnrm2(nd, velocity, 1),
            cblas_dnrm2(nd, d_reaction, 1) / cblas_dnrm2(nd, reaction, 1), residu_LS1_m,
            residu_LS1_nd, residu_LS1_ndplus1, residu_LS2_m, residu_LS2_nd, residu_LS2_ndplus1,
            residu_refine);
      } else {
        numerics_printf_verbose(
            -1,
            "| %2i| %7.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | "
            "%.1e | %.1e | %.1e | %.1e %.1e %.1e |",
            iteration, relgap, pinfeas, dinfeas, udotr, proj_error, complem_1, complem_2,
            barr_param, alpha_primal, sigma,
            cblas_dnrm2(m, d_globalVelocity, 1) / cblas_dnrm2(m, globalVelocity, 1),
            cblas_dnrm2(nd, d_velocity, 1) / cblas_dnrm2(nd, velocity, 1),
            cblas_dnrm2(nd, d_reaction, 1) / cblas_dnrm2(nd, reaction, 1), residu_LS2_m,
            residu_LS2_nd, residu_LS2_ndplus1);
      }

      // 10. Update variables
      cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
      cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);
      cblas_daxpy(n, alpha_dual, d_t, 1, t, 1);
      cblas_daxpy(n, alpha_dual, d_t_prime, 1, t_prime, 1);

      if (block_1) {
        block_1 = NM_free(block_1);
      }
      if (block_2) {
        block_2 = NM_free(block_2);
      }
      if (arrowMat_u1) {
        arrowMat_u1 = NM_free(arrowMat_u1);
      }
      if (arrowMat_u2) {
        arrowMat_u2 = NM_free(arrowMat_u2);
      }
      if (arrowMat_r1) {
        arrowMat_r1 = NM_free(arrowMat_r1);
      }
      if (arrowMat_r2) {
        arrowMat_r2 = NM_free(arrowMat_r2);
      }
      if (Z) {
        Z = NM_free(Z);
      }
      if (ZJT) {
        ZJT = NM_free(ZJT);
      }

      if (Qp_bar) {
        Qp_bar = NM_free(Qp_bar);
      }
      if (Qp_tilde) {
        Qp_tilde = NM_free(Qp_tilde);
      }
      if (Qpinv_bar) {
        Qpinv_bar = NM_free(Qpinv_bar);
      }
      if (Qpinv_tilde) {
        Qpinv_tilde = NM_free(Qpinv_tilde);
      }
      if (Qp2_bar) {
        Qp2_bar = NM_free(Qp2_bar);
      }
      if (Qp2_tilde) {
        Qp2_tilde = NM_free(Qp2_tilde);
      }
      if (Qpinv2_bar) {
        Qpinv2_bar = NM_free(Qpinv2_bar);
      }
      if (Qpinv2_tilde) {
        Qpinv2_tilde = NM_free(Qpinv2_tilde);
      }

      if (Qinv) {
        Qinv = NM_free(Qinv);
      }
      if (Qinv2) {
        Qinv2 = NM_free(Qinv2);
      }
      if (JQinv) {
        JQinv = NM_free(JQinv);
      }
      if (JQinvT) {
        JQinvT = NM_free(JQinvT);
      }
      if (JQinv2) {
        JQinv2 = NM_free(JQinv2);
      }
      if (JQJ) {
        JQJ = NM_free(JQJ);
      }

      if (P_inv) {
        P_inv = NM_free(P_inv);
      }
      if (P_invT) {
        P_invT = NM_free(P_invT);
      }
      if (HMHP) HMHP = NM_free(HMHP);
      if (PHMHP) PHMHP = NM_free(PHMHP);
      if (P_inv_F) {
        P_inv_F = NM_free(P_inv_F);
      }
      if (PinvH) {
        PinvH = NM_free(PinvH);
      }
      if (PinvH_T) {
        PinvH_T = NM_free(PinvH_T);
      }

      if (chol_U) chol_U = NM_free(chol_U);
      if (chol_L) chol_L = cs_spfree(chol_L);
      // if (chol_U_csc) chol_U_csc = cs_spfree(chol_U_csc); // already by NM_free(chol_U);
      if (chol_UT_csc) chol_UT_csc = cs_spfree(chol_UT_csc);

      if (Jac) Jac = NM_free(Jac);

      if (jacobian_is_nan | NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) |
          NV_isnan(reaction, nd)) {
        hasNotConverged = 2;
        break;
      }

      iteration++;
    }  // end of while loop

    if (hasNotConverged && options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
                               SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_AFTER) {
      options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] =
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES;
      reinit = 1;
    } else
      break;

  }  // end while (refinement_after)

  /* -------------------------- Return to original variables -------------------------- */
  // TO DO: If we need this point, plz uncomment
  // NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->original_point->velocity);
  // NM_gemv(1.0, P_mu, reaction, 0.0, data->original_point->reaction);
  // cblas_dcopy(m, globalVelocity, 1, data->original_point->globalVelocity, 1);

  // options->dparam[SICONOS_DPARAM_RESIDU] = full_error; //NV_max(error, 4);
  options->dparam[SICONOS_DPARAM_RESIDU] =
      fmax(pinfeas, fmax(dinfeas, fmin(udotr, proj_error)));
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;

  clock_t t2 = clock();
  long clk_tck = CLOCKS_PER_SEC;

  /* writing data in a Matlab file */
  // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  // {
  // char matlab_file_name[256];
  // sprintf(matlab_file_name,"sigma_nc-%d-.m", problem->numberOfContacts);
  // matlab_file = fopen(matlab_file_name, "w");
  // printInteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d, n, m,
  // (double)(t2-t1)/(double)clk_tck, matlab_file); fclose(matlab_file);
  // }

  if (iden) {
    free(iden);
    iden = NULL;
  }
  if (rhs_save) {
    free(rhs_save);
    rhs_save = NULL;
  }
  if (rhs_dx) {
    free(rhs_dx);
    rhs_dx = NULL;
  }

  if (Minv) {
    Minv = NM_free(Minv);
  }
  if (HMinv) {
    HMinv = NM_free(HMinv);
  }
  if (HMinvHt) {
    HMinvHt = NM_free(HMinvHt);
  }
  if (minus_M) {
    minus_M = NM_free(minus_M);
  }
  if (minus_H) {
    minus_H = NM_free(minus_H);
  }
  if (minus_Ht) {
    minus_Ht = NM_free(minus_Ht);
  }
  if (Ht) {
    Ht = NM_free(Ht);
  }
  if (H_origin) {
    H_origin = NM_free(H_origin);
  }
  if (H) {
    H = NM_free(H);
  }
  if (J) {
    J = NM_free(J);
  }
  if (Jt) {
    Jt = NM_free(Jt);
  }
  if (identity) {
    identity = NM_free(identity);
  }

  if (internal_allocation) {
    grfc3d_IPM_free(problem, options);
  }

  options->solverData = (double *)malloc(sizeof(double));
  double *projerr_ptr = (double *)options->solverData;
  *projerr_ptr = proj_error;

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) fclose(iterates);

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE]) {
    fprintf(iterates_python, "it = %lu\n", iteration + 1);
    fclose(iterates_python);
  }
  *info = hasNotConverged;

}  // end of grfc3d_IPM

/* initialize solver (allocate memory) */
void grfc3d_IPM_init(GlobalRollingFrictionContactProblem *problem, SolverOptions *options) {
  size_t m = problem->M->size0;
  size_t nd = problem->H->size1;
  size_t d = problem->dimension;  // d must be 5 because of rolling friction problem
  size_t n = problem->numberOfContacts;

  size_t n_dminus2 = n * (d - 2);
  size_t n_dplus1 = n * (d + 1);

  if (!options->dWork) {
    switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv:
        options->dWork = (double *)calloc(m + nd + n_dplus1, sizeof(double));
        options->dWorkSize = m + nd + n_dplus1;
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ:
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH:
        options->dWork = (double *)calloc(m + nd, sizeof(double));
        options->dWorkSize = m + nd;
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ:
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
        options->dWork = (double *)calloc(nd, sizeof(double));
        options->dWorkSize = nd;
        break;

      default:
        printf("ERROR of options->dWork allocation in grfc3d_IPM_init\n");
    }
  }

  /* ------------- initialize starting point ------------- */
  options->solverData = (Grfc3d_IPM_data *)malloc(sizeof(Grfc3d_IPM_data));
  Grfc3d_IPM_data *data = (Grfc3d_IPM_data *)options->solverData;

  /* --------- allocate memory for IPM point ----------- */
  data->starting_point = (IPM_point *)malloc(sizeof(IPM_point));

  /* 1. v */
  data->starting_point->globalVelocity = (double *)calloc(m, sizeof(double));
  for (size_t i = 0; i < m; ++i) data->starting_point->globalVelocity[i] = 0.01;

  /* 2. u */
  data->starting_point->velocity = (double *)calloc(nd, sizeof(double));
  for (size_t i = 0; i < nd; ++i) {
    data->starting_point->velocity[i] = 0.001;
    if (i % d == 0) data->starting_point->velocity[i] = 3.0;
  }

  /* 3. r */
  data->starting_point->reaction = (double *)calloc(nd, sizeof(double));
  for (size_t i = 0; i < nd; ++i) {
    data->starting_point->reaction[i] = 0.04;
    if (i % d == 0) data->starting_point->reaction[i] = 0.5;
  }

  /* original point which is not changed by the matrix P_mu */
  data->original_point = (IPM_point *)malloc(sizeof(IPM_point));
  data->original_point->globalVelocity = (double *)calloc(m, sizeof(double));
  data->original_point->velocity = (double *)calloc(nd, sizeof(double));
  data->original_point->reaction = (double *)calloc(nd, sizeof(double));

  /* ------ initialize the change of variable matrix P_mu ------- */
  data->P_mu = (IPM_change_of_variable *)malloc(sizeof(IPM_change_of_variable));
  data->P_mu->mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->mat, nd);
  data->P_mu->mat->matrix2->origin = NSM_TRIPLET;
  for (size_t i = 0; i < nd; ++i) {
    if (i % d == 0)
      NM_entry(data->P_mu->mat, i, i, 1.);
    else if (i % d == 1 || i % d == 2)
      NM_entry(data->P_mu->mat, i, i, problem->mu[(int)(i / d)]);
    else
      NM_entry(data->P_mu->mat, i, i, problem->mu_r[(int)(i / d)]);
    // NM_entry(data->P_mu->mat, i, i, problem->mu[(int)(i/d)]/1000);
  }

  /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
  data->P_mu->inv_mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->inv_mat, nd);
  data->P_mu->inv_mat->matrix2->origin = NSM_TRIPLET;
  for (size_t i = 0; i < nd; ++i) {
    if (i % d == 0)
      NM_entry(data->P_mu->inv_mat, i, i, 1.);
    else if (i % d == 1 || i % d == 2)
      NM_entry(data->P_mu->inv_mat, i, i, 1.0 / problem->mu[(int)(i / d)]);
    else
      NM_entry(data->P_mu->inv_mat, i, i, 1.0 / problem->mu_r[(int)(i / d)]);
    // NM_entry(data->P_mu->inv_mat, i, i, 1.0/(problem->mu[(int)(i/d)])/1000);
  }

  /* ------ initial parameters initialization ---------- */
  data->internal_params = (IPM_internal_params *)malloc(sizeof(IPM_internal_params));
  data->internal_params->alpha_primal = 1.0;
  data->internal_params->alpha_dual = 1.0;
  data->internal_params->sigma = 0.1;
  data->internal_params->barr_param = 1.0;

  /* ----- temporary vaults initialization ------- */
  data->tmp_vault_m = (double **)malloc(5 * sizeof(double *));
  for (size_t i = 0; i < 5; ++i) data->tmp_vault_m[i] = (double *)calloc(m, sizeof(double));

  data->tmp_vault_n = (double **)malloc(5 * sizeof(double *));
  for (size_t i = 0; i < 5; ++i) data->tmp_vault_n[i] = (double *)calloc(n, sizeof(double));

  data->tmp_vault_n_dminus2 = (double **)malloc(30 * sizeof(double *));
  for (size_t i = 0; i < 30; ++i)
    data->tmp_vault_n_dminus2[i] = (double *)calloc(n_dminus2, sizeof(double));

  data->tmp_vault_nd = (double **)malloc(10 * sizeof(double *));
  for (size_t i = 0; i < 10; ++i) data->tmp_vault_nd[i] = (double *)calloc(nd, sizeof(double));

  data->tmp_vault_n_dplus1 = (double **)malloc(2 * sizeof(double *));
  for (size_t i = 0; i < 2; ++i)
    data->tmp_vault_n_dplus1[i] = (double *)calloc(n_dplus1, sizeof(double));

}  // end of grfc3d_IPM_init

/* deallocate memory */
void grfc3d_IPM_free(GlobalRollingFrictionContactProblem *problem, SolverOptions *options) {
  if (options->dWork) {
    free(options->dWork);
    options->dWork = NULL;
    options->dWorkSize = 0;
  }

  if (options->solverData) {
    Grfc3d_IPM_data *data = (Grfc3d_IPM_data *)options->solverData;

    free(data->starting_point->globalVelocity);
    data->starting_point->globalVelocity = NULL;

    free(data->starting_point->velocity);
    data->starting_point->velocity = NULL;

    free(data->starting_point->reaction);
    data->starting_point->reaction = NULL;

    free(data->starting_point);
    data->starting_point = NULL;

    free(data->original_point->globalVelocity);
    data->original_point->globalVelocity = NULL;

    free(data->original_point->velocity);
    data->original_point->velocity = NULL;

    free(data->original_point->reaction);
    data->original_point->reaction = NULL;

    free(data->original_point);
    data->original_point = NULL;

    NM_clear(data->P_mu->mat);
    free(data->P_mu->mat);
    data->P_mu->mat = NULL;

    NM_clear(data->P_mu->inv_mat);
    free(data->P_mu->inv_mat);
    data->P_mu->inv_mat = NULL;

    free(data->P_mu);
    data->P_mu = NULL;

    free(data->internal_params);
    data->internal_params = NULL;

    for (size_t i = 0; i < 5; ++i) free(data->tmp_vault_m[i]);
    free(data->tmp_vault_m);
    data->tmp_vault_m = NULL;

    for (size_t i = 0; i < 5; ++i) free(data->tmp_vault_n[i]);
    free(data->tmp_vault_n);
    data->tmp_vault_n = NULL;

    for (size_t i = 0; i < 30; ++i) free(data->tmp_vault_n_dminus2[i]);
    free(data->tmp_vault_n_dminus2);
    data->tmp_vault_n_dminus2 = NULL;

    for (size_t i = 0; i < 10; ++i) free(data->tmp_vault_nd[i]);
    free(data->tmp_vault_nd);
    data->tmp_vault_nd = NULL;

    for (size_t i = 0; i < 2; ++i) free(data->tmp_vault_n_dplus1[i]);
    free(data->tmp_vault_n_dplus1);
    data->tmp_vault_n_dplus1 = NULL;
  }

  free(options->solverData);
  options->solverData = NULL;
}  // end of grfc3d_IPM_free

/* setup default solver parameters */
void grfc3d_IPM_set_default(SolverOptions *options) {
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  /* 0: convex case;  1: non-smooth case */
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = 0;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
      SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_AFTER;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] =
      SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY] =
      SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] =
      SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] =
  // SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE] = 0;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-5;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09;  // 0.09

}  // end of grfc3d_IPM_set_default
