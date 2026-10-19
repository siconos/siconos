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

/* solver based on Interior Point Method (IPM) for friction-contact 3D problem based on an AVI
 * reformulation Vincent Acary, Paul Armand, Hoang Minh NGUYEN. High-accuracy computation of
 * rolling friction contact problems. 2022. https://hal.inria.fr/hal-03741048
 */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "CSparseMatrix.h"
#include "FrictionContact_options.h"
#include "JordanAlgebra.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "NumericsVector.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "fc3d_short_names.h"
#include "float.h"
#include "gfc3d_Solvers.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "projectionOnCone.h"
#include "safe_casts.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include <time.h>

#include "gfc3d_ipm.h"
#include "siconos_debug.h"

/* ------------------------- Helper functions implementation ------------------------------ */

/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double.
 */
double getStepLength(const double* const x, const double* const dx, const size_t vecSize,
                     const size_t varsCount, const double gamma) {
  size_t dimension = (int)(vecSize / varsCount);
  size_t pos;
  float_type aL, bL, cL, dL, alphaL;
  double min_alpha;

  min_alpha = 1e20;  // 1.0;

  for (size_t i = 0; i < varsCount; ++i) {
    pos = i * dimension;
    aL = dnrm2l(dimension - 1, dx + pos + 1);
    aL = (dx[pos] - aL) * (dx[pos] + aL);
    bL = x[pos] * dx[pos];
    for (size_t k = 1; k < dimension; bL -= x[pos + k] * dx[pos + k], k++);
    cL = dnrm2l(dimension - 1, x + pos + 1);
    cL = (x[pos] - cL) * (x[pos] + cL);
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

/* Returns the primal constraint vector for global fricprob: out = velocity - H x
 * globalVelocity - w */
/* and the relative 2-norm of this vector: |out|/max{|velocity|, |H x globalVelocity|, |w|} */
void primalResidual(const double* velocity, NumericsMatrix* H, const double* globalVelocity,
                    const double* w, double* out, double* rnorm, const double tol) {
  size_t nd = H->size0;
  double rn;

  /* The memory for the result vectors should be allocated using calloc
   * since H is a sparse matrix. In other case the behaviour will be undefined.*/
  //  double *Hv = (double*)calloc(nd, sizeof(double));
  // double *u_minus_Hv = (double*)calloc(nd, sizeof(double));

  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  rn = cblas_dnrm2(nd, out, 1);
  cblas_daxpy(nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(nd, -1.0, w, 1, out, 1);
  rn = fmax(rn, cblas_dnrm2(nd, velocity, 1));
  rn = fmax(rn, cblas_dnrm2(nd, w, 1));
  *rnorm = (rn > tol ? cblas_dnrm2(nd, out, 1) / rn : cblas_dnrm2(nd, out, 1));

  /* *rnorm = cblas_dnrm2(nd, out, 1);  */
  // printf("rn = %e, tol = %e\n", rn, tol);
}

void primalResidual_type(const double* velocity, NumericsMatrix* H,
                         const double* globalVelocity, const double* w, double* out,
                         double* rnorm, const double tol, const int type) {
  size_t nd = H->size0;
  double rn;

  /* The memory for the result vectors should be allocated using calloc
   * since H is a sparse matrix. In other case the behaviour will be undefined.*/
  //  double *Hv = (double*)calloc(nd, sizeof(double));
  // double *u_minus_Hv = (double*)calloc(nd, sizeof(double));

  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  rn = NV_norm_type(nd, out, type);
  cblas_daxpy(nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(nd, -1.0, w, 1, out, 1);
  rn = fmax(rn, NV_norm_type(nd, velocity, type));
  rn = fmax(rn, NV_norm_type(nd, w, type));
  *rnorm = (rn > tol ? NV_norm_type(nd, out, type) : NV_norm_type(nd, out, type));
}

/* Returns the dual constraint vector for global fricprob ( M*globalVelocity - f - H'*reaction
 * ) */
void dualResidual(NumericsMatrix* M, const double* globalVelocity, NumericsMatrix* H,
                  const double* reaction, const double* f, double* out, double* rnorm,
                  const double tol) {
  double m = H->size1;
  double* HTr = (double*)calloc(m, sizeof(double));
  double rn;

  NM_gemv(1.0, M, globalVelocity, 0.0, out);
  rn = cblas_dnrm2(m, out, 1);
  cblas_daxpy(m, -1.0, f, 1, out, 1);
  NM_tgemv(1.0, H, reaction, 0.0, HTr);
  cblas_daxpy(m, -1.0, HTr, 1, out, 1);
  rn = fmax(rn, cblas_dnrm2(m, f, 1));
  rn = fmax(rn, cblas_dnrm2(m, HTr, 1));
  *rnorm = (rn > tol ? cblas_dnrm2(m, out, 1) / rn : cblas_dnrm2(m, out, 1));
  /* *rnorm = cblas_dnrm2(m, out, 1); */
  free(HTr);
}

void dualResidual_type(NumericsMatrix* M, const double* globalVelocity, NumericsMatrix* H,
                       const double* reaction, const double* f, double* out, double* rnorm,
                       const double tol, const int type) {
  double m = H->size1;
  double* HTr = (double*)calloc(m, sizeof(double));
  double rn;

  NM_gemv(1.0, M, globalVelocity, 0.0, out);
  rn = NV_norm_type(m, out, type);
  cblas_daxpy(m, -1.0, f, 1, out, 1);
  NM_tgemv(1.0, H, reaction, 0.0, HTr);
  cblas_daxpy(m, -1.0, HTr, 1, out, 1);
  rn = fmax(rn, NV_norm_type(m, f, type));
  rn = fmax(rn, NV_norm_type(m, HTr, type));
  *rnorm = (rn > tol ? NV_norm_type(m, out, type) : NV_norm_type(m, out, type));
  free(HTr);
}

double NV_norm_type(const size_t vecSize, const double* const vec, const int type) {
  double norm = -1;

  if (type == NORM_2)  // L-2 norm
  {
    norm = cblas_dnrm2(vecSize, vec, 1);
  }

  else if (type == NORM_INF)  // L-inf norm
  {
    int maxIndex = cblas_idamax(vecSize, vec, 1);
    norm = fabs(vec[maxIndex]);
  }

  else {
    fprintf(stderr, "NV_norm_type: type = %d is undefined.\n", type);
    return -1;
  }

  return norm;
}

double xdoty_type(const size_t varsCount, const size_t vecSize, const double* x,
                  const double* y, const int type) {
  double xdoty = -1;

  if (type == NORM_2)
    xdoty = cblas_ddot(vecSize, x, 1, y, 1);

  else if (type == NORM_INF)
    for (size_t i = 0; i < varsCount; i++) {
      xdoty = fmax(xdoty, fabs(cblas_ddot(3, x + i * 3, 1, y + i * 3, 1)));
    }

  else if (type == NORM_2_INF) {
    double* xioyi = (double*)calloc(3, sizeof(double));
    for (size_t i = 0; i < varsCount; i++) {
      JA_prod(x + i * 3, y + i * 3, 3, 1, xioyi);
      xdoty = fmax(xdoty, cblas_dnrm2(3, xioyi, 1));
    }
    free(xioyi);
  }

  else {
    fprintf(stderr, "xdoty_type: type = %d is undefined.\n", type);
    return -1;
  }

  return xdoty;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product
 * velocity o reaction  */
double complemResidualNorm(const double* const velocity, const double* const reaction,
                           const size_t vecSize, const size_t varsCount) {
  double* resid = (double*)calloc(vecSize, sizeof(double));
  JA_prod(velocity, reaction, vecSize, varsCount, resid);
  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  free(resid);
  return norm2;
}

/* Returns the type-norm of the complementarity residual vector = type-norm of the Jordan
 * product velocity o reaction  */
double complemResidualNorm_type(const double* const velocity, const double* const reaction,
                                const size_t vecSize, const size_t varsCount, const int type) {
  double* resid = NULL;
  double norm = -1;
  size_t d = vecSize / varsCount;

  if (type == NORM_2) {
    resid = (double*)calloc(vecSize, sizeof(double));
    JA_prod(velocity, reaction, vecSize, varsCount, resid);
    norm = cblas_dnrm2(vecSize, resid, 1);
  }

  else if (type == NORM_2_INF) {
    resid = (double*)calloc(d, sizeof(double));

    for (size_t i = 0; i < vecSize; i += d) {
      JA_prod(velocity + i, reaction + i, d, 1, resid);
      norm = fmax(norm, cblas_dnrm2(d, resid, 1));
    }
  }

  else {
    numerics_error_log("complemResidualNorm_type", "unknown norm type");
    return NAN;
  }

  if (resid) free(resid);
  return norm;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product
 * (Qp*velocity) o (Qp_inv * reaction)  */
double complemResidualNorm_p(const double* const velocity, const double* const reaction,
                             const size_t vecSize, const size_t varsCount) {
  double* resid = (double*)calloc(vecSize, sizeof(double));
  double* u_p = (double*)calloc(vecSize, sizeof(double));
  double* r_p = (double*)calloc(vecSize, sizeof(double));
  double* a = (double*)calloc(vecSize, sizeof(double));
  double* b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(velocity, reaction, vecSize, varsCount, a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(velocity, b, vecSize, varsCount, a);

  Qx50y(a, velocity, vecSize, varsCount, u_p);
  Qx05y(a, reaction, vecSize, varsCount, r_p);
  JA_prod(u_p, r_p, vecSize, varsCount, resid);

  double norm2 = cblas_dnrm2(vecSize, resid, 1);

  free(resid);
  free(u_p);
  free(r_p);
  free(a);
  free(b);

  return norm2;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product
 * (Qp*velocity) o (Qp_inv * reaction)  */
/* This computation is done with the formula "F" */
double complemResidualNorm_p_F(NumericsMatrix* Qp, NumericsMatrix* Qpinv,
                               const double* const velocity, const double* const reaction,
                               const size_t vecSize, const size_t varsCount) {
  double* resid = (double*)calloc(vecSize, sizeof(double));
  double* u_p = (double*)calloc(vecSize, sizeof(double));
  double* r_p = (double*)calloc(vecSize, sizeof(double));

  NM_gemv(1.0, Qp, velocity, 0.0, u_p);
  NM_gemv(1.0, Qpinv, reaction, 0.0, r_p);
  JA_prod(u_p, r_p, vecSize, varsCount, resid);
  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  free(resid);
  free(u_p);
  free(r_p);
  return norm2;
}

/* computation of the duality gap  */
/* dualgap = gapVal / (1 + (abs(primal value) + abs(dual value))/2) */
double dualGap(NumericsMatrix* M, const double* f, const double* w,
               const double* globalVelocity, const double* reaction, const size_t nd,
               const size_t m) {
  double* Mv = (double*)calloc(m, sizeof(double));
  double vMv, pval, dval;

  NM_gemv(0.5, M, globalVelocity, 0.0, Mv);
  vMv = cblas_ddot(m, globalVelocity, 1, Mv, 1);
  free(Mv);
  pval = vMv - cblas_ddot(m, f, 1, globalVelocity, 1);
  dval = -vMv - cblas_ddot(nd, w, 1, reaction, 1);
  return (pval - dval) / (1 + (fabs(pval) + fabs(dval)) / 2);
}

/* Rel gap = gapVal / (1 + abs(primal value) + abs(dual value)) */
double relGap(NumericsMatrix* M, const double* f, const double* w,
              const double* globalVelocity, const double* reaction, const size_t nd,
              const size_t m, const double gapVal) {
  double* Mv = (double*)calloc(m, sizeof(double));
  double vMv, pval, dval;

  NM_gemv(0.5, M, globalVelocity, 0.0, Mv);
  vMv = cblas_ddot(m, globalVelocity, 1, Mv, 1);
  free(Mv);
  pval = vMv - cblas_ddot(m, f, 1, globalVelocity, 1);
  dval = -vMv - cblas_ddot(nd, w, 1, reaction, 1);
  return gapVal / (1 + fabs(pval) + fabs(dval));
}

/* Computation of the projection error |r - proj(r-u)|/max{|r|, |u|} */
double projectionError(const double* velocity, const double* reaction, const size_t nc,
                       const double tol) {
  double worktmp[3];
  double out = 0.0;
  double norm_u, norm_r, relative_scaling;

  for (size_t ic = 0; ic < nc; ic++) {
    worktmp[0] = reaction[3 * ic] - velocity[3 * ic];
    worktmp[1] = reaction[3 * ic + 1] - velocity[3 * ic + 1];
    worktmp[2] = reaction[3 * ic + 2] - velocity[3 * ic + 2];
    projectionOnCone(worktmp, 1.0);
    worktmp[0] = reaction[3 * ic] - worktmp[0];
    worktmp[1] = reaction[3 * ic + 1] - worktmp[1];
    worktmp[2] = reaction[3 * ic + 2] - worktmp[2];
    out += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }
  out = sqrt(out);
  norm_u = cblas_dnrm2(3 * nc, velocity, 1);
  norm_r = cblas_dnrm2(3 * nc, reaction, 1);
  relative_scaling = fmax(norm_u, norm_r);
  relative_scaling = 1.;
  // if(relative_scaling > tol)
  if (relative_scaling > DBL_EPSILON) out = out / relative_scaling;
  return out;
}

/* Computation of the projection error |r - proj(r-u)|_inf */
double projectionError_norm_infinity_conic(const double* velocity, const double* reaction,
                                           const size_t nc) {
  double worktmp[3];
  double out = -1.0;

  for (size_t ic = 0; ic < nc; ic++) {
    worktmp[0] = reaction[3 * ic] - velocity[3 * ic];
    worktmp[1] = reaction[3 * ic + 1] - velocity[3 * ic + 1];
    worktmp[2] = reaction[3 * ic + 2] - velocity[3 * ic + 2];
    projectionOnCone(worktmp, 1.0);
    worktmp[0] = reaction[3 * ic] - worktmp[0];
    worktmp[1] = reaction[3 * ic + 1] - worktmp[1];
    worktmp[2] = reaction[3 * ic + 2] - worktmp[2];
    out =
        fmax(out, worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
  }

  return sqrt(out);
}

/* Computation of the projection error max { |ri - proj(ri-ui)| / max { |ri|, |ui| } } */
double projectionError_relative_norm_infinity_conic(const double* velocity,
                                                    const double* reaction, const size_t nc) {
  double worktmp[3];
  double out = -1.0;
  double nr = -1.0, nu = -1.0, error;
  int idx = 0;

  for (size_t ic = 0; ic < nc; ic++) {
    idx = 3 * ic;

    nr = cblas_dnrm2(3, reaction + idx, 1);
    nu = cblas_dnrm2(3, velocity + idx, 1);

    worktmp[0] = reaction[idx] - velocity[idx];
    worktmp[1] = reaction[idx + 1] - velocity[idx + 1];
    worktmp[2] = reaction[idx + 2] - velocity[idx + 2];
    projectionOnCone(worktmp, 1.0);
    worktmp[0] = reaction[idx] - worktmp[0];
    worktmp[1] = reaction[idx + 1] - worktmp[1];
    worktmp[2] = reaction[idx + 2] - worktmp[2];

    error = cblas_dnrm2(3, worktmp, 1) / fmax(nr, nu);

    out = fmax(out, error);
  }

  return out;
}

/* Computation of the projection error on dual cone |u - proj(u-r)|_inf */
double projectionError_dual_norm_infinity_conic(const double* velocity, const double* reaction,
                                                const size_t nc) {
  double worktmp[3];
  double out = 0.0;

  for (size_t ic = 0; ic < nc; ic++) {
    worktmp[0] = velocity[3 * ic] - reaction[3 * ic];
    worktmp[1] = velocity[3 * ic + 1] - reaction[3 * ic + 1];
    worktmp[2] = velocity[3 * ic + 2] - reaction[3 * ic + 2];
    projectionOnDualCone(worktmp, 1.0);
    worktmp[0] = velocity[3 * ic] - worktmp[0];
    worktmp[1] = velocity[3 * ic + 1] - worktmp[1];
    worktmp[2] = velocity[3 * ic + 2] - worktmp[2];
    out =
        fmax(out, worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
  }

  return sqrt(out);
}

/* Computation of the projection error with recomputation of u and v */
double projectionError_based_reaction_norm_infinity_conic(NumericsMatrix* H, NumericsMatrix* M,
                                                          const double* f, const double* w,
                                                          const double* reaction,
                                                          const size_t varsCount,
                                                          const int on_dual_cone) {
  size_t m = M->size0;
  size_t d = 3;
  size_t nd = varsCount * d;

  double* globalVelocity = (double*)calloc(m, sizeof(double));
  double* velocity = (double*)calloc(nd, sizeof(double));

  // Re-compute v = M\(H'*r+f)
  cblas_dcopy(m, f, 1, globalVelocity, 1);
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  int info_cholesky = NM_Cholesky_solve(NM_preserve(M), globalVelocity, 1);
  if (info_cholesky) {
    numerics_error_log("projectionError_based_reaction_norm_infinity_conic",
                       "NM_Cholesky_solve failed, we switch to LU");
    int info_lu = NM_LU_solve(NM_preserve(M), globalVelocity, 1);
    if (info_lu) {
      numerics_error_log("projectionError_based_reaction_norm_infinity_conic",
                         "NM_LU_solve failed");
    }
  }
  // Re-compute u = Hv + w + phi(Hv + w)
  cblas_dcopy(nd, w, 1, velocity, 1);
  NM_gemv(1, H, globalVelocity, 1.0, velocity);
  for (size_t i = 0; i < nd; i += d) {
    velocity[i] += cblas_dnrm2(2, velocity + i + 1, 1);
  }

  double out = -1;

  if (on_dual_cone) {
    out = projectionError_dual_norm_infinity_conic(velocity, reaction, varsCount);
  } else {
    out = projectionError_norm_infinity_conic(velocity, reaction, varsCount);
  }

  free(globalVelocity);
  free(velocity);

  return out;
}

void setErrorArray(double* error, const double pinfeas, const double dinfeas,
                   const double udotr, const double dualgap, const double complem,
                   const double projerr) {
  error[0] = pinfeas;
  error[1] = dinfeas;
  error[2] = udotr;
  error[3] = dualgap;
  error[4] = complem;
  error[5] = projerr;
}

/* Return the 2-norm of the difference between two vectors */
double norm2VecDiff(const double* vec1, const double* vec2, const size_t vecSize) {
  double* vecDiff;
  double nvd;
  vecDiff = (double*)calloc(vecSize, sizeof(double));
  cblas_dcopy(vecSize, vec1, 1, vecDiff, 1);
  cblas_daxpy(vecSize, -1.0, vec2, 1, vecDiff, 1);
  nvd = cblas_dnrm2(vecSize, vecDiff, 1);
  free(vecDiff);
  return nvd;
}

/* Returns the product Q_{p}*H where p is the NT vector related to the pair (x,y) and H is the
 * matrix of the linear constraint. */
static NumericsMatrix* QNTpH(const double* const x, const double* const y, NumericsMatrix* H,
                             const size_t vecSize, const size_t varsCount) {
  double* a = (double*)calloc(vecSize, sizeof(double));
  double* b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(x, y, vecSize, varsCount, a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);

  NumericsMatrix* QpH = NM_new();

  // Qx50y(a, z, vecSize, varsCount, out);
  NM_types storage = H->storageType;
  switch (storage) {
      /* case NM_DENSE: */
      /*   cblas_dgemv(CblasColMajor, CblasNoTrans, sizeY, sizeX, alpha, A->dense, sizeY, x,
       * 1, beta, y, 1); */
      /*   break; */
      /* /\* SparseBlock storage *\/ */
      /* case NM_SPARSE_BLOCK: */
      /*   SBM_gemv_3x3(sizeX, sizeY, A->sparse_block, x, y); */
      /*   break; */
      /* coordinate */
    case NM_SPARSE: {
      CSparseMatrix* H_csc = NM_csc(H);
      CSparseMatrix* QpH_csc =
          cs_spalloc(H->size0, H->size1, H_csc->nzmax, 1, 0); /* allocate result */

      CS_INT i, p, *Hp, *Hi;
      Hp = H_csc->p;
      Hi = H_csc->i;

      CS_INT *QpHp, *QpHi;
      QpHp = QpH_csc->p;

      CS_ENTRY *Hx, *QpHx;
      Hx = H_csc->x;

      // size_t dimension = (int)(vecSize / varsCount);  // should be used to work
      // properly
      //                                                       // on different dimension of
      //                                                       cones.

      double* z_beta = b;
      for (CS_INT k = 0; k < H->size0; k++) z_beta[k] = 0;

      int* beta = (int*)malloc(H->size0 / 3 * sizeof(int));
      for (CS_INT k = 0; k < H->size0 / 3; k++) beta[k] = -1;

      CS_INT nz = 0;
      for (CS_INT k = 0; k < H->size1; k++) {
        /* search for beta and z_beta that are non null in the column k of H */
        CS_INT n_beta = -1, beta_old = -1;
        for (p = Hp[k]; p < Hp[k + 1]; p++) {
          i = Hi[p];
          z_beta[i] = Hx[p];
          CS_INT beta_current = i / 3;
          if (beta_old != beta_current) {
            n_beta++;
            beta_old = beta_current;
          }
          beta[n_beta] = beta_current;
        }
        n_beta++;

        QpHp[k] = nz; /* column k of QpH starts here */

        /* reallocate if needed */
        if ((nz + H->size0) > QpH_csc->nzmax &&
            !cs_sprealloc(QpH_csc, 2 * (QpH_csc->nzmax) + H->size0)) {
          free(beta);
          return NULL; /* out of memory */
        }
        QpHi = QpH_csc->i;
        QpHx = QpH_csc->x; /* C->i and C->x may be reallocated */

        /* multiplication and storage */
        for (int b = 0; b < n_beta; b++) {
          CS_INT alpha = beta[b];

          Qx50y(&a[alpha * 3], &z_beta[alpha * 3], 3, 1, QpHx + nz);

          /* store out in QpH */
          QpHi[nz++] = alpha * 3;
          QpHi[nz++] = alpha * 3 + 1;
          QpHi[nz++] = alpha * 3 + 2;

          z_beta[alpha * 3] = 0.0;
          z_beta[alpha * 3 + 1] = 0.0;
          z_beta[alpha * 3 + 2] = 0.0;
          beta[b] = -1;
        }
      }  // end loop k

      QpHp[H->size1] = nz;
      cs_sprealloc(QpH_csc, 0);

      QpH->storageType = H->storageType;
      numericsSparseMatrix(QpH)->csc = QpH_csc;
      QpH->size0 = (int)QpH->sparse->csc->m;
      QpH->size1 = (int)QpH->sparse->csc->n;
      numericsSparseMatrix(QpH)->origin = NSM_CSC;
      free(beta);
      break;
    } break;
    default:
      fprintf(stderr, "Numerics, GFC3D IPM, QNTpH failed, unknown storage type for H.\n");
      return NULL;
  }
  free(a);
  free(b);
  return QpH;
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
static void printDataProbMatlabFile(NumericsMatrix* M, double* f, NumericsMatrix* H, double* w,
                                    int d, int n, int m, double* mu, FILE* file) {
  fprintf(file, "d = %3i;\n", d);
  fprintf(file, "n = %6i;\n", n);
  fprintf(file, "m = %6i;\n", m);

  fprintf(file, "M = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(M), 0, file);
  fprintf(file, "];\n");
  fprintf(file, "M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");

  fprintf(file, "H = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, file);
  fprintf(file, "];\n");
  fprintf(file, "H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");

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
}

// /* print iteres under a Matlab format in a file */
// /* iteration = index of the iteration
//    v = global velocity
//    u = velocity
//    r = reaction
//    d = dimenion
//    n = number of contact points
//    m = number of degrees of freedom
// */
// static void printIteresProbMatlabFile(int iteration, double* v, double* u, double* r,
//                                       double* dv, double* du, double* dr, int d, int n, int
//                                       m, FILE* file) {
//   fprintf(file, "v(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < m; i++) {
//     fprintf(file, "%8.20e, ", v[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "u(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", u[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "r(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", r[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "dv(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < m; i++) {
//     fprintf(file, "%8.20e, ", dv[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "du(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", du[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "dr(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", dr[i]);
//   }
//   fprintf(file, "];\n");

//   // fprintf(file,"pinfeas(%3i) = %20.16e;\n",iteration+1,pinfeas);
//   // fprintf(file,"dinfeas(%3i) = %20.16e;\n",iteration+1,dinfeas);
//   // fprintf(file,"udotr(%3i) = %20.16e;\n",iteration+1,udotr);
//   // fprintf(file,"residu(%3i) = %20.16e;\n",iteration+1,totalresidual);
//   // fprintf(file,"phiu(%3i) = %20.16e;\n",iteration+1,norm_phiu);
//   // // fprintf(file,"alpha_diff_phi(%3i) = %20.16e;\n",iteration+1,alpha_diff_phi);
//   // fprintf(file,"alpha_u(%3i) = %20.16e;\n",iteration+1,alpha_diff_phi);

//   return;
// }

// static void printVectorMatlabFile(int iteration, double* vec, int vecSize, FILE* file) {
//   fprintf(file, "vector(%4i,:) = [", iteration + 1);
//   for (int i = 0; i < vecSize; i++) {
//     fprintf(file, "%24.16e, ", vec[i]);
//   }
//   fprintf(file, "];\n");
//   return;
// }

void classify_BNRT(const double* velocity, const double* reaction, const size_t vecSize,
                   const size_t varsCount, int* nB, int* nN, int* nR, int* nT) {
  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    assert(0);
  }

  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns;
  *nB = *nN = *nR = *nT = 0;
  for (size_t i = 0; i < varsCount; i++) {
    somme[0] = velocity[3 * i] + reaction[3 * i];
    somme[1] = velocity[3 * i + 1] + reaction[3 * i + 1];
    somme[2] = velocity[3 * i + 2] + reaction[3 * i + 2];
    dur[0] = velocity[3 * i] - reaction[3 * i];
    dur[1] = velocity[3 * i + 1] - reaction[3 * i + 1];
    dur[2] = velocity[3 * i + 2] - reaction[3 * i + 2];

    ns = somme[0] - cblas_dnrm2(2, somme + 1, 1);
    // ndur = dur[0] - cblas_dnrm2(2, dur + 1, 1);
    if (ns > cesp * cblas_dnrm2(3, somme, 1)) {
      if (dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        *nN += 1;
      else if (-dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        *nB += 1;
      else
        *nR += 1;
      // if (dur[0] >= cblas_dnrm2(2,dur+1,1))
      // {
      //   *nN +=1;
      //   printf("N:%d  ",i);
      // }
      // else if (-dur[0] >= cblas_dnrm2(2,dur+1,1))
      // {
      //   *nB +=1;
      //   printf("B:%d  ",i);
      // }
      // else
      // {
      //   *nR +=1;
      //   printf("R:%d  ",i);
      // }
    }
  }
  // printf("\n");
  *nT = varsCount - *nB - *nN - *nR;
}

/* Return the classification BNRT of the input: orignal u (uN;uT) and r (rN; rT)
 * We need first a change of velocity : u_tilde = (uN + mu*|uT|; mu*uT)
 * these u_tilde and r belong to friction cones
 */
void classify_BNRT_velocity_original(const double* mu, const double* velocity,
                                     const double* reaction, const size_t vecSize,
                                     const size_t varsCount, int* nB, int* nN, int* nR,
                                     int* nT) {
  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    assert(0);
  }

  double* velocity_no_mu =
      (double*)calloc(vecSize, sizeof(double));  // = (uN + mu*|uT|; mu*uT)
  double* reaction_no_mu =
      (double*)calloc(vecSize, sizeof(double));  // = (rN          ; rT/mu)

  for (size_t i = 0; i < vecSize; i += d) {
    velocity_no_mu[i] = velocity[i] + mu[i / d] * cblas_dnrm2(2, velocity + i + 1, 1);
    velocity_no_mu[i + 1] = mu[i / d] * velocity[i + 1];
    velocity_no_mu[i + 2] = mu[i / d] * velocity[i + 2];

    reaction_no_mu[i] = reaction[i];
    reaction_no_mu[i + 1] = reaction[i + 1] / mu[i / d];
    reaction_no_mu[i + 2] = reaction[i + 2] / mu[i / d];
  }

  classify_BNRT(velocity_no_mu, reaction_no_mu, vecSize, varsCount, nB, nN, nR, nT);

  free(velocity_no_mu);
  free(reaction_no_mu);
}

/* Return the classification BNRT of the input: modified u (uN + mu*|uT|; uT) and r (rN; rT)
 * We need first a change of velocity : u_tilde = (uN + mu*|uT|; mu*uT)
 * these u_tilde and r belong to friction cones
 */
void classify_BNRT_velocity_modified(const double* mu, const double* velocity,
                                     const double* reaction, const size_t vecSize,
                                     const size_t varsCount, int* nB, int* nN, int* nR,
                                     int* nT) {
  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    assert(0);
  }

  double* velocity_no_mu =
      (double*)calloc(vecSize, sizeof(double));  // = (uN + mu*|uT|; mu*uT)
  double* reaction_no_mu =
      (double*)calloc(vecSize, sizeof(double));  // = (rN          ; rT/mu)

  for (size_t i = 0; i < vecSize; i += d) {
    velocity_no_mu[i] = velocity[i];
    velocity_no_mu[i + 1] = mu[i / d] * velocity[i + 1];
    velocity_no_mu[i + 2] = mu[i / d] * velocity[i + 2];

    reaction_no_mu[i] = reaction[i];
    reaction_no_mu[i + 1] = reaction[i + 1] / mu[i / d];
    reaction_no_mu[i + 2] = reaction[i + 2] / mu[i / d];
  }

  classify_BNRT(velocity_no_mu, reaction_no_mu, vecSize, varsCount, nB, nN, nR, nT);

  free(velocity_no_mu);
  free(reaction_no_mu);
}

// By defaut, array setR is set to 0. If a cone is classified in R, the corresponding index in
// setR is set to 1.
void classify_indices_R(const double* velocity, const double* reaction, const size_t vecSize,
                        const size_t varsCount, int* nR, int* setR) {
  // Note that setR must be allocated varsCount memory space
  // Reset setR
  for (size_t i = 0; i < varsCount; i++) {
    setR[i] = 0;
  }
  *nR = 0;

  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    assert(0);
  }

  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns;

  for (size_t i = 0; i < varsCount; i++) {
    somme[0] = velocity[3 * i] + reaction[3 * i];
    somme[1] = velocity[3 * i + 1] + reaction[3 * i + 1];
    somme[2] = velocity[3 * i + 2] + reaction[3 * i + 2];
    dur[0] = velocity[3 * i] - reaction[3 * i];
    dur[1] = velocity[3 * i + 1] - reaction[3 * i + 1];
    dur[2] = velocity[3 * i + 2] - reaction[3 * i + 2];

    ns = somme[0] - cblas_dnrm2(2, somme + 1, 1);
    // ndur = dur[0] - cblas_dnrm2(2, dur + 1, 1);
    if (ns > cesp * cblas_dnrm2(3, somme, 1)) {
      if (dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        continue;
      else if (-dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        continue;
      else {
        setR[i] = 1;
        (*nR)++;
      }
    }
  }

  // for (int i = 0; i < varsCount; i++)
  // {
  //   printf("classify_indices_R: setR[%d] = %d\n", i, setR[i]);
  // }
}

// cl: classify
void printBlockVec(double* vec, int vecSize, int sizeBlock, int cl) {
  if (cl > 0) {
    for (int i = 0; i < vecSize; i++) {
      if (i % sizeBlock == 0) printf("\n%4i-%4i: ", i, i + sizeBlock - 1);
      printf(" %.8e,", vec[i]);
      if ((i + 1) % sizeBlock == 0)
        printf("\t\t x0 - |xb| = %.2e",
               vec[((i + 1) / sizeBlock - 1) * sizeBlock] -
                   cblas_dnrm2(sizeBlock - 1, vec + ((i + 1) / sizeBlock - 1) * sizeBlock + 1,
                               1));
    }
  }

  else {
    for (int i = 0; i < vecSize; i++) {
      if (i % sizeBlock == 0) printf("\n%4d-%4d: ", i, i + sizeBlock - 1);
      printf(" %.15e", vec[i]);
    }
  }
  return;
}

/* --------------------------- Interior-point method implementation
 * ------------------------------ */
/*
 * Implementation contains the following functions:
 *  - gfc3d_IPM_init - initialize solver (allocate memory)
 *  - gfc3d_IPM_free - deallocate memory
 *  - gfc3d_IPM_setDefaultSolverOptions - setup default solver parameters
 *  - gfc3d_IPM - optimization method
 */
void gfc3d_IPM_init(GlobalFrictionContactProblem* problem, SolverOptions* options,
                    size_t work_size) {
  size_t m = problem->M->size0;
  size_t nd = problem->H->size1;
  size_t d = problem->dimension;

  if (!options->dWork || options->dWorkSize < work_size) {
    double* p = (double*)realloc(options->dWork, work_size * sizeof(double));
    if (!p) {
      numerics_error_log("gfc3d init", "bad alloc");
    }
    options->dWork = p;
    options->dWorkSize = work_size;
  }

  /* ------------- initialize starting point ------------- */
  options->solverData = (Gfc3d_IPM_init_data*)malloc(sizeof(Gfc3d_IPM_init_data));
  Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)options->solverData;

  /* --------- allocate memory for tmp point ----------- */
  data->tmp_point = (IPM_tmp_point*)malloc(sizeof(IPM_tmp_point));
  data->tmp_point->t_globalVelocity = (double*)calloc(m, sizeof(double));
  data->tmp_point->t_velocity = (double*)calloc(nd, sizeof(double));
  data->tmp_point->t_reaction = (double*)calloc(nd, sizeof(double));

  /* 1. v */
  data->starting_point = (IPM_starting_point*)malloc(sizeof(IPM_starting_point));
  data->starting_point->globalVelocity = (double*)calloc(m, sizeof(double));
  for (size_t i = 0; i < m; ++i) data->starting_point->globalVelocity[i] = 0.01;

  /* 2. u */
  data->starting_point->velocity = (double*)calloc(nd, sizeof(double));
  for (size_t i = 0; i < nd; ++i) {
    data->starting_point->velocity[i] = 0.01;
    if (i % d == 0) data->starting_point->velocity[i] = 0.1;
  }

  /* 3. r */
  data->starting_point->reaction = (double*)calloc(nd, sizeof(double));
  for (size_t i = 0; i < nd; ++i) {
    data->starting_point->reaction[i] = 0.01;                 // 0.0351;
    if (i % d == 0) data->starting_point->reaction[i] = 0.1;  // 0.2056;
  }

  /* ------ initialize the change of variable matrix P_mu ------- */
  data->P_mu = (IPM_change_of_variable*)malloc(sizeof(IPM_change_of_variable));
  data->P_mu->mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->mat, nd);
  data->P_mu->mat->sparse->origin = NSM_TRIPLET;
  for (size_t i = 0; i < nd; ++i)
    if (i % d == 0) /* NM_entry(data->P_mu->mat, i, i, 1. / problem->mu[(int)(i/d)]); */
      NM_entry(data->P_mu->mat, i, i, 1.);
    else
      /* NM_entry(data->P_mu->mat, i, i, 1.); */
      NM_entry(data->P_mu->mat, i, i, problem->mu[(int)(i / d)]);

  /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
  data->P_mu->inv_mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->inv_mat, nd);
  data->P_mu->inv_mat->sparse->origin = NSM_TRIPLET;
  for (size_t i = 0; i < nd; ++i)
    if (i % d == 0) /* NM_entry(data->P_mu->inv_mat, i, i, problem->mu[(int)(i/d)]); */
      NM_entry(data->P_mu->inv_mat, i, i, 1.);
    else
      /* NM_entry(data->P_mu->inv_mat, i, i, 1.); */
      NM_entry(data->P_mu->inv_mat, i, i, 1.0 / problem->mu[(int)(i / d)]);
  /* ------ initial parameters initialization ---------- */
  data->internal_params = (IPM_internal_params*)malloc(sizeof(IPM_internal_params));
  data->internal_params->alpha_primal = 1.0;
  data->internal_params->alpha_dual = 1.0;
  data->internal_params->sigma = 0.1;
  data->internal_params->barr_param = 1.0;

  /* ----- temporary vaults initialization ------- */
  data->tmp_vault_nd = (double**)malloc(17 * sizeof(double*));
  for (size_t i = 0; i < 17; ++i) data->tmp_vault_nd[i] = (double*)calloc(nd, sizeof(double));

  data->tmp_vault_m = (double**)malloc(2 * sizeof(double*));
  for (size_t i = 0; i < 2; ++i) data->tmp_vault_m[i] = (double*)calloc(m, sizeof(double));
}

// /* check the solution of the linear system A*x = b  */
// /* return |A*x[p:p+n-1]-b|/|b[p:p+n-1]|            */
// static double relative_error_linear_system_solution(NumericsMatrix* const A, const double*
// x,
//                                                     const double* b, int b_size, int p,
//                                                     int n) {
//   double out;
//   double norm_b = cblas_dnrm2(n, b + p, 1);
//   double* rhs = (double*)calloc(b_size, sizeof(double));
//   cblas_dcopy(b_size, b, 1, rhs, 1);
//   NM_gemv(1.0, A, x, -1.0, rhs);
//   out = cblas_dnrm2(n, rhs + p, 1);
//   printf("out=%9.2e\n", out);
//   out = out / norm_b;
//   free(rhs);
//   return out;
// }

void gfc3d_IPM_free(GlobalFrictionContactProblem* problem, SolverOptions* options) {
  if (options->solverData) {
    Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)options->solverData;

    if (data->starting_point) {
      free(data->starting_point->globalVelocity);
      data->starting_point->globalVelocity = NULL;

      free(data->starting_point->velocity);
      data->starting_point->velocity = NULL;

      free(data->starting_point->reaction);
      data->starting_point->reaction = NULL;

      free(data->starting_point);
    }
    data->starting_point = NULL;

    if (data->P_mu) {
      data->P_mu->mat = NM_free(data->P_mu->mat);

      data->P_mu->inv_mat = NM_free(data->P_mu->inv_mat);

      free(data->P_mu);
    }
    data->P_mu = NULL;

    if (data->tmp_vault_nd) {
      for (size_t i = 0; i < 17; ++i) free(data->tmp_vault_nd[i]);
      free(data->tmp_vault_nd);
    }
    data->tmp_vault_nd = NULL;

    if (data->tmp_vault_m) {
      for (size_t i = 0; i < 2; ++i) free(data->tmp_vault_m[i]);
      free(data->tmp_vault_m);
    }
    data->tmp_vault_m = NULL;

    if (data->tmp_point) {
      free(data->tmp_point->t_globalVelocity);
      data->tmp_point->t_globalVelocity = NULL;

      free(data->tmp_point->t_velocity);
      data->tmp_point->t_velocity = NULL;

      free(data->tmp_point->t_reaction);
      data->tmp_point->t_reaction = NULL;

      free(data->tmp_point);
    }
    data->tmp_point = NULL;

    if (data->internal_params) free(data->internal_params);
    data->internal_params = NULL;

    free(options->solverData);
  }
  options->solverData = NULL;
}

/* Prepare the problem matrices: symmetrize M, optionally force sparse storage, and compute the
 * block-diagonal inverse Minv needed by the 1x1 QPH formulation. Outputs M (may point to
 * problem->M or a fresh copy), the original problem->M if it was symmetrized (NULL otherwise),
 * H_tilde and Minv. */
static void gfc3d_IPM_prepare_problem(GlobalFrictionContactProblem* problem,
                                      SolverOptions* options, size_t m, size_t nd,
                                      NumericsMatrix** M_out, NumericsMatrix** M_original_out,
                                      NumericsMatrix** H_tilde_out,
                                      NumericsMatrix** Minv_out) {
  /* symmetrization of the matrix M */
  *M_original_out = NULL;
  if (!(NM_is_symmetric(problem->M))) {
    printf("#################### SYMMETRIZATION ####################\n");
    NumericsMatrix* MT = NM_transpose(problem->M);
    NumericsMatrix* MSym = NM_add(1 / 2., problem->M, 1 / 2., MT);
    *M_original_out = problem->M; /* keep the user-provided matrix for restoration */
    problem->M = MSym;
    NM_free(MT);
  }

  /* if SICONOS_FRICTION_3D_FORCED_SPARSE_STORAGE is set, force a copy into NM_SPARSE
   * storageType */
  DEBUG_PRINTF("problem->M->storageType : %i\n", problem->M->storageType);
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_SPARSE_STORAGE] ==
          SICONOS_FRICTION_3D_FORCED_SPARSE_STORAGE &&
      problem->M->storageType == NM_SPARSE_BLOCK) {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    printf(
        "\n\n\n######################### FORCE SPARSE STORAGE "
        "#########################\n\n\n");
    *M_out = NM_copy_to_sparse(problem->M, DBL_EPSILON);
  } else {
    *M_out = problem->M;
  }

  *Minv_out = NULL;
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] ==
      SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH) {
    size_t block_number_of_M = to_size_t((*M_out)->size0) / 3;
    size_t* blocksizes_of_M = (size_t*)malloc(block_number_of_M * sizeof(size_t));
    for (size_t i = 0; i < block_number_of_M; i++) *(blocksizes_of_M + i) = 3;
    *Minv_out = NM_inverse_diagonal_block_matrix(*M_out, block_number_of_M, blocksizes_of_M);
    free(blocksizes_of_M);
  }

  DEBUG_PRINTF("problem->M->storageType : %i\n", problem->H->storageType);
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_SPARSE_STORAGE] ==
          SICONOS_FRICTION_3D_FORCED_SPARSE_STORAGE &&
      problem->H->storageType == NM_SPARSE_BLOCK) {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    *H_tilde_out = NM_copy_to_sparse(NM_transpose(problem->H), DBL_EPSILON);
  } else {
    *H_tilde_out = NM_transpose(problem->H);
  }

  numerics_printf_verbose(1, "nnz M = %8zu density = %9.4f, \t size = (%i, %i)",
                          NM_nnz(problem->M),
                          NM_nnz(problem->M) / 1.0 / problem->M->size0 / problem->M->size1,
                          problem->M->size0, problem->M->size1);
  numerics_printf_verbose(1, "nnz H = %8zu density = %9.4f, \t mu = %.2e", NM_nnz(problem->H),
                          NM_nnz(problem->H) / 1.0 / nd / m, problem->mu[0]);
}

/* Apply the friction-coefficient scaling (P_mu) to H_tilde and w, then build the derived
 * matrices needed by the IPM iterations (-H, -H^T, H^T, -M). */
static void gfc3d_IPM_scale_problem(NumericsMatrix* M, NumericsMatrix* H_tilde,
                                    GlobalFrictionContactProblem* problem,
                                    Gfc3d_IPM_init_data* data, size_t nd,
                                    NumericsMatrix** H_out, double** w_out,
                                    NumericsMatrix** minus_H_out,
                                    NumericsMatrix** minus_Ht_out, NumericsMatrix** Ht_out,
                                    NumericsMatrix** minus_M_out) {
  NumericsMatrix* P_mu = data->P_mu->mat;

  double* w_tilde = problem->b;
  double* w = data->tmp_vault_nd[0];

  NumericsMatrix* H = NM_multiply(P_mu, H_tilde);
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  NumericsMatrix* minus_H = NM_copy(H);
  NM_scal(-1.0, minus_H);
  NumericsMatrix* minus_Ht = NM_transpose(minus_H);
  NumericsMatrix* Ht = NM_transpose(H);

  NumericsMatrix* minus_M = NM_copy(M);
  NM_scal(-1.0, minus_M);

  *H_out = H;
  *w_out = w;
  *minus_H_out = minus_H;
  *minus_Ht_out = minus_Ht;
  *Ht_out = Ht;
  *minus_M_out = minus_M;
}

/* Initialize the IPM iterate (reaction, velocity, globalVelocity) with a default interior
 * point. */
static void gfc3d_IPM_compute_starting_point(size_t m, size_t nd, size_t d, NumericsMatrix* M,
                                             NumericsMatrix* H, const double* f,
                                             double* reaction, double* velocity,
                                             double* globalVelocity) {
  for (size_t i = 0; i < nd; i++)
    if (i % d == 0)
      reaction[i] = 0.1;
    else
      reaction[i] = 0.01;

  for (size_t i = 0; i < m; i++) globalVelocity[i] = f[i];
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  NM_Cholesky_solve(NM_preserve(M), globalVelocity, 1);

  for (size_t i = 0; i < nd; i++)
    if (i % d == 0)
      velocity[i] = 0.1;
    else
      velocity[i] = 0.01;
}

/* Bundle of resources that must be released at the end of gfc3d_IPM and gfc3d_IPM_fixed.
 * Keeping them in one place guarantees a single cleanup path for both functions. */
typedef struct {
  GlobalFrictionContactProblem* problem;
  SolverOptions* options;
  int internal_allocation;
  int hasNotConverged;
  int restore_scaling;
  int write_matlab_footer;
  size_t m;
  size_t nd;
  size_t iteration;
  double* velocity;
  double* reaction;
  NumericsMatrix* M;
  NumericsMatrix* M_original;
  double* globalVelocity;
  double* s;
  double* diff_fixp_vec;
  double totalresidual;
  NumericsMatrix* P_mu;
  NumericsMatrix* P_mu_inv;
  NumericsMatrix* H_tilde;
  NumericsMatrix* H;
  NumericsMatrix* minus_H;
  NumericsMatrix* Minv;
  NumericsMatrix* minus_M;
  NumericsMatrix* minus_Ht;
  NumericsMatrix* eye_nd;
  NumericsMatrix* Ht;
  FILE* iterates;

  /* extra workspace used by gfc3d_IPM */
  double* r_rhs;
  double* r_rhs_2;
  double* r_dv;
  double* r_du;
  double* r_dr;
  double* r_adu;
  double* r_adr;
  double* r_dr_a;
  double* r_du_a;
  double* r_dv_a;
  double* Hvw;
  double* rhs_2;
  double* sr_rhs;
  double* sr_rhs_2;
  double* phiu;
  double* d_globalVelocity;
  double* d_velocity;
  double* d_reaction;
  short* a_velo;
  short* a_reac;
} gfc3d_IPM_cleanup_data;

/* Single cleanup routine for gfc3d_IPM and gfc3d_IPM_fixed. */
static void gfc3d_IPM_cleanup(gfc3d_IPM_cleanup_data* cd, int* info) {
  Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)cd->options->solverData;

  if (cd->restore_scaling) {
    NM_gemv(1.0, cd->P_mu_inv, cd->velocity, 0.0, data->tmp_point->t_velocity);
    cblas_dcopy(cd->nd, data->tmp_point->t_velocity, 1, cd->velocity, 1);
    NM_gemv(1.0, cd->P_mu, cd->reaction, 0.0, data->tmp_point->t_reaction);
    cblas_dcopy(cd->nd, data->tmp_point->t_reaction, 1, cd->reaction, 1);
  }

  cd->options->dparam[SICONOS_DPARAM_RESIDU] = cd->totalresidual;
  cd->options->iparam[SICONOS_IPARAM_ITER_DONE] = (int)cd->iteration;

  if (cd->globalVelocity) {
    cblas_dcopy(cd->nd, cd->velocity, 1, data->starting_point->velocity, 1);
    cblas_dcopy(cd->nd, cd->reaction, 1, data->starting_point->reaction, 1);
    cblas_dcopy(cd->m, cd->globalVelocity, 1, data->starting_point->globalVelocity, 1);
  }

  if (cd->s) free(cd->s);
  if (cd->internal_allocation) {
    gfc3d_IPM_free(cd->problem, cd->options);
  }

  cd->options->solverData = (double*)malloc(sizeof(double));

  /* Free any working copy of M that is not the problem's own matrix (e.g. a
   * sparse-storage copy), then restore the original user-provided M if
   * symmetrization replaced it. */
  if (cd->M && cd->M != cd->problem->M) {
    cd->M = NM_free(cd->M);
  }
  if (cd->M_original) {
    cd->problem->M = NM_free(cd->problem->M);
    cd->problem->M = cd->M_original;
  }

  if (cd->H_tilde) cd->H_tilde = NM_free(cd->H_tilde);
  if (cd->minus_H) cd->minus_H = NM_free(cd->minus_H);
  if (cd->H) cd->H = NM_free(cd->H);
  if (cd->Minv) cd->Minv = NM_free(cd->Minv);
  if (cd->minus_M) cd->minus_M = NM_free(cd->minus_M);
  if (cd->minus_Ht) cd->minus_Ht = NM_free(cd->minus_Ht);
  if (cd->eye_nd) cd->eye_nd = NM_free(cd->eye_nd);
  if (cd->Ht) cd->Ht = NM_free(cd->Ht);

  if (cd->options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
    if (cd->write_matlab_footer) fprintf(cd->iterates, "];\n\n");
    fclose(cd->iterates);
  }

  if (cd->r_rhs) free(cd->r_rhs);
  if (cd->r_rhs_2) free(cd->r_rhs_2);
  if (cd->r_dv) free(cd->r_dv);
  if (cd->r_du) free(cd->r_du);
  if (cd->r_dr) free(cd->r_dr);
  if (cd->r_adu) free(cd->r_adu);
  if (cd->r_adr) free(cd->r_adr);
  if (cd->r_dr_a) free(cd->r_dr_a);
  if (cd->r_du_a) free(cd->r_du_a);
  if (cd->r_dv_a) free(cd->r_dv_a);
  if (cd->Hvw) free(cd->Hvw);
  if (cd->a_velo) free(cd->a_velo);
  if (cd->a_reac) free(cd->a_reac);
  if (cd->rhs_2) free(cd->rhs_2);
  if (cd->sr_rhs) free(cd->sr_rhs);
  if (cd->sr_rhs_2) free(cd->sr_rhs_2);

  if (cd->phiu) free(cd->phiu);

  if (cd->d_globalVelocity) free(cd->d_globalVelocity);
  if (cd->d_velocity) free(cd->d_velocity);
  if (cd->d_reaction) free(cd->d_reaction);
  if (cd->diff_fixp_vec) free(cd->diff_fixp_vec);

  *info = cd->hasNotConverged;
}

/* Workspace used by gfc3d_IPM and gfc3d_IPM_fixed.  Grouping the state in one struct lets us
 * split the per-linear-system-formulation switch bodies into separate functions without
 * passing dozens of arguments. */
typedef struct {
  /* problem dimensions */
  size_t m, nd, d, n;
  int type;
  double tol;
  size_t max_iter;

  /* handles */
  GlobalFrictionContactProblem* problem;
  SolverOptions* options;
  Gfc3d_IPM_init_data* data;

  /* algorithm parameters */
  double sgmp1, sgmp2, sgmp3;
  double gmmp1, gmmp2;
  double gmm;

  /* matrices */
  NumericsMatrix *M, *H_tilde, *H, *minus_H, *minus_Ht, *Ht, *minus_M, *Minv, *eye_nd;
  NumericsMatrix *J, *JR, *Qp2;

  /* vectors */
  double *w, *w_tilde, *f;
  double *reaction, *velocity, *globalVelocity;
  double *primalConstraint, *dualConstraint, *complemConstraint;
  double *rhs, *rhs_2;
  double *v_plus_dv, *r_plus_dr, *dvdr_jprod;
  double *r_rhs, *r_rhs_2, *sr_rhs, *sr_rhs_2;
  double *d_globalVelocity, *d_velocity, *d_reaction;
  double* Hvw;
  double* p2;
  double *velocity_t, *d_velocity_t, *d_reaction_t, *velocity_t_inv;
  double *d_s, *s;
  double* phiu;
  double *r_dv, *r_du, *r_dr, *r_dv_a, *r_du_a, *r_dr_a, *r_adu, *r_adr;
  short *a_velo, *a_reac;

  /* gfc3d_IPM_fixed specific */
  double* diff_fixp_vec;
  double *r_Qp_u, *r_Qp_u_inv;

  /* iteration state */
  int internal_allocation;
  int hasNotConverged;
  size_t iteration;
  int jacobian_is_nan;
  double pinfeas, dinfeas, complem, udotr, dualgap, projerr;
  double totalresidual, totalresidual_Jor, totalresidual_mu, diff_fixp;
  double barr_param, barr_param_a, sigma, e;
  double alpha_primal, alpha_dual, alpha_tmp;
  double LS_norm_p, LS_norm_d, LS_norm_c;
  double norm_f, norm_w;
  long J_nzmax, JR_nzmax;
  size_t H_nzmax;
  int nB, nN, nR, nT;
  double nub, unitur;

  /* misc */
  FILE* iterates;
  char fws;
  char* str;
  char* strToken;
} gfc3d_IPM_workspace;

/* ---- per-linear-system-formulation step helpers ---- */

static void gfc3d_IPM_predictor_3X3_NOSCAL(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  int iteration = (int)ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // t1 = clock();
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2 * nd, m + 2 * nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + (3 * d - 2) * n + 3 * d * n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->sparse->origin = NSM_TRIPLET;

      NumericsMatrix* arrow_r = Arrow_repr(reaction, nd, n);
      NumericsMatrix* arrow_u = Arrow_repr(velocity, nd, n);

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);
      NM_insert(J, arrow_r, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u, m, m + nd);

      /* regularization */
      // NM_insert(J, NM_scalar(nd, -fmin(1e-7,barr_param)), m + nd, m + nd);

      NM_free(arrow_r);
      NM_free(arrow_u);

      cblas_dcopy(m + 2 * nd, rhs, 1, rhs_2, 1);
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(m + 2 * nd, sizeof(double));
        cblas_dcopy(m + 2 * nd, rhs, 1, rhs_save, 1);
        // NM_LU_refine(J, rhs, rhs_save, 1, tol, 10, 0);
        double residu;
        NM_LU_refine(J, rhs, 1e-14, 10, &residu);
        free(rhs_save);
      } else
        NM_LU_solve(J, rhs, 1);

      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      JA_prod(velocity, reaction, nd, n, complemConstraint);

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, complemConstraint, 1, rhs + m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs + m + nd, 1);

      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      cblas_dcopy(m + 2 * nd, rhs, 1, rhs_2, 1);

      NM_LU_solve(J, rhs, 1);

      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs + m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs + m + nd, 1, d_reaction, 1);

      // t2 = clock();
      // total_time += (double)(t2 - t1) / (double)clk_tck;
      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_predictor_3X3_QP2(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2 * nd, m + 2 * nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->sparse->origin = NSM_TRIPLET;

      Nesterov_Todd_vector(2, velocity, reaction, nd, n, p2);
      Qp2 = QRmat(p2, nd, n);

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);
      NM_insert(J, Qp2, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, eye_nd, m, m + nd);

      if (Qp2) Qp2 = NM_free(Qp2);

      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, reaction, 1, rhs + m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs + m + nd, 1);

      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      cblas_dcopy(m + 2 * nd, rhs, 1, rhs_2, 1);

      NM_LDLT_solve(J, rhs, 1);

      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs + m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs + m + nd, 1, d_reaction, 1);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_predictor_3X3_QPH(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2 * nd, m + 2 * nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->sparse->origin = NSM_TRIPLET;

      NumericsMatrix* minusQpH = QNTpH(velocity, reaction, H, nd, n);
      NM_scal(-1.0, minusQpH);
      NumericsMatrix* minusQpHt = NM_transpose(minusQpH);

      NM_insert(J, M, 0, 0);
      NM_insert(J, minusQpH, m + nd, 0);
      NM_insert(J, eye_nd, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, minusQpHt, 0, m + nd);
      NM_insert(J, eye_nd, m, m + nd);

      NM_free(minusQpH);
      NM_free(minusQpHt);

      double* r_check = (double*)calloc(nd, sizeof(double));
      QNTpinvz(velocity, reaction, reaction, nd, n, r_check);

      double* Qp_primalConstraint = (double*)calloc(nd, sizeof(double));
      QNTpz(velocity, reaction, primalConstraint, nd, n, Qp_primalConstraint);
      cblas_dcopy(nd, Qp_primalConstraint, 1, primalConstraint, 1);

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, r_check, 1, rhs + m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs + m + nd, 1);
      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      cblas_dcopy(m + 2 * nd, rhs, 1, rhs_2, 1);

      free(r_check);
      free(Qp_primalConstraint);
      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(m + 2 * nd, sizeof(double));
        cblas_dcopy(m + 2 * nd, rhs, 1, rhs_save, 1);
        NM_LDLT_refine(J, rhs, rhs_save, 1, 1e-14, 10, 0);
        free(rhs_save);
      } else {
        NM_LDLT_solve(J, rhs, 1);
      }

      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      QNTpinvz(velocity, reaction, rhs + m, nd, n, d_velocity);
      QNTpz(velocity, reaction, rhs + m + nd, nd, n, d_reaction);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_predictor_2X2_QP2(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // First linear linear system
      JR = NM_create(NM_SPARSE, m + nd, m + nd);
      JR_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
      NM_triplet_alloc(JR, JR_nzmax);
      JR->sparse->origin = NSM_TRIPLET;

      Nesterov_Todd_vector(3, velocity, reaction, nd, n, p2);
      Qp2 = QRmat(p2, nd, n);

      NM_insert(JR, minus_M, 0, 0);
      NM_insert(JR, Ht, 0, m);
      NM_insert(JR, H, m, 0);
      NM_insert(JR, Qp2, m, m);

      if (Qp2) Qp2 = NM_free(Qp2);

      jacobian_is_nan = NM_isnan(JR);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      NV_insert(r_rhs, m + nd, dualConstraint, m, 0);
      NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);  // Hvw <- H*v
      cblas_daxpy(nd, 1.0, w, 1, Hvw,
                  1);  // Hvw <- H*v + w. Hvw will be used for the computation of r_du
      NV_insert(r_rhs, m + nd, Hvw, nd, m);
      cblas_dscal(nd, -1.0, r_rhs + m, 1);

      cblas_dcopy(m + nd, r_rhs, 1, r_rhs_2, 1);

      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, r_rhs, 1);

      cblas_dcopy(m, r_rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, r_rhs + m, 1, d_reaction, 1);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity,
                  1);  // d_velocity <- H*(v + dv) + w  - u

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_predictor_2X2_QPH(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // t1 = clock();
      // First linear linear system
      JR = NM_create(NM_SPARSE, m + nd, m + nd);
      JR_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
      NM_triplet_alloc(JR, JR_nzmax);
      JR->sparse->origin = NSM_TRIPLET;

      NumericsMatrix* QpH = QNTpH(velocity, reaction, H, nd, n);
      NumericsMatrix* QpHt = NM_transpose(QpH);

      NM_insert(JR, minus_M, 0, 0);
      NM_insert(JR, QpH, m, 0);
      NM_insert(JR, QpHt, 0, m);
      NM_insert(JR, eye_nd, m, m);

      if (QpH) QpH = NM_free(QpH);
      if (QpHt) QpHt = NM_free(QpHt);

      jacobian_is_nan = NM_isnan(JR);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      NV_insert(r_rhs, m + nd, dualConstraint, m, 0);
      NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);
      cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);  // Hvw will be used for the computation of r_du
      QNTpz(velocity, reaction, Hvw, nd, n, r_rhs + m);
      cblas_dscal(nd, -1.0, r_rhs + m, 1);

      cblas_dcopy(m + nd, r_rhs, 1, r_rhs_2, 1);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(m + nd, sizeof(double));
        cblas_dcopy(m + nd, r_rhs, 1, rhs_save, 1);
        NM_LDLT_refine(JR, r_rhs, rhs_save, 1, 1e-14, 10, 0);
        free(rhs_save);
      } else {
        NM_LDLT_solve(JR, r_rhs, 1);
      }

      cblas_dcopy(m, r_rhs, 1, d_globalVelocity, 1);
      QNTpz(velocity, reaction, r_rhs + m, nd, n, d_reaction);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity,
                  1);  // d_velocity <- H*(v + dv) + w  - u

      // t2 = clock();
      // total_time += (double)(t2 - t1) / (double)clk_tck;

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_predictor_1X1_QPH(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // First linear linear system
      JR = NM_create(NM_SPARSE, nd, nd);
      // JR_nzmax = 2*H_nzmax + nd;
      // NM_triplet_alloc(JR, JR_nzmax);
      // JR->sparse->origin = NSM_TRIPLET;

      NumericsMatrix* JR_a = NULL;
      NumericsMatrix* JR_b = NULL;

      NumericsMatrix* QpH = QNTpH(velocity, reaction, H, nd, n);
      NumericsMatrix* QpHt = NM_transpose(QpH);

      JR_a = NM_multiply(Minv, QpHt);
      JR_b = NM_multiply(QpH, JR_a);
      JR = NM_add(1.0, JR_b, 1.0, eye_nd);
      JR_a = NM_free(JR_a);
      JR_b = NM_free(JR_b);

      jacobian_is_nan = NM_isnan(JR);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      //      if (iteration == 0) printf("NNZ1X1_QPH %zu\n",NM_nnz(JR));

      double* fHr = (double*)calloc(m, sizeof(double));
      double* MfHr = (double*)calloc(m, sizeof(double));
      double* vdv = (double*)calloc(m, sizeof(double));
      double* rdr = (double*)calloc(nd, sizeof(double));
      double* w_h = (double*)calloc(nd, sizeof(double));

      cblas_dcopy(m, f, 1, fHr, 1);
      NM_gemv(1.0, Ht, reaction, 1.0, fHr);  // fHr <- H'*r + f
      NM_gemv(1.0, Minv, fHr, 0.0, MfHr);    // MfHr <- M^{-1}* (H'*r + f)
      NM_gemv(-1.0, QpH, MfHr, 0.0, sr_rhs);
      QNTpz(velocity, reaction, w, nd, n, w_h);
      cblas_daxpy(nd, -1.0, w_h, 1, sr_rhs, 1);

      cblas_dcopy(nd, sr_rhs, 1, sr_rhs_2, 1);

      /* NM_Cholesky_solve(JR, sr_rhs, 1); */

      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, sr_rhs, 1);

      QNTpz(velocity, reaction, sr_rhs, nd, n, d_reaction);

      NV_add(reaction, d_reaction, nd, rdr);
      cblas_dcopy(m, f, 1, MfHr, 1);
      NM_gemv(1.0, Ht, rdr, 1.0, MfHr);
      NM_gemv(1.0, Minv, MfHr, 0.0, d_globalVelocity);
      cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity, 1);

      NV_add(globalVelocity, d_globalVelocity, m, vdv);
      cblas_dcopy(nd, w, 1, d_velocity, 1);
      NM_gemv(1.0, H, vdv, 1.0, d_velocity);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      free(fHr);
      free(MfHr);
      free(vdv);
      free(rdr);
      free(w_h);
      QpH = NM_free(QpH);
      QpHt = NM_free(QpHt);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_corrector_3X3_NOSCAL(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // t1 = clock();

      JA_prod(d_velocity, d_reaction, nd, n, dvdr_jprod);
      cblas_daxpy(nd, -1.0, dvdr_jprod, 1, rhs_2 + m, 1);
      for (size_t k = 0; k < nd; rhs_2[m + k] += 2 * sigma * barr_param, k += d);

      double* rhs_tmp = (double*)calloc(m + 2 * nd, sizeof(double));
      cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_tmp, 1);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(m + 2 * nd, sizeof(double));
        cblas_dcopy(m + 2 * nd, rhs_tmp, 1, rhs_save, 1);
        // NM_LU_refine(J, rhs, rhs_save, 1, tol, 10, 0);
        double residu;
        NM_LU_refine(J, rhs_tmp, 1e-14, 10, &residu);
        free(rhs_save);
      } else
        NM_LU_solve(J, rhs_tmp, 1);

      cblas_dcopy(m, rhs_tmp, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs_tmp + m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs_tmp + m + nd, 1, d_reaction, 1);
      NM_gemv(1.0, J, rhs_tmp, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2 + m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2 + m + nd, 1);

      // t2 = clock();
      // total_time += (double)(t2 - t1) / (double)clk_tck;

      free(rhs_tmp);
      // free(sol);

      J = NM_free(J);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_corrector_3X3_QP2(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // Second linear linear system
      QNTpz(velocity, reaction, d_velocity, nd, n, d_velocity_t);
      QNTpinvz(velocity, reaction, d_reaction, nd, n, d_reaction_t);
      JA_prod(d_velocity_t, d_reaction_t, nd, n, dvdr_jprod);
      QNTpz(velocity, reaction, velocity, nd, n, velocity_t);
      JA_inv(velocity_t, nd, n, velocity_t_inv);
      JA_prod(velocity_t_inv, dvdr_jprod, nd, n, d_velocity_t);
      QNTpz(velocity, reaction, d_velocity_t, nd, n, dvdr_jprod);
      cblas_daxpy(nd, -1.0, dvdr_jprod, 1, rhs_2 + m, 1);

      JA_inv(velocity, nd, n, d_velocity_t);
      for (size_t k = 0; k < nd;
           rhs_2[m + k] += 2 * sigma * barr_param * d_velocity_t[k], k++);

      double* rhs_tmp = (double*)calloc(m + 2 * nd, sizeof(double));
      cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_tmp, 1);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(m + 2 * nd, sizeof(double));
        cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_save, 1);
        NM_LDLT_refine(J, rhs_2, rhs_save, 1, 1e-14, 10, 0);
        free(rhs_save);
      } else {
        NM_LDLT_solve(J, rhs_2, 1);
      }

      cblas_dcopy(m, rhs_2, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs_2 + m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs_2 + m + nd, 1, d_reaction, 1);

      NM_gemv(1.0, J, rhs_2, -1.0, rhs_tmp);

      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp + m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_tmp + m + nd, 1);

      free(rhs_tmp);

      J = NM_free(J);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_corrector_3X3_QPH(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // Second linear linear system
      JA_prod(rhs + m, rhs + m + nd, nd, n, dvdr_jprod);
      for (size_t k = 0; k < nd; dvdr_jprod[k] -= 2 * sigma * barr_param, k += d);

      QNTpz(velocity, reaction, velocity, nd, n, velocity_t);
      JA_inv(velocity_t, nd, n, velocity_t_inv);
      JA_prod(velocity_t_inv, dvdr_jprod, nd, n, d_velocity_t);
      cblas_daxpy(nd, -1.0, d_velocity_t, 1, rhs_2 + m, 1);

      double* rhs_tmp = (double*)calloc(m + 2 * nd, sizeof(double));
      cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_tmp, 1);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(m + 2 * nd, sizeof(double));
        cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_save, 1);
        NM_LDLT_refine(J, rhs_2, rhs_save, 1, tol, 10, 0);
        free(rhs_save);
      } else {
        NM_LDLT_solve(J, rhs_2, 1);
      }

      NM_gemv(1.0, J, rhs_2, -1.0, rhs_tmp);

      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp + m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_tmp + m + nd, 1);

      free(rhs_tmp);

      cblas_dcopy(m, rhs_2, 1, d_globalVelocity, 1);
      QNTpinvz(velocity, reaction, rhs_2 + m, nd, n, d_velocity);
      QNTpz(velocity, reaction, rhs_2 + m + nd, nd, n, d_reaction);

      if (J) J = NM_free(J);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_corrector_2X2_QP2(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // Second linear linear system
      double* r_Qp_u = (double*)calloc(nd, sizeof(double));
      double* r_Qp_du = (double*)calloc(nd, sizeof(double));
      double* r_dudr = (double*)calloc(nd, sizeof(double));
      double* r_ududr = (double*)calloc(nd, sizeof(double));
      double* r_Qpinv_dr = (double*)calloc(nd, sizeof(double));

      QNTpz(velocity, reaction, d_velocity, nd, n, r_Qp_du);        // du_affine_hat
      QNTpinvz(velocity, reaction, d_reaction, nd, n, r_Qpinv_dr);  // dr_affine_check
      JA_prod(r_Qp_du, r_Qpinv_dr, nd, n,
              r_dudr);  // Jordan product du_affine__hat o dr_affine_check
      for (size_t k = 0; k < nd; r_dudr[k] -= 2 * sigma * barr_param,
                  k += d);  // du_affine_hat o dr_affine_check - 2*sigma*mu*e
      QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);  // u_hat
      Jxinvprody(r_Qp_u, r_dudr, nd, n,
                 r_ududr);  // Jordan product u_hat_inv o (du_affine__hat o dr_affine_check)
      QNTpinvz(velocity, reaction, r_ududr, nd, n,
               r_ududr);  // Qp_inv * u_hat_inv o du_affine_hat o dr_affine_check -
                          // 2*sigma*mu*u_hatinv)
      cblas_daxpy(nd, -1.0, r_ududr, 1, r_rhs_2 + m, 1);

      double* rhs_tmp = (double*)calloc(m + nd, sizeof(double));
      cblas_dcopy(m + nd, r_rhs_2, 1, rhs_tmp, 1);

      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, r_rhs_2, 1);

      cblas_dcopy(m, r_rhs_2, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, r_rhs_2 + m, 1, d_reaction, 1);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      NM_gemv(1.0, JR, r_rhs_2, -1.0, rhs_tmp);
      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp + m, 1);

      free(rhs_tmp);

      if (JR) JR = NM_free(JR);
      free(r_Qp_u);
      free(r_Qp_du);
      free(r_dudr);
      free(r_ududr);
      free(r_Qpinv_dr);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_corrector_2X2_QPH(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // t1 = clock();

      // Second linear linear system
      double* r_Qp_u = (double*)calloc(nd, sizeof(double));
      double* r_Qp_du = (double*)calloc(nd, sizeof(double));
      double* r_dudr = (double*)calloc(nd, sizeof(double));
      double* r_ududr = (double*)calloc(nd, sizeof(double));

      QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);  // r_Qp_u <- u_hat  (Qp formula)
      cblas_dcopy(nd, r_Qp_u, 1, r_Qp_du, 1);              // r_Qp_du <- u_hat
      cblas_daxpy(nd, 1.0, r_rhs + m, 1, r_Qp_du, 1);      // r_Qp_du <- dr_check + u_hat
      cblas_dscal(nd, -1.0, r_Qp_du,
                  1);  // r_Qp_du <- -(dr_check + u_hat) = du_hat.
                       // Here we use the following formula satisfied by the affine step
                       //    u_hat o dr_check + r_check o du_hat = - u_hat o r_check
                       // Then we use the fact that for the NT scaling we have: u_hat =
                       // r_check. Therefore du_hat = -dr_check - r_check = -dr_check - u_hat

      JA_prod(r_Qp_du, r_rhs + m, nd, n, r_dudr);  // r_dudr <- du_hat o dr_check

      for (size_t k = 0; k < nd; r_dudr[k] -= 2 * sigma * barr_param,
                  k += d);  // r_dudr <- du_hat o dr_check - 2*sigma*mu*e
      Jxinvprody(r_Qp_u, r_dudr, nd, n,
                 r_ududr);  // r_ududr <- u_hat^{-1} o  (du_hat o dr_check - 2*sigma*mu*e)

      cblas_daxpy(nd, -1.0, r_ududr, 1, r_rhs_2 + m,
                  1);  // r_rhs_2+m <- -Qp*(Hv+w) - u_hat^{-1} o (du_hat o dr_check) +
                       // 2*sigma*mu*u_hat^{-1}

      double* rhs_tmp = (double*)calloc(m + nd, sizeof(double));
      cblas_dcopy(m + nd, r_rhs_2, 1, rhs_tmp, 1);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(m + nd, sizeof(double));
        cblas_dcopy(m + nd, r_rhs_2, 1, rhs_save, 1);
        NM_LDLT_refine(JR, r_rhs_2, rhs_save, 1, 1e-14, 10, 0);
        free(rhs_save);
      } else {
        NM_LDLT_solve(JR, r_rhs_2, 1);
      }

      cblas_dcopy(m, r_rhs_2, 1, d_globalVelocity, 1);
      QNTpz(velocity, reaction, r_rhs_2 + m, nd, n, d_reaction);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      NM_gemv(1.0, JR, r_rhs_2, -1.0, rhs_tmp);
      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp + m, 1);
      free(rhs_tmp);

      // t2 = clock();
      // total_time += (double)(t2 - t1) / (double)clk_tck;

      if (JR) JR = NM_free(JR);
      free(r_Qp_u);
      free(r_Qp_du);
      free(r_dudr);
      free(r_ududr);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

static void gfc3d_IPM_corrector_1X1_QPH(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // Second linear linear system
      double* vdv = (double*)calloc(m, sizeof(double));
      double* rdr = (double*)calloc(nd, sizeof(double));
      double* rhs_cor = (double*)calloc(nd, sizeof(double));
      double* MfHr = (double*)calloc(m, sizeof(double));

      QNTpz(velocity, reaction, d_velocity, nd, n, d_velocity_t);
      JA_prod(d_velocity_t, sr_rhs, nd, n, dvdr_jprod);
      for (size_t k = 0; k < nd; dvdr_jprod[k] -= 2 * sigma * barr_param, k += d);
      QNTpz(velocity, reaction, velocity, nd, n, velocity_t);
      JA_inv(velocity_t, nd, n, velocity_t_inv);
      JA_prod(velocity_t_inv, dvdr_jprod, nd, n, rhs_cor);
      cblas_daxpy(nd, -1.0, rhs_cor, 1, sr_rhs_2, 1);

      double* rhs_tmp = (double*)calloc(nd, sizeof(double));
      cblas_dcopy(nd, sr_rhs_2, 1, rhs_tmp, 1);

      /* NM_Cholesky_solve(JR, sr_rhs_2, 1); */

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        double* rhs_save = (double*)calloc(nd, sizeof(double));
        cblas_dcopy(nd, sr_rhs_2, 1, rhs_save, 1);
        NM_LDLT_refine(JR, sr_rhs_2, rhs_save, 1, 1e-10, 10, 0);
        free(rhs_save);
      } else {
        NM_LDLT_solve(JR, sr_rhs_2, 1);
      }

      QNTpz(velocity, reaction, sr_rhs_2, nd, n, d_reaction);

      NV_add(reaction, d_reaction, nd, rdr);
      cblas_dcopy(m, f, 1, MfHr, 1);
      NM_gemv(1.0, Ht, rdr, 1.0, MfHr);
      NM_gemv(1.0, Minv, MfHr, 0.0, d_globalVelocity);
      cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity, 1);

      NV_add(globalVelocity, d_globalVelocity, m, vdv);
      cblas_dcopy(nd, w, 1, d_velocity, 1);
      NM_gemv(1.0, H, vdv, 1.0, d_velocity);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      NM_gemv(1.0, JR, sr_rhs_2, -1.0, rhs_tmp);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp, 1);
      free(rhs_tmp);

      JR = NM_free(JR);
      free(vdv);
      free(rdr);
      free(rhs_cor);
      free(MfHr);

      break;
    }
  }

  /* ---- copy back modified state ---- */
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->hasNotConverged = hasNotConverged;
  ws->J = J;
  ws->JR = JR;
  ws->Qp2 = Qp2;
  ws->J_nzmax = J_nzmax;
  ws->JR_nzmax = JR_nzmax;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->dualgap = dualgap;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->barr_param = barr_param;
  ws->barr_param_a = barr_param_a;
  ws->sigma = sigma;
  ws->e = e;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->alpha_tmp = alpha_tmp;
  ws->gmm = gmm;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->nB = nB;
  ws->nN = nN;
  ws->nR = nR;
  ws->nT = nT;
  ws->nub = nub;
  ws->unitur = unitur;
}

void gfc3d_IPM(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, double* restrict globalVelocity, int* restrict info,
               SolverOptions* restrict options) {
  // verbose = 3;

  // int type = NORM_2;
  int type = NORM_INF;

  // the size of the problem detection
  size_t m = problem->M->size0;
  size_t nd = problem->H->size1;
  size_t d = problem->dimension;
  size_t n = problem->numberOfContacts;

  NumericsMatrix* M = NULL;
  NumericsMatrix* M_original = NULL;
  NumericsMatrix* H_tilde = NULL;

  /* globalFrictionContact_display(problem); */

  NumericsMatrix* Minv = NULL;
  gfc3d_IPM_prepare_problem(problem, options, m, nd, &M, &M_original, &H_tilde, &Minv);

  // initialize solver if it is not set
  int internal_allocation = 0;
  size_t work_size = m + nd + nd;
  if (!options->dWork || (options->dWorkSize < work_size)) {
    gfc3d_IPM_init(problem, options, work_size);
    internal_allocation = 1;
  }

  Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)options->solverData;
  NumericsMatrix* P_mu = data->P_mu->mat;
  NumericsMatrix* P_mu_inv = data->P_mu->inv_mat;

  double* w_tilde = problem->b;
  double* f = problem->q;

  NumericsMatrix* H = NULL;
  double* w = NULL;
  NumericsMatrix* minus_H = NULL;
  NumericsMatrix* minus_Ht = NULL;
  NumericsMatrix* Ht = NULL;
  NumericsMatrix* minus_M = NULL;
  gfc3d_IPM_scale_problem(M, H_tilde, problem, data, nd, &H, &w, &minus_H, &minus_Ht, &Ht,
                          &minus_M);

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;
  double alpha_tmp = 0.0;

  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  /* COMPUTATION OF A NEW STARTING POINT */
  gfc3d_IPM_compute_starting_point(m, nd, d, M, H, f, reaction, velocity, globalVelocity);

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  size_t max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  size_t iteration = 0;
  double pinfeas = 1e300;
  //  double pinfeas_new = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  // double complem_p = 1e300;
  double dualgap = 1e300;
  double udotr = 1e300;
  double projerr = 1e300;
  double error[6];
  double totalresidual = 1e300, totalresidual_Jor = 1e300;
  // double diff_fixp = 1e300, nub = 1e300;

  double* primalConstraint = data->tmp_vault_nd[1];
  double* dualConstraint = data->tmp_vault_m[0];
  double* complemConstraint = data->tmp_vault_nd[2];

  double gmm = gmmp1 + gmmp2;
  double barr_param_a = 0.0, e = 0.0;  //, barr_param_fixed;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  double* d_globalVelocity = (double*)calloc(m, sizeof(double));
  double* d_velocity = (double*)calloc(nd, sizeof(double));
  double* d_reaction = (double*)calloc(nd, sizeof(double));
  // double* dw = (double*)calloc(nd, sizeof(double));
  //  double* u_plus_phiu = (double*)calloc(nd, sizeof(double));
  double* d_s = (double*)calloc(n, sizeof(double));
  double* s = (double*)calloc(n, sizeof(double));
  for (size_t i = 0; i < n; i++) s[i] = cblas_dnrm2(2, velocity + i * d + 1, 1);

  /* double *tmpsol = (double*)calloc(m+2*nd,sizeof(double)); // temporary solution */

  double* rhs = options->dWork;
  double* rhs_2 = (double*)calloc(m + 2 * nd, sizeof(double));
  // double* gv_plus_dgv = data->tmp_vault_m[1];
  // double* vr_jprod = data->tmp_vault_nd[3];
  double* v_plus_dv = data->tmp_vault_nd[4];
  double* r_plus_dr = data->tmp_vault_nd[5];
  //  double* vr_prod_sub_iden = data->tmp_vault_nd[6];
  double* dvdr_jprod = data->tmp_vault_nd[7];

  // double * r_p = (double*)calloc(nd,sizeof(double));                          // scaling
  // vector p
  // NumericsMatrix* r_Qp = NULL;  // matrix Qp
  // NumericsMatrix *QpH = NM_create(H->storageType, H->size0, H->size1);        // store the
  // matrix Qp*H
  double* r_rhs = (double*)calloc(m + nd, sizeof(double));
  double* r_rhs_2 = (double*)calloc(m + nd, sizeof(double));
  double* sr_rhs = (double*)calloc(nd, sizeof(double));
  double* sr_rhs_2 = (double*)calloc(nd, sizeof(double));
  double* r_dv = (double*)calloc(m, sizeof(double));
  double* r_dr = (double*)calloc(nd, sizeof(double));
  double* r_du = (double*)calloc(nd, sizeof(double));
  double* r_dv_a = (double*)calloc(m, sizeof(double));
  double* r_dr_a = (double*)calloc(nd, sizeof(double));
  double* r_du_a = (double*)calloc(nd, sizeof(double));
  double* r_adu = (double*)calloc(nd, sizeof(double));
  double* r_adr = (double*)calloc(nd, sizeof(double));
  // double r_alpha_p, r_alpha_d; /* primal and dual steplengths */
  // double r_alpha_primal, r_alpha_dual;
  // double r_mu, r_mu_a;    /* duality gap, affine duality gap */
  NumericsMatrix* JR = 0; /* Reduced Jacobian with NT scaling */

  long JR_nzmax = 0;
  double* Hvw = (double*)calloc(nd, sizeof(double));
  // double err = 1e300;
  char fws = ' '; /* finish without scaling */

  /* list of active constraints : = 0 if x_0 <= epsilon, = 1 if lambda_2 <= epsilon , = 3
   * either */
  short* a_velo = (short*)calloc(n, sizeof(short));
  short* a_reac = (short*)calloc(n, sizeof(short));

  /* norm of the residuals of teh second linear system */
  double LS_norm_p = 0;  // primal feasibility
  double LS_norm_d = 0;  // dual feaqsibility
  double LS_norm_c = 0;  // complementarity

  // Timer
  // long clk_tck = CLOCKS_PER_SEC;
  // clock_t t1, t2;
  //  double total_time = 0;  //, max_time = 5.;

  int nB, nN, nR, nT;
  nB = nN = nR = nT = 0;

  NumericsMatrix* J = 0;
  long J_nzmax = 0;

  size_t H_nzmax = NM_nnz(H);

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] ==
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES) {
    numerics_printf_verbose(1, "---- GFC3D - IPM - Problem information");
    numerics_printf_verbose(1, "---- GFC3D - IPM - 1-norm of M = %g norm of f = %g ",
                            NM_norm_1(M), norm_f);
    numerics_printf_verbose(1, "---- GFC3D - IPM - inf-norm of M = %g ", NM_norm_inf(M));

    numerics_printf_verbose(1, "---- GFC3D - IPM - 1-norm of H = %g norm of w = %g ",
                            NM_norm_1(problem->H), norm_w);
    numerics_printf_verbose(1, "---- GFC3D - IPM - inf-norm of H = %g ",
                            NM_norm_inf(problem->H));
    numerics_printf_verbose(1, "---- GFC3D - IPM - M is symmetric = %i ", NM_is_symmetric(M));

    numerics_printf_verbose(1, "---- GFC3D - IPM - M size = (%i, %i) ", M->size0, M->size1);
    numerics_printf_verbose(1, "---- GFC3D - IPM - H size = (%i, %i) ", problem->H->size0,
                            problem->H->size1);
  }

  /* ---- IPM iterations ---- */
  numerics_printf_verbose(1, "problem dimensions n, nd x m: %1i, %6i x %-6i", n, nd, m);
  switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL: {
      numerics_printf_verbose(1, "LS solution: 3x3 no scaling\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2: {
      numerics_printf_verbose(1, "LS solution: 3x3 NT scaling with Qp2\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH: {
      numerics_printf_verbose(1, "LS solution: 3x3 NT scaling with QpH\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2: {
      numerics_printf_verbose(1, "LS solution: 2x2 NT scaling with Qp2\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH: {
      // numerics_printf_verbose(1,"LS solution: 2x2 NT scaling with QpH\n");
      if (type == NORM_2)
        numerics_printf_verbose(
            -1, "LS solution: 2x2 NT scaling with QpH,\t stopping test: norm 2\n");
      else if (type == NORM_INF)
        numerics_printf_verbose(
            -1, "LS solution: 2x2 NT scaling with QpH,\t stopping test: norm inf\n");
      else
        numerics_printf_verbose(1, "LS solution: 2x2 NT scaling with QpH\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH: {
      numerics_printf_verbose(1, "LS solution: 1x1 NT scaling with QpH\n");
      break;
    }
    default: {
      printf("ERROR\n");
    }
  }

  numerics_printf_verbose(
      1,
      "| it  |  dualgap | pinfeas | dinfeas |  <u, r> | complem | prj err | barpram | alpha_p "
      "| alpha_d |  sigma  |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");
  numerics_printf_verbose(
      1,
      "---------------------------------------------------------------------------------------"
      "--------------------------------------------------------------------------------");

  // double* p = data->tmp_vault_nd[8];
  double* p2 = data->tmp_vault_nd[9];
  // double* pinv = data->tmp_vault_nd[10];
  NumericsMatrix* Qp = NULL;
  NumericsMatrix* Qpinv = NULL;
  // NumericsMatrix* F = NULL;
  // NumericsMatrix* Finv = NULL;
  // NumericsMatrix* tmpmat = NULL;
  // NumericsMatrix* Qp_F = NULL;
  NumericsMatrix* Qp2 = NULL;
  //  NumericsMatrix* F2 = NULL;
  NumericsMatrix* eye_nd = NM_eye(nd);
  // NumericsMatrix* eye_phiu = NULL;
  double* velocity_t = data->tmp_vault_nd[11];
  double* d_velocity_t = data->tmp_vault_nd[12];
  double* d_reaction_t = data->tmp_vault_nd[13];
  double* velocity_t_inv = data->tmp_vault_nd[14];
  // double* Qp_velocity_t_inv = data->tmp_vault_nd[15];
  //  double * tmp1 = data->tmp_vault_nd[16];
  FILE* iterates = NULL;
  // FILE* matrixH;

  char* str = (char*)malloc(200);
  if (problem->name) {
    strcpy(str, problem->name);
  } else {
    strcpy(str, "foo_");
  }

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

  char matlab_name[100];
  snprintf(matlab_name, sizeof(matlab_name), "iterates.m");
  free(str);
  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
    iterates = fopen(matlab_name, "w");
    printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates);
  }

  // ComputeErrorGlobalPtr computeError = NULL;
  // computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error_convex;
  // if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
  // {
  //   computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;
  // }
  // else
  // {
  //   computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error_convex;
  // }
  /* check the full criterion */
  // double norm_q = cblas_dnrm2(m, problem->q, 1);
  // double norm_b = cblas_dnrm2(nd, problem->b, 1);

  double* phiu = (double*)calloc(n, sizeof(double));
  // double phiu_val = 0.0, alpha_diff_phi = 0.0;
  // int load_starting_point = 0, save_sol_point = 0;
  // FILE* sol_file = NULL;
  // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n ;
  /* ---- Pack all workspace into a single struct for the per-LS-form step functions ---- */
  gfc3d_IPM_workspace ws = {0};
  ws.m = m;
  ws.nd = nd;
  ws.d = d;
  ws.n = n;
  ws.type = type;
  ws.tol = tol;
  ws.max_iter = max_iter;
  ws.problem = problem;
  ws.options = options;
  ws.data = data;
  ws.sgmp1 = sgmp1;
  ws.sgmp2 = sgmp2;
  ws.sgmp3 = sgmp3;
  ws.gmmp1 = gmmp1;
  ws.gmmp2 = gmmp2;
  ws.gmm = gmm;
  ws.M = M;
  ws.H_tilde = H_tilde;
  ws.H = H;
  ws.minus_H = minus_H;
  ws.minus_Ht = minus_Ht;
  ws.Ht = Ht;
  ws.minus_M = minus_M;
  ws.Minv = Minv;
  ws.eye_nd = eye_nd;
  ws.J = J;
  ws.JR = JR;
  ws.Qp2 = Qp2;
  ws.w = w;
  ws.w_tilde = w_tilde;
  ws.f = f;
  ws.reaction = reaction;
  ws.velocity = velocity;
  ws.globalVelocity = globalVelocity;
  ws.primalConstraint = primalConstraint;
  ws.dualConstraint = dualConstraint;
  ws.complemConstraint = complemConstraint;
  ws.rhs = rhs;
  ws.rhs_2 = rhs_2;
  ws.v_plus_dv = v_plus_dv;
  ws.r_plus_dr = r_plus_dr;
  ws.dvdr_jprod = dvdr_jprod;
  ws.r_rhs = r_rhs;
  ws.r_rhs_2 = r_rhs_2;
  ws.sr_rhs = sr_rhs;
  ws.sr_rhs_2 = sr_rhs_2;
  ws.d_globalVelocity = d_globalVelocity;
  ws.d_velocity = d_velocity;
  ws.d_reaction = d_reaction;
  ws.Hvw = Hvw;
  ws.p2 = p2;
  ws.velocity_t = velocity_t;
  ws.d_velocity_t = d_velocity_t;
  ws.d_reaction_t = d_reaction_t;
  ws.velocity_t_inv = velocity_t_inv;
  ws.d_s = d_s;
  ws.s = s;
  ws.phiu = phiu;
  ws.r_dv = r_dv;
  ws.r_du = r_du;
  ws.r_dr = r_dr;
  ws.r_dv_a = r_dv_a;
  ws.r_du_a = r_du_a;
  ws.r_dr_a = r_dr_a;
  ws.r_adu = r_adu;
  ws.r_adr = r_adr;
  ws.a_velo = a_velo;
  ws.a_reac = a_reac;
  ws.internal_allocation = internal_allocation;
  ws.hasNotConverged = hasNotConverged;
  ws.iteration = iteration;
  ws.pinfeas = pinfeas;
  ws.dinfeas = dinfeas;
  ws.complem = complem;
  ws.udotr = udotr;
  ws.dualgap = dualgap;
  ws.projerr = projerr;
  ws.totalresidual = totalresidual;
  ws.totalresidual_Jor = totalresidual_Jor;
  ws.barr_param = barr_param;
  ws.barr_param_a = barr_param_a;
  ws.sigma = sigma;
  ws.e = e;
  ws.alpha_primal = alpha_primal;
  ws.alpha_dual = alpha_dual;
  ws.alpha_tmp = alpha_tmp;
  ws.LS_norm_p = LS_norm_p;
  ws.LS_norm_d = LS_norm_d;
  ws.LS_norm_c = LS_norm_c;
  ws.norm_f = norm_f;
  ws.norm_w = norm_w;
  ws.J_nzmax = J_nzmax;
  ws.JR_nzmax = JR_nzmax;
  ws.H_nzmax = H_nzmax;
  ws.nB = nB;
  ws.nN = nN;
  ws.nR = nR;
  ws.nT = nT;
  ws.iterates = iterates;
  ws.fws = fws;
  ws.str = str;
  ws.strToken = strToken;

  while (iteration < max_iter) {
    // t1 = clock();
    /** Correction of w to take into account the dependence on the tangential velocity */
    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] ==
        SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_INTERNAL) {
      for (size_t i = 0; i < nd; ++i) {
        if (i % d == 0) {
          w[i] = w_tilde[i] + cblas_dnrm2(2, velocity + i + 1, 1);
        }
      }
    } else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] ==
               SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_EXTERNAL) {
      primalResidual_type(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol,
                          type);
      dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol,
                        type);
      udotr = xdoty_type(n, nd, velocity, reaction, type);
      complem = complemResidualNorm_type(velocity, reaction, nd, n, NORM_2_INF);

      if (fmax(pinfeas, fmax(dinfeas, complem)) < tol) {
        numerics_printf_verbose(
            1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |", iteration, fws,
            dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param);
        iteration++;

        for (size_t i = 0; i < nd; ++i) {
          if (i % d == 0) {
            w[i] = w_tilde[i] + cblas_dnrm2(2, velocity + i + 1, 1);
          }
        }
        printf("\n");
        numerics_printf_verbose(1,
                                "| it  |  dualgap | pinfeas | dinfeas |  <u, r> | complem | "
                                "prj err | barpram | alpha_p | alpha_d |  sigma  |  |dv|   |  "
                                "|du|   |  |dr|   | ls prim | ls dual | ls comp |");
        numerics_printf_verbose(1,
                                "-------------------------------------------------------------"
                                "-------------------------------------------------------------"
                                "---------------------------------------------");
      }
    }

    /* Computation of the values of
     - primal residual: u - H*v - w
     - dual residual: M*v - f - H'*r
     - duality gap: u'*r
     - true duality gap: (value_of_primal - value_of_dual)
     - complementarity: u o r
     - projection error: r - proj(r-u)
    */

    // Stopping test using norm type
    primalResidual_type(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol, type);
    dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol, type);
    udotr = xdoty_type(n, nd, velocity, reaction, type);
    dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
    // complem = complemResidualNorm(velocity, reaction, nd, n);
    complem = complemResidualNorm_type(velocity, reaction, nd, n, NORM_2_INF);

    // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n / 3.;
    // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0) {
      barr_param = complemResidualNorm(velocity, reaction, nd, n) / n;
      // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;
    } else {
      barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n / 3.;
    }

    /* Keep the workspace copy of the barrier parameter in sync; the corrector (and some
     * predictors) read it from the workspace. */
    ws.barr_param = barr_param;

    // t2 = clock();
    //    total_time += (double)(t2 - t1) / (double)clk_tck;

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0) {
      projerr = projectionError_norm_infinity_conic(velocity, reaction, n);
    } else
      projerr = projectionError(velocity, reaction, n, tol);

    setErrorArray(error, pinfeas, dinfeas, udotr, dualgap, complem, projerr);

    // check exit condition
    totalresidual = fmax(fmax(error[0], error[1]), error[2]);
    totalresidual_Jor = fmax(fmax(pinfeas, dinfeas), complem);

    // if (total_time > max_time)
    if (0 == 1) {
      // hasNotConverged = 5;
      // numerics_printf_verbose(1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e
      // |",
      //                         iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem,
      //                         projerr, barr_param);
      // printf("\nREQUEST TIMEOUT\n");
      // double unitur;
      // for (int i = 0; i < n; i++) {
      //   unitur = cblas_ddot(3, velocity + 3 * i, 1, reaction + 3 * i, 1);
      //   if (unitur < 0) printf("UR NEGATIF %9.2e\n", unitur);
      // }
      // break;
    }

    if (alpha_primal < 1e-10) {
      printf("\nSMALL Step-length\n");
      break;
    }

    // double worktmp[3];
    // if ( totalresidual <= tol )
    if (totalresidual_Jor <= tol)
    // if ( projerr <= tol )
    // if ( fmax(totalresidual,projerr) <= tol )
    // if ( fmax(fmax(totalresidual,projerr),complem) <= tol )
    {
      numerics_printf_verbose(1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                              iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem,
                              projerr, barr_param);

      double unitur;
      for (size_t i = 0; i < n; i++) {
        unitur = cblas_ddot(3, velocity + 3 * i, 1, reaction + 3 * i, 1);
        // unitur = cblas_ddot(3, u_plus_phiu+3*i, 1, reaction+3*i, 1);
        if (unitur < 0) printf("UR NEGATIF %9.2e\n", unitur);
      }

      hasNotConverged = 0;
      if (Qp) Qp = NM_free(Qp);
      if (Qpinv) Qpinv = NM_free(Qpinv);
      break;
    }

    /*  Build the Jacobian matrix without reducing the linear system.
     *
     *
     *  In the case where the NT scaling is not used, the matrix to factorize is the following:
     *
     *         m     nd       nd
     *      |  M     0      -H^T  | m
     *      |                     |
     *  J = |  0   Arw(r)  Arw(u) | nd
     *      |                     |
     *      | -H     I        0   | nd
     *
     *
     *  In the case where the NT scaling is used, the matrix to factorize is the following:
     *
     *         m     nd       nd
     *      |  M     0      -H^T  | m
     *      |                     |
     *  J = |  0     Qp2      I   | nd
     *      |                     |
     *      | -H     I        0   | nd
     *
     *  where Qp2 is the square of the quadratic representation of the Nesterov-Todd scaling
     * vector p
     *
     */
    /* Build the Jacobian matrix corresponding to a reduced linear system. The variable du
     (velocity) is eliminated.
     *
     * In this case, NT scaling is always performed, but two posible ways can be done to solve
     the linear system.
     *
     * In a first approach, the reduced matrix is of the following form:
     *
     *         m       nd
     *      |  -M      H^T    | m
     * JR = |                 |
     *      |   H     Qp^{-2} | nd
     *
     * where Qp^2 is the square of the quadratic representation of the vector p.
     *
     * A second possibility, which deals with a better contioned matrix, is to factorize the
     following reduced matrix:
     *
     *         m       nd
     *      |  -M     QpH^T  | m
     * JR = |                |
     *      |  QpH     I     | nd
     *
     *  where QpH = Qp * H.
     *
     *  The matrix QpH is computed by means of the function QNTpH.
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

       without NT scaling
       rhs = -
       [ M*v - H'*r - f ]  m         dualConstraint
       [     u o r      ]  nd        complemConstraint
       [  u - H*v - w   ]  nd        primalConstraint

       with NT scaling
       rhs = -
       [ M*v - H'*r - f ]  m         dualConstraint
       [       r        ]  nd        complemConstraint
       [  u - H*v - w   ]  nd        primalConstraint
    */

    /* Building the rhs for the first linear system in case of a reduced linear system  */
    /* In the case where Qp2 is inside the matrix, the rhs is the of the form           */
    /*         [ M*v - f - H'*r                ]  m         dualConstraint              */
    /* r_rhs = [                               ]                                        */
    /*         [ u - H*v - w - Qp_inv*r_check  ]  nd                                    */
    /*                                                                                  */
    /* In the case where QpH is inside the matrix, the rhs is the of the form           */
    /*        [ M*v - f - H'*r ]  m         dualConstraint                             */
    /* r_rhs = [                ]                                                       */
    /*         [ -Qp*(H*v + w)  ]  nd                                                   */

    int jacobian_is_nan = 0;
    ws.jacobian_is_nan = 0;
    switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
        gfc3d_IPM_predictor_3X3_NOSCAL(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
        gfc3d_IPM_predictor_3X3_QP2(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH:
        gfc3d_IPM_predictor_3X3_QPH(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2:
        gfc3d_IPM_predictor_2X2_QP2(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH:
        gfc3d_IPM_predictor_2X2_QPH(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
        gfc3d_IPM_predictor_1X1_QPH(&ws);
        break;

      default: {
        printf("ERROR 1\n");
      }
    }

    /* ---- sync main-loop locals from workspace ---- */
    jacobian_is_nan = ws.jacobian_is_nan;
    J = ws.J;
    JR = ws.JR;
    Qp2 = ws.Qp2;
    J_nzmax = ws.J_nzmax;
    JR_nzmax = ws.JR_nzmax;

    if (jacobian_is_nan) {
      hasNotConverged = 2;
      if (J) J = NM_free(J);
      if (JR) JR = NM_free(JR);
      break;
    }

    /* computing the affine step-length */
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    /* ----- Predictor step of Mehrotra ----- */
    cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
    cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0) {
      barr_param_a = complemResidualNorm(v_plus_dv, r_plus_dr, nd, n) / n;
      // barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n;
    } else
      barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n / 3.;

    /* computing the centralization parameter */
    e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(fmin(alpha_primal, alpha_dual), 2)) : sgmp3;
    sigma = fmin(1.0, pow(barr_param_a / barr_param, e));

    /* Computing the corrector step of Mehrotra algorithm */

    /*      without NT scaling */
    /*      rhs = - */
    /*      [         M*v - H'*r - f             ]  m         dualConstraint */
    /*      [ u o r + da^u o da^r - 2*sigma*mu*e ]  nd        complemConstraint */
    /*      [         u - H*v - w                ]  nd        primalConstraint */

    /*      with NT scaling */
    /*      rhs = - */
    /*      [ M*v - H'*r - f ]  m         dualConstraint */
    /*      [       r        ]  nd        complemConstraint */
    /*      [  u - H*v - w   ]  nd        primalConstraint */

    /* ---- sync intermediate values into workspace before corrector ---- */
    ws.alpha_primal = alpha_primal;
    ws.alpha_dual = alpha_dual;
    ws.gmm = gmm;
    ws.barr_param_a = barr_param_a;
    ws.e = e;
    ws.sigma = sigma;

    switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
        gfc3d_IPM_corrector_3X3_NOSCAL(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
        gfc3d_IPM_corrector_3X3_QP2(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH:
        gfc3d_IPM_corrector_3X3_QPH(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2:
        gfc3d_IPM_corrector_2X2_QP2(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH:
        gfc3d_IPM_corrector_2X2_QPH(&ws);
        break;

      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
        gfc3d_IPM_corrector_1X1_QPH(&ws);
        break;

      default: {
        printf("ERROR 2\n");
      }
    }

    /* ---- sync main-loop locals from workspace ---- */
    J = ws.J;
    JR = ws.JR;
    Qp2 = ws.Qp2;
    LS_norm_p = ws.LS_norm_p;
    LS_norm_d = ws.LS_norm_d;
    LS_norm_c = ws.LS_norm_c;

    if (Qp) {
      Qp = NM_free(Qp);
    }
    if (Qpinv) {
      Qpinv = NM_free(Qpinv);
    }

    // t1 = clock();

    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);
    alpha_tmp = alpha_dual;

    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length at the next iteration */

    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    // t2 = clock();
    // total_time += (double)(t2 - t1) / (double)clk_tck;

    numerics_printf_verbose(1,
                            "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e "
                            "| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr,
                            barr_param, alpha_primal, alpha_tmp, sigma,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]), LS_norm_p,
                            LS_norm_d, LS_norm_c);
    // numerics_printf_verbose(1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |
    // %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
    //                       iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem,
    //                       projerr, barr_param, alpha_primal, alpha_tmp, sigma,
    //                       cblas_dnrm2(m, d_globalVelocity, 1),
    //                       cblas_dnrm2(nd, d_velocity, 1),
    //                       cblas_dnrm2(nd, d_reaction, 1),
    //     LS_norm_p, LS_norm_d, LS_norm_c);

    /* ----- Update variables ----- */
    if (NV_isnan(d_globalVelocity, m) | NV_isnan(d_velocity, nd) | NV_isnan(d_reaction, nd)) {
      hasNotConverged = 2;
      break;
    }

    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction,
    //   d_globalVelocity, d_velocity, d_reaction, d, n, m, iterates);

    for (size_t i = 0; i < n; i++) {
      d_s[i] = cblas_ddot(2, velocity + i * d + 1, 1, d_velocity + i * d + 1, 1) /
                   cblas_dnrm2(2, velocity + i * d + 1, 1) -
               s[i] + cblas_dnrm2(2, velocity + i * d + 1, 1);
    }

    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);

    if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd)) {
      hasNotConverged = 2;
      break;
    }

    iteration++;

  }  // while loop
  free(d_s);
  /* Checking strict complementarity */
  /* For each cone i from 1 to n, one checks if u+r is in the interior of the Lorentz cone */
  /* One first computes the 3 dimensional vector somme = (u+r)/norm(u+r) */
  /* Then one checks if somme[0] > sqrt(somme[1]^2 + somme[2]^2) + ceps */
  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns;
  int nsc = 0;
  nN = nB = nR = 0;
  for (size_t i = 0; i < n; i++) {
    somme[0] = velocity[3 * i] + reaction[3 * i];
    somme[1] = velocity[3 * i + 1] + reaction[3 * i + 1];
    somme[2] = velocity[3 * i + 2] + reaction[3 * i + 2];
    dur[0] = velocity[3 * i] - reaction[3 * i];
    dur[1] = velocity[3 * i + 1] - reaction[3 * i + 1];
    dur[2] = velocity[3 * i + 2] - reaction[3 * i + 2];

    ns = somme[0] - cblas_dnrm2(2, somme + 1, 1);
    //    ndur = dur[0] - cblas_dnrm2(2, dur + 1, 1);
    if (ns > cesp * cblas_dnrm2(3, somme, 1)) {
      nsc += 1;
      if (dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        nN += 1;
      else if (-dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        nB += 1;
      else
        nR += 1;
    } else
      printf("cone %zu %9.2e %9.2e\n", i, somme[0], cblas_dnrm2(2, somme + 1, 1));
  }
  if (nsc < (int)n)
    printf("Ratio of Strict complementarity solutions: %4i / %4zu = %4.2f\n", nsc, n,
           (double)nsc / n);
  else
    printf("Strict complementarity satisfied: %4i / %4zu  %9.2e %9.2e  %4i %4i %4i\n", nsc, n,
           ns, ns / cblas_dnrm2(3, s, 1), nB, nN, nR);
  free(s);
  if (hasNotConverged == 0)
    printf("test IPM: success, ");
  else
    printf("test IPM: failure, ");

  printf("%.1e  %.1e,     iter=%3zu, n=%4zu, \t%s \n", totalresidual, projerr, iteration, n,
         strToken);

  printf("\n\n");
  // classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);

  /* ----- cleanup: restore scaling and free memory ------ */
  gfc3d_IPM_cleanup_data cd = {0};
  cd.problem = problem;
  cd.options = options;
  cd.internal_allocation = internal_allocation;
  cd.hasNotConverged = hasNotConverged;
  cd.restore_scaling = 1;
  cd.write_matlab_footer = 0;
  cd.m = m;
  cd.nd = nd;
  cd.iteration = iteration;
  cd.velocity = velocity;
  cd.reaction = reaction;
  cd.totalresidual = totalresidual;
  cd.M = M;
  cd.M_original = M_original;
  cd.P_mu = P_mu;
  cd.P_mu_inv = P_mu_inv;
  cd.H_tilde = H_tilde;
  cd.H = H;
  cd.minus_H = minus_H;
  cd.Minv = Minv;
  cd.minus_M = minus_M;
  cd.minus_Ht = minus_Ht;
  cd.eye_nd = eye_nd;
  cd.Ht = Ht;
  cd.iterates = iterates;
  cd.r_rhs = r_rhs;
  cd.r_rhs_2 = r_rhs_2;
  cd.r_dv = r_dv;
  cd.r_du = r_du;
  cd.r_dr = r_dr;
  cd.r_adu = r_adu;
  cd.r_adr = r_adr;
  cd.r_dr_a = r_dr_a;
  cd.r_du_a = r_du_a;
  cd.r_dv_a = r_dv_a;
  cd.Hvw = Hvw;
  cd.a_velo = a_velo;
  cd.a_reac = a_reac;
  cd.rhs_2 = rhs_2;
  cd.sr_rhs = sr_rhs;
  cd.sr_rhs_2 = sr_rhs_2;
  cd.phiu = phiu;
  cd.d_globalVelocity = d_globalVelocity;
  cd.d_velocity = d_velocity;
  cd.d_reaction = d_reaction;

  gfc3d_IPM_cleanup(&cd, info);
}

static void gfc3d_IPM_fixed_step(gfc3d_IPM_workspace* ws) {
  /* problem dimensions, handles, matrices, vectors and scalars */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const int type = ws->type;
  const double tol = ws->tol;
  SolverOptions* options = ws->options;
  GlobalFrictionContactProblem* problem = ws->problem;
  Gfc3d_IPM_init_data* data = ws->data;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H_tilde = ws->H_tilde;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* Ht = ws->Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* Minv = ws->Minv;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* JR = ws->JR;
  NumericsMatrix* Qp2 = ws->Qp2;
  double* w = ws->w;
  double* w_tilde = ws->w_tilde;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* v_plus_dv = ws->v_plus_dv;
  double* r_plus_dr = ws->r_plus_dr;
  double* dvdr_jprod = ws->dvdr_jprod;
  double* r_rhs = ws->r_rhs;
  double* r_rhs_2 = ws->r_rhs_2;
  double* sr_rhs = ws->sr_rhs;
  double* sr_rhs_2 = ws->sr_rhs_2;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* Hvw = ws->Hvw;
  double* p2 = ws->p2;
  double* velocity_t = ws->velocity_t;
  double* d_velocity_t = ws->d_velocity_t;
  double* d_reaction_t = ws->d_reaction_t;
  double* velocity_t_inv = ws->velocity_t_inv;
  double* d_s = ws->d_s;
  double* s = ws->s;
  double* phiu = ws->phiu;
  double* r_dv = ws->r_dv;
  double* r_du = ws->r_du;
  double* r_dr = ws->r_dr;
  double* r_dv_a = ws->r_dv_a;
  double* r_du_a = ws->r_du_a;
  double* r_dr_a = ws->r_dr_a;
  double* r_adu = ws->r_adu;
  double* r_adr = ws->r_adr;
  short* a_velo = ws->a_velo;
  short* a_reac = ws->a_reac;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* r_Qp_u = ws->r_Qp_u;
  double* r_Qp_u_inv = ws->r_Qp_u_inv;
  size_t iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  int jacobian_is_nan = ws->jacobian_is_nan;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double dualgap = ws->dualgap;
  double projerr = ws->projerr;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double totalresidual_mu = ws->totalresidual_mu;
  double diff_fixp = ws->diff_fixp;
  double barr_param = ws->barr_param;
  double barr_param_a = ws->barr_param_a;
  double sigma = ws->sigma;
  double e = ws->e;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double alpha_tmp = ws->alpha_tmp;
  double gmm = ws->gmm;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double norm_f = ws->norm_f;
  double norm_w = ws->norm_w;
  long J_nzmax = ws->J_nzmax;
  long JR_nzmax = ws->JR_nzmax;
  size_t H_nzmax = ws->H_nzmax;
  int nB = ws->nB;
  int nN = ws->nN;
  int nR = ws->nR;
  int nT = ws->nT;
  double nub = ws->nub;
  double unitur = ws->unitur;
  FILE* iterates = ws->iterates;
  char fws = ws->fws;
  char* str = ws->str;
  char* strToken = ws->strToken;

#pragma GCC diagnostic pop

  switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S]) {
    case SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_EXTERNAL: {
      if (iteration == 0) {
        for (size_t i = 0; i < nd; i += d) {
          nub = cblas_dnrm2(2, velocity + i + 1, 1);
          w[i] = w_tilde[i] + nub;
          s[i / d] = nub;
        }
      }

      primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas,
                     tol);  // primalConstraint = u - H*v -w
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas,
                   tol);  // dualConstraint =  M*v - H'*r - f
      dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
      JA_prod(velocity, reaction, nd, n, complemConstraint);
      for (size_t k = 0; k < nd; complemConstraint[k] -= 2. * barr_param, k += d);
      complem = cblas_dnrm2(nd, complemConstraint, 1);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1);
      projerr = projectionError(velocity, reaction, n, tol);
      for (size_t i = 0; i < nd; i += d) {
        nub = cblas_dnrm2(2, velocity + i + 1, 1);
        diff_fixp_vec[i / d] = s[i / d] - nub;
      }
      diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);

      totalresidual_mu = fmax(fmax(pinfeas, dinfeas), complem);
      totalresidual = fmax(totalresidual_mu, diff_fixp);

      if (totalresidual <= tol)  // SUCCESS => need to quit the main loop
      {
        numerics_printf_verbose(
            1, "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |", iteration,
            dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param);

        for (size_t i = 0; i < n; i++) {
          unitur = cblas_ddot(3, velocity + 3 * i, 1, reaction + 3 * i, 1);
          if (unitur < 0) printf("Cone %zu: UR NEGATIF ui'ri = %9.2e\n", i, unitur);
        }

        hasNotConverged = 0;
      }

      else if (totalresidual_mu <= tol)  // Finish solving ls (1 time) => update s after
      {
        numerics_printf_verbose(
            1, "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |", iteration,
            dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param);

        printf("\n");
        numerics_printf_verbose(1,
                                "| it  |  dualgap | pinfeas | dinfeas |  <u, r> |  |s-ub| | "
                                "complem | prj err | barpram |  2*n*mu |  alpha  |  |dv|   "
                                "|  |du|   |  |dr|   | ls prim | ls dual | ls comp |");

        // update s
        for (size_t i = 0; i < nd; i += d) {
          nub = cblas_dnrm2(2, velocity + i + 1, 1);
          w[i] = w_tilde[i] + nub;
          s[i / d] = nub;
        }
      }

      break;
    }

    case SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_INTERNAL: {
      // update s
      for (size_t i = 0; i < nd; i += d) {
        nub = cblas_dnrm2(2, velocity + i + 1, 1);
        w[i] = w_tilde[i] + nub;
        s[i / d] = nub;
      }
      diff_fixp = 0.;

      primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas,
                     tol);  // primalConstraint = u - H*v -w
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas,
                   tol);  // dualConstraint =  M*v - H'*r - f
      dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
      JA_prod(velocity, reaction, nd, n, complemConstraint);
      for (size_t k = 0; k < nd; complemConstraint[k] -= 2. * barr_param, k += d);
      complem = cblas_dnrm2(nd, complemConstraint, 1);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1);
      projerr = projectionError(velocity, reaction, n, tol);

      totalresidual = fmax(fmax(pinfeas, dinfeas), complem);

      if (totalresidual <= tol)  // SUCCESS => need to quit the main loop
      {
        numerics_printf_verbose(
            1, "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |", iteration,
            dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param);

        for (size_t i = 0; i < n; i++) {
          unitur = cblas_ddot(3, velocity + 3 * i, 1, reaction + 3 * i, 1);
          if (unitur < 0) printf("Cone %zu: UR NEGATIF ui'ri = %9.2e\n", i, unitur);
        }

        hasNotConverged = 0;
      }
      break;
    }

    default: {
      printf("gfc3d_IPM_fixed: update_s strategy is unknown. \n\n");
    }
  }

  /* ---- copy back modified state ---- */
  ws->hasNotConverged = hasNotConverged;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->dualgap = dualgap;
  ws->udotr = udotr;
  ws->projerr = projerr;
  ws->totalresidual = totalresidual;
  ws->totalresidual_mu = totalresidual_mu;
  ws->diff_fixp = diff_fixp;
  ws->nub = nub;
  ws->unitur = unitur;
}

void gfc3d_IPM_fixed(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
                     double* restrict velocity, double* restrict globalVelocity,
                     int* restrict info, SolverOptions* restrict options, double barr_param,
                     const char* problem_name) {
  size_t m = problem->M->size0;
  size_t nd = problem->H->size1;
  size_t d = problem->dimension;
  size_t n = problem->numberOfContacts;

  size_t no_m = 0, no_nd = 0;

  NumericsMatrix* M = NULL;
  NumericsMatrix* M_original = NULL;
  NumericsMatrix* H_tilde = NULL;
  NumericsMatrix* eye_nd = NM_eye(nd);

  /* symmetrization of the matrix M */
  if (!(NM_is_symmetric(problem->M))) {
    printf("#################### SYMMETRIZATION ####################\n");
    NumericsMatrix* MT = NM_transpose(problem->M);
    NumericsMatrix* MSym = NM_add(1 / 2., problem->M, 1 / 2., MT);
    M_original = problem->M; /* keep the user-provided matrix for restoration */
    problem->M = MSym;
    NM_free(MT);
  }

  M = problem->M;
  H_tilde = NM_transpose(problem->H);

  NumericsMatrix* minus_M = NM_copy(M);
  NM_scal(-1.0, minus_M);

  // initialize solver if it is not set
  int internal_allocation = 0;
  size_t work_size = m + nd + nd;
  if (!options->dWork || (options->dWorkSize < work_size)) {
    gfc3d_IPM_init(problem, options, work_size);
    internal_allocation = 1;
  }

  Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)options->solverData;
  NumericsMatrix* P_mu = data->P_mu->mat;
  // NumericsMatrix* P_mu_inv = data->P_mu->inv_mat;

  double* w_tilde = problem->b;
  double* w = data->tmp_vault_nd[no_nd++];
  double* f = problem->q;

  NumericsMatrix* H = NM_multiply(P_mu, H_tilde);
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;

  // Copy the initial data into executive vars
  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  // // computation of the global velocity vector: v = M\(H'*r+f)
  // for (size_t  i = 0; i<m; i++) globalVelocity[i] = f[i];
  // NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  // NM_Cholesky_solve(NM_preserve(M), globalVelocity, 1);

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  int iteration = 0;
  double pinfeas = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  double dualgap = 1e300;
  double udotr = 1e300;
  double projerr = 1e300;
  double unitur = 1e300;
  double totalresidual = 1e300, totalresidual_mu = 1e300;
  double diff_fixp = 1e300, nub = 1e300;
  double* diff_fixp_vec = (double*)calloc(n, sizeof(double));

  double* primalConstraint = data->tmp_vault_nd[no_nd++];
  double* dualConstraint = data->tmp_vault_m[no_m++];
  double* complemConstraint = data->tmp_vault_nd[no_nd++];

  double gmm = gmmp1 + gmmp2;
  // double barr_param_a, e;
  // double norm_f = cblas_dnrm2(m, f, 1);
  // double norm_w = cblas_dnrm2(nd, w, 1);

  double* d_globalVelocity = data->tmp_vault_m[no_m++];
  double* d_velocity = data->tmp_vault_nd[no_nd++];
  double* d_reaction = data->tmp_vault_nd[no_nd++];
  double* s = (double*)calloc(n, sizeof(double));

  double* rhs = options->dWork;
  double* rhs_2 = (double*)calloc(m + 2 * nd, sizeof(double));

  double* Hvw = data->tmp_vault_nd[no_nd++];
  double* r_Qp_u = data->tmp_vault_nd[no_nd++];
  double* r_Qp_u_inv = data->tmp_vault_nd[no_nd++];

  // double norm_q = cblas_dnrm2(m, problem->q, 1);
  // double norm_b = cblas_dnrm2(nd, problem->b, 1);

  /* norm of the residuals of teh second linear system */
  double LS_norm_d = 0;  // dual feaqsibility
  double LS_norm_c = 0;  // complementarity

  NumericsMatrix* J = NULL;
  long J_nzmax = 0;
  size_t H_nzmax = NM_nnz(H);

  numerics_printf_verbose(1, "LS solution: 2x2 NT scaling with QpH, mu fixed, s fixed\n");

  numerics_printf_verbose(
      1,
      "| it  |  dualgap | pinfeas | dinfeas |  <u, r> |  |s-ub| | complem | prj err | barpram "
      "|  2*n*mu |  alpha  |  |dv|   |  |du|   |  |dr|   | ls dual | ls comp |");
  numerics_printf_verbose(
      1,
      "---------------------------------------------------------------------------------------"
      "-----------------------------------------------------------------------");

  FILE* iterates = NULL;

  char* str = (char*)malloc(200);
  strcpy(str, problem_name);
  const char* separators = "/";
  char* strToken = strtok(str, separators);
  for (int i = 0; i < 5; i++) {
    if (strToken != NULL) strToken = strtok(NULL, separators);
  }
  strToken = strtok(strToken, ".");

  char matlab_name[100];

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
    snprintf(matlab_name, sizeof(matlab_name), "iterates_muF_sF.m");
    iterates = fopen(matlab_name, "a+");
    fprintf(iterates, "data(end+1).name = \"%s\";\n", strToken);
    fprintf(iterates, "data(end).val = [\n");
  }

  /* ---- Pack all workspace into a single struct for the fixed step function ---- */
  gfc3d_IPM_workspace ws = {0};
  ws.m = m;
  ws.nd = nd;
  ws.d = d;
  ws.n = n;
  ws.type = NORM_2;
  ws.tol = tol;
  ws.max_iter = (size_t)max_iter;
  ws.problem = problem;
  ws.options = options;
  ws.data = data;
  ws.gmmp1 = gmmp1;
  ws.gmmp2 = gmmp2;
  ws.gmm = gmm;
  ws.M = M;
  ws.H_tilde = H_tilde;
  ws.H = H;
  ws.minus_M = minus_M;
  ws.eye_nd = eye_nd;
  ws.J = J;
  ws.w_tilde = w_tilde;
  ws.w = w;
  ws.f = f;
  ws.reaction = reaction;
  ws.velocity = velocity;
  ws.globalVelocity = globalVelocity;
  ws.primalConstraint = primalConstraint;
  ws.dualConstraint = dualConstraint;
  ws.complemConstraint = complemConstraint;
  ws.d_globalVelocity = d_globalVelocity;
  ws.d_velocity = d_velocity;
  ws.d_reaction = d_reaction;
  ws.s = s;
  ws.rhs = rhs;
  ws.rhs_2 = rhs_2;
  ws.Hvw = Hvw;
  ws.r_Qp_u = r_Qp_u;
  ws.r_Qp_u_inv = r_Qp_u_inv;
  ws.diff_fixp_vec = diff_fixp_vec;
  ws.internal_allocation = internal_allocation;
  ws.hasNotConverged = hasNotConverged;
  ws.iteration = (size_t)iteration;
  ws.pinfeas = pinfeas;
  ws.dinfeas = dinfeas;
  ws.complem = complem;
  ws.udotr = udotr;
  ws.dualgap = dualgap;
  ws.projerr = projerr;
  ws.totalresidual = totalresidual;
  ws.totalresidual_mu = totalresidual_mu;
  ws.diff_fixp = diff_fixp;
  ws.barr_param = barr_param;
  ws.alpha_primal = alpha_primal;
  ws.alpha_dual = alpha_dual;
  ws.LS_norm_d = LS_norm_d;
  ws.LS_norm_c = LS_norm_c;
  ws.J_nzmax = J_nzmax;
  ws.H_nzmax = H_nzmax;
  ws.nub = nub;
  ws.unitur = unitur;
  ws.iterates = iterates;
  ws.str = str;
  ws.strToken = strToken;

  while (iteration < max_iter) {
    gfc3d_IPM_fixed_step(&ws);

    /* ---- sync main-loop locals from workspace ---- */
    hasNotConverged = ws.hasNotConverged;
    pinfeas = ws.pinfeas;
    dinfeas = ws.dinfeas;
    complem = ws.complem;
    udotr = ws.udotr;
    dualgap = ws.dualgap;
    projerr = ws.projerr;
    totalresidual = ws.totalresidual;
    totalresidual_mu = ws.totalresidual_mu;
    diff_fixp = ws.diff_fixp;
    nub = ws.nub;
    unitur = ws.unitur;

    if (hasNotConverged == 0) {
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
        fprintf(
            iterates,
            "%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
            iteration, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal,
            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]), LS_norm_d, LS_norm_c);
      }
      break;
    }

    J = NM_create(NM_SPARSE, m + nd, m + nd);
    J_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
    NM_triplet_alloc(J, J_nzmax);
    J->sparse->origin = NSM_TRIPLET;

    NumericsMatrix* QpH = QNTpH(velocity, reaction, H, nd, n);
    NumericsMatrix* QpHt = NM_transpose(QpH);

    NM_insert(J, minus_M, 0, 0);
    NM_insert(J, QpH, m, 0);
    NM_insert(J, QpHt, 0, m);
    NM_insert(J, eye_nd, m, m);

    if (QpH) QpH = NM_free(QpH);
    if (QpHt) QpHt = NM_free(QpHt);

    if (NM_isnan(J)) {
      numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
      hasNotConverged = 2;
      if (J) J = NM_free(J);
      break;
    }

    QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);  // r_Qp_u <- u_hat
    JA_inv(r_Qp_u, nd, n, r_Qp_u_inv);                   // r_Qp_u_inv <- u_hat_inv
    cblas_dscal(nd, -2. * barr_param, r_Qp_u_inv, 1);    // r_Qp_u_inv <= -2 * mu * u_hat_inv

    cblas_dcopy(m, dualConstraint, 1, rhs, 1);
    NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);
    cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);  // Hvw will be used for the computation of r_du
    QNTpz(velocity, reaction, Hvw, nd, n, rhs + m);
    cblas_daxpy(nd, 1.0, r_Qp_u_inv, 1, rhs + m, 1);
    cblas_dscal(nd, -1.0, rhs + m, 1);

    cblas_dcopy(m + nd, rhs, 1, rhs_2, 1);

    NSM_linearSolverParams(J)->solver = NSM_HSL;
    NM_LDLT_solve(J, rhs, 1);

    cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
    QNTpz(velocity, reaction, rhs + m, nd, n, d_reaction);
    NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
    cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
    cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u

    double* rhs_tmp = (double*)calloc(m + nd, sizeof(double));
    cblas_dcopy(m + nd, rhs_2, 1, rhs_tmp, 1);
    NM_gemv(1.0, J, rhs, -1.0, rhs_tmp);
    LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
    LS_norm_c = cblas_dnrm2(nd, rhs_tmp + m, 1);

    free(rhs_tmp);

    if (J) J = NM_free(J);

    /* computing the affine step-length */
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    numerics_printf_verbose(1,
                            "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | "
                            "%.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem,
                            projerr, barr_param, 2 * n * barr_param, alpha_primal,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]), LS_norm_d,
                            LS_norm_c);

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
      fprintf(iterates,
              "%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
              iteration, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal,
              fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
              fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
              fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]), LS_norm_d, LS_norm_c);
    }

    /* ----- Update variables ----- */

    if (NV_isnan(d_globalVelocity, m) | NV_isnan(d_velocity, nd) | NV_isnan(d_reaction, nd)) {
      hasNotConverged = 2;
      break;
    }
    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);

    iteration++;
  }  // while loop

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
    //    ndur = dur[0] - cblas_dnrm2(2, dur + 1, 1);
    if (ns > cesp * cblas_dnrm2(3, somme, 1)) {
      nsc += 1;
      if (dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        nN += 1;
      else if (-dur[0] >= cblas_dnrm2(2, dur + 1, 1))
        nB += 1;
      else
        nR += 1;
    } else
      printf("cone %zu %9.2e %9.2e\n", i, somme[0], cblas_dnrm2(2, somme + 1, 1));
  }
  if (nsc < (int)n)
    printf("Ratio of Strict complementarity solutions: %4i / %4zu = %4.2f\n", nsc, n,
           (double)nsc / n);
  else
    printf("Strict complementarity satisfied: %4i / %4zu  %9.2e %9.2e  %4i %4i %4i\n", nsc, n,
           ns, ns / cblas_dnrm2(3, s, 1), nB, nN, nR);

  /* ----- cleanup: save starting point and free memory ------ */
  gfc3d_IPM_cleanup_data cd = {0};
  cd.problem = problem;
  cd.options = options;
  cd.internal_allocation = internal_allocation;
  cd.hasNotConverged = hasNotConverged;
  cd.m = m;
  cd.nd = nd;
  cd.iteration = iteration;
  cd.velocity = velocity;
  cd.reaction = reaction;
  cd.globalVelocity = globalVelocity;
  cd.s = s;
  cd.diff_fixp_vec = diff_fixp_vec;
  cd.totalresidual = totalresidual;
  cd.M = M;
  cd.M_original = M_original;
  cd.H_tilde = H_tilde;
  cd.H = H;
  cd.minus_M = minus_M;
  cd.eye_nd = eye_nd;
  cd.iterates = iterates;
  cd.write_matlab_footer = 1;
  cd.restore_scaling = 0;

  gfc3d_IPM_cleanup(&cd, info);
}

/* Solver registration system */
#include "numerics_errors.h"
#include "solver_registry.h"

void gfc3d_ipm_set_default(SolverOptions* options) {
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = 0;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_EXTERNAL;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] =
      SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_INTERNAL;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
      SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH;

  // options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] =
  // SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_NO;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] =
  // SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] =
      SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES;

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09;  // 0.095
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers GFC3D_IPM in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int gfc3d_ipm_init_wrap(void* problem, SolverOptions* options) {
  /* set_default is now called by the registry before init_wrap */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int gfc3d_ipm_solve_wrap(void* problem, double* reaction, double* velocity,
                                double* globalVelocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  gfc3d_IPM((GlobalFrictionContactProblem*)problem, reaction, velocity, globalVelocity, &info,
            options);
  return info;
}

static void gfc3d_ipm_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER_3VAR(GFC3D_IPM, "GFC3D_IPM",
                     "Interior Point Method for 3D Global Friction Contact",
                     gfc3d_ipm_init_wrap, gfc3d_ipm_solve_wrap, gfc3d_ipm_free_wrap,
                     NULL,                  /* error function */
                     gfc3d_ipm_set_default, /* set_default */
                     1000,                  /* default_max_iter */
                     1e-8,                  /* default_tol */
                     0 /* is_local_solver */);
