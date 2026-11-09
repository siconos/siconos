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

#include "gfc3d_ipm_internal.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

#include "JordanAlgebra.h"
#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include "SiconosBlas.h"
#include "fc3d_short_names.h"
#include "gfc3d_ipm.h"
#include "safe_casts.h"

/* Returns the array of step-lengths to the boundary reduced by a factor gamma. Uses long
 * double.
 */
double* array_getStepLength(const double* const x, const double* const dx,
                            const size_t vecSize, const size_t varsCount, const double gamma) {
  size_t dimension = vecSize / varsCount;
  size_t pos;
  float_type aL, bL, cL, dL, alphaL;
  double* arr_alpha = (double*)calloc(varsCount, sizeof(double));

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

    arr_alpha[i] = alphaL;
  }

  return arr_alpha;
}

/* Returns the primal constraint vector for global fricprob: out = velocity - H x
 * globalVelocity - w - phi(s) and the relative 2-norm of this vector: |out|/max{|velocity|, |H
 * x globalVelocity|, |w|, |phi(s)|} */
void primalResidual_s(const double* velocity, NumericsMatrix* H, const double* globalVelocity,
                      const double* w, const double* s, double* out, double* rnorm,
                      const double tol) {
  const size_t nd = to_size_t(H->size0);
  const blasint blas_nd = to_blasint(nd);
  double rn;

  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  rn = cblas_dnrm2(blas_nd, out, 1);
  cblas_daxpy(blas_nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(blas_nd, -1.0, w, 1, out, 1);

  for (size_t i = 0; i < nd; i += 3) out[i] -= s[i / 3];

  rn = fmax(rn, cblas_dnrm2(blas_nd, velocity, 1));
  rn = fmax(rn, cblas_dnrm2(blas_nd, w, 1));
  rn = fmax(rn, cblas_dnrm2(blas_nd / 3, s, 1));
  *rnorm = (rn > tol ? cblas_dnrm2(blas_nd, out, 1) : cblas_dnrm2(blas_nd, out, 1));
}

void primalResidual_s_type(const double* velocity, NumericsMatrix* H,
                           const double* globalVelocity, const double* w, const double* s,
                           double* out, double* rnorm, const double tol, const int type) {
  const size_t nd = to_size_t(H->size0);
  const blasint blas_nd = to_blasint(nd);
  double rn;

  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  rn = NV_norm_type(nd, out, type);
  cblas_daxpy(blas_nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(blas_nd, -1.0, w, 1, out, 1);

  for (size_t i = 0; i < nd; i += 3) out[i] -= s[i / 3];

  rn = fmax(rn, NV_norm_type(nd, velocity, type));
  rn = fmax(rn, NV_norm_type(nd, w, type));
  rn = fmax(rn, NV_norm_type(nd / 3, s, type));
  *rnorm = (rn > tol ? NV_norm_type(nd, out, type) : NV_norm_type(nd, out, type));
}
