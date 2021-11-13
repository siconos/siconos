/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include "JordanAlgebra.h"
#include "cblas.h"
#include "math.h"
#include "NumericsVector.h"
#include <float.h>               // for DBL_EPSILON, DBL_MAX

//#define DEBUG_MESSAGES
#include "siconos_debug.h"


NumericsMatrix* Arrow_repr(const double* const vec, const unsigned int vecSize, const size_t varsCount)
{
  /* validation */
  if(vecSize % varsCount != 0)
  {
    fprintf(stderr, "Arrow_repr: %zu variables can not be extracted from vector of size %d.\n", varsCount, vecSize);
    exit(EXIT_FAILURE);
  }

  size_t dimension = (size_t)(vecSize / varsCount);
  if(dimension < 2)
  {
    fprintf(stderr, "Arrow_repr: The dimension of variables can not be less than 2 but given %zu.\n", dimension);
    exit(EXIT_FAILURE);
  }

  NumericsMatrix * Arw_mat = NM_create(NM_SPARSE, vecSize, vecSize);
  size_t nzmax = (dimension * 3 - 2) * varsCount;
  NM_triplet_alloc(Arw_mat, nzmax);
  NM_fill(Arw_mat, NM_SPARSE, vecSize, vecSize, Arw_mat->matrix2);

  /* Arrow matirx filling */
  size_t pos;
  for(size_t i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;
    NM_entry(Arw_mat, pos, pos, vec[pos]);

    for(size_t j = 1; j < dimension; ++j)
    {
      NM_entry(Arw_mat, pos, pos + j, vec[pos + j]);
      NM_entry(Arw_mat, pos + j, pos, vec[pos + j]);
      NM_entry(Arw_mat, pos + j, pos + j, vec[pos]);
    }
  }
  return Arw_mat;
}


NumericsMatrix* Reflect_mat(const unsigned int size, NM_types type)
{
  NumericsMatrix * Refl_mat = NM_create(type, size, size);

  if(type == NM_SPARSE)
  {
    NM_triplet_alloc(Refl_mat, size);
    NM_fill(Refl_mat, NM_SPARSE, size, size, Refl_mat->matrix2);
  }

  NM_entry(Refl_mat, 0, 0, 1.0);
  for(unsigned int i = 1; i < size; ++i)
    NM_entry(Refl_mat, i, i, -1.0);
  return Refl_mat;
}

NumericsMatrix* Quad_repr(const double* const vec, const unsigned int vecSize, const size_t varsCount)
{
  size_t dimension = (size_t)(vecSize / varsCount);
  NumericsMatrix* out = NM_create(NM_SPARSE, vecSize, vecSize);
  NM_triplet_alloc(out, dimension * dimension * varsCount);
  //NM_fill(out, NM_SPARSE, vecSize, vecSize, out->matrix2);

  NumericsMatrix* quad_tmp = NM_create(NM_DENSE, dimension, dimension);

  NumericsMatrix* R = Reflect_mat(dimension, NM_DENSE);

  double * dets = (double*)malloc(varsCount * sizeof(double));
  JA_det(vec, vecSize, varsCount, dets);

  for(int i = 0; i < vecSize; i += dimension)
  {
    NV_dott(vec + i, vec + i, dimension, quad_tmp);

    for(int j = 0; j < dimension; ++j)
    {
      for(int k = 0; k < dimension; ++k)
        quad_tmp->matrix0[j+k*quad_tmp->size0] *= 2.0;
      NM_entry(quad_tmp, j, j, NM_get_value(quad_tmp, j, j) - dets[(int)(i / dimension)] * NM_get_value(R, j, j));
    }
    NM_insert(out, quad_tmp, i, i);
  }

  NM_clear(R);
  NM_clear(quad_tmp);
  free(R);
  free(quad_tmp);
  free(dets);

  return out;
}


void NesterovToddVector(const double* const vec1, const double* const vec2,
                        const unsigned int vecSize, const size_t varsCount, double * out)
{
  DEBUG_BEGIN("NesterovToddVector(...)\n");
  NumericsMatrix* quad_repr;
  double * x05 = (double*)calloc(vecSize, sizeof(double));
  double * Qx05y = (double*)calloc(vecSize, sizeof(double));
  double * Qx05yi = (double*)calloc(vecSize, sizeof(double));
  double * _p = (double*)calloc(vecSize, sizeof(double));

  assert(!NV_isnan(vec1,vecSize));
  assert(!NV_isnan(vec2,vecSize));

  JA_sqrt(vec1, vecSize, varsCount, x05);
  quad_repr = Quad_repr(x05, vecSize, varsCount);
  NM_gemv(1.0, quad_repr, vec2, 0.0, Qx05y);
  JA_sqrt_inv(Qx05y, vecSize, varsCount, Qx05yi);
  NM_gemv(1.0, quad_repr, Qx05yi, 0.0, _p);
  JA_sqrt_inv(_p, vecSize, varsCount, out);

  //assert(!NV_isnan(out,vecSize));


  free(_p);
  free(Qx05yi);
  free(Qx05y);
  free(x05);
  DEBUG_END("NesterovToddVector(...)\n");
}


double * JA_iden(const unsigned int vecSize, const size_t varsCount)
{
  if(vecSize % varsCount != 0)
  {
    fprintf(stderr, "JA_iden: %zu variables can not be extracted from vector of size %u.\n", varsCount, vecSize);
    exit(EXIT_FAILURE);
  }
  double * out = (double*)calloc(vecSize, sizeof(double));
  unsigned int dimension = (int)(vecSize / varsCount);
  for(unsigned int i = 0; i < vecSize; i += dimension)
    out[i] = 1.0;
  return out;
}


/* Jordan product of two vectors */
void JA_prod(const double * const vec1, const double * const vec2, const unsigned int vecSize, const int varsCount, double * out)
{
  assert(vec1);
  assert(vec2);

  unsigned int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  for(int i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;
    out[pos] = cblas_ddot(dimension, vec1 + pos, 1, vec2 + pos, 1);
    for(unsigned int j = 1; j < dimension; j++)
      out[pos + j] = vec1[pos] * vec2[pos + j] + vec2[pos] * vec1[pos + j];
  }
}


/* Returns the eigenvalues of each element in the vector. */
void JA_eigenvals(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out)
{
  DEBUG_BEGIN("JA_eigenvals(...)\n");
  unsigned int dimension = (int)(vecSize / varsCount);
  register unsigned int pos;

  for(size_t i = 0; i < 2*varsCount; i += 2)
  {
    pos = (i / 2.) * dimension;
    //assert(!isnan(vec[pos]));
    out[i] = vec[pos] + NV_norm_2(vec + pos + 1, dimension - 1);
  }

  for(size_t i = 1; i < 2*varsCount; i += 2)
  {
    pos = ((i - 1) / 2.) * dimension;
    out[i] = vec[pos] - NV_norm_2(vec + pos + 1, dimension - 1);
  }
  DEBUG_EXPR(NV_display(vec,vecSize););
  DEBUG_EXPR(NV_display(out,2*varsCount););
  DEBUG_END("JA_eigenvals(...)\n");
}

void JA_eigenvecs(const double * const vec, const unsigned int vecSize, const size_t varsCount, double ** out)
{
  //const double EPS = 1e-12;
  unsigned int dimension = (int)(vecSize / varsCount);
  register unsigned int pos;
  double xi_bar_norm;

  for(size_t i = 0; i < 2*varsCount; i += 2)
  {
    pos = (i / 2.) * dimension;
    xi_bar_norm = NV_norm_2(vec + pos + 1, dimension - 1);
    out[i][0] = 0.5;
    //NV_const_add(vec + pos + 1, dimension - 1, 1. / (2 * xi_bar_norm + EPS), 0, out[i] + 1);
    NV_const_add(vec + pos + 1, dimension - 1, 1. / (2 * xi_bar_norm), 0, out[i] + 1);
  }

  for(size_t i = 1; i < 2*varsCount; i += 2)
  {
    pos = ((i - 1) / 2.) * dimension;
    xi_bar_norm = NV_norm_2(vec + pos + 1, dimension - 1);
    out[i][0] = 0.5;
    //NV_const_add(vec + pos + 1, dimension - 1, - 1. / (2 * xi_bar_norm + EPS), 0, out[i] + 1);
    NV_const_add(vec + pos + 1, dimension - 1, - 1. / (2 * xi_bar_norm), 0, out[i] + 1);
  }
}

void JA_sqrt(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int pos;
  unsigned int dimension = (int)(vecSize / varsCount);
  double * eigenvals = (double*)malloc(2 * varsCount * sizeof(double));
  double ** eigenvecs = (double**)malloc(2 * varsCount * sizeof(double*));
  for(unsigned int i = 0; i < 2 * varsCount; ++i)
    eigenvecs[i] = (double*)calloc(dimension, sizeof(double));

  double *tmp_vec1 = (double*)malloc(dimension * sizeof(double));
  double *tmp_vec2 = (double*)malloc(dimension * sizeof(double));
  double sqrt_eigenval1, sqrt_eigenval2;

  JA_eigenvals(vec, vecSize, varsCount, eigenvals);
  JA_eigenvecs(vec, vecSize, varsCount, eigenvecs);

  for(size_t i = 0; i < 2 * varsCount; i += 2)
  {
    assert(eigenvals[i] >=0);
    assert(eigenvals[i+1] >=0);
    sqrt_eigenval1 = sqrt(eigenvals[i]);
    sqrt_eigenval2 = sqrt(eigenvals[i + 1]);
    pos = (i / 2) * dimension;
    NV_const_add(eigenvecs[i], dimension, sqrt_eigenval1, 0, tmp_vec1);
    NV_const_add(eigenvecs[i + 1], dimension, sqrt_eigenval2, 0, tmp_vec2);
    NV_add(tmp_vec1, tmp_vec2, dimension, out + pos);
  }

  free(eigenvals);
  for(size_t i = 0; i < 2 * varsCount; ++i)
    free(eigenvecs[i]);
  free(eigenvecs);
  free(tmp_vec1);
  free(tmp_vec2);
}


void JA_sqrt_inv(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int pos;
  unsigned int dimension = (int)(vecSize / varsCount);
  double * eigenvals = (double*)malloc(2 * varsCount * sizeof(double));
  double ** eigenvecs = (double**)malloc(2 * varsCount * sizeof(double*));
  for(size_t i = 0; i < 2 * varsCount; ++i)
    eigenvecs[i] = (double*)calloc(dimension, sizeof(double));

  double *tmp_vec1 = (double*)malloc(dimension * sizeof(double));
  double *tmp_vec2 = (double*)malloc(dimension * sizeof(double));
  double sqrt_eigenval1, sqrt_eigenval2;

  JA_eigenvals(vec, vecSize, varsCount, eigenvals);
  JA_eigenvecs(vec, vecSize, varsCount, eigenvecs);

  for(size_t i = 0; i < 2 * varsCount; i += 2)
  {
    //printf("eigenvals[i] = %e \t, eigenvals[i+1] = %e\n ", eigenvals[i], eigenvals[i + 1]);
    //assert(eigenvals[i]>0);
    //assert(eigenvals[i+1]>0);
    sqrt_eigenval1 = 1. / sqrt(eigenvals[i]);
    sqrt_eigenval2 = 1. / sqrt(eigenvals[i + 1]);
    pos = (i / 2) * dimension;
    NV_const_add(eigenvecs[i], dimension, sqrt_eigenval1, 0, tmp_vec1);
    NV_const_add(eigenvecs[i + 1], dimension, sqrt_eigenval2, 0, tmp_vec2);
    NV_add(tmp_vec1, tmp_vec2, dimension, out + pos);
  }

  free(eigenvals);
  for(size_t i = 0; i < 2 * varsCount; ++i)
    free(eigenvecs[i]);
  free(eigenvecs);
  free(tmp_vec1);
  free(tmp_vec2);
}

void JA_power2(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int pos;
  unsigned int dimension = (int)(vecSize / varsCount);
  double * eigenvals = (double*)malloc(2 * varsCount * sizeof(double));
  double ** eigenvecs = (double**)malloc(2 * varsCount * sizeof(double*));
  for(size_t i = 0; i < 2 * varsCount; ++i)
    eigenvecs[i] = (double*)calloc(dimension, sizeof(double));

  double *tmp_vec1 = (double*)malloc(dimension * sizeof(double));
  double *tmp_vec2 = (double*)malloc(dimension * sizeof(double));
  double sqrt_eigenval1, sqrt_eigenval2;

  JA_eigenvals(vec, vecSize, varsCount, eigenvals);
  JA_eigenvecs(vec, vecSize, varsCount, eigenvecs);

  for(size_t i = 0; i < 2 * varsCount; i += 2)
  {
    sqrt_eigenval1 = eigenvals[i] * eigenvals[i];
    sqrt_eigenval2 = eigenvals[i + 1] * eigenvals[i + 1];
    pos = (i / 2) * dimension;
    NV_const_add(eigenvecs[i], dimension, sqrt_eigenval1, 0, tmp_vec1);
    NV_const_add(eigenvecs[i + 1], dimension, sqrt_eigenval2, 0, tmp_vec2);
    NV_add(tmp_vec1, tmp_vec2, dimension, out + pos);
  }

  free(eigenvals);
  for(size_t i = 0; i < 2 * varsCount; ++i)
    free(eigenvecs[i]);
  free(eigenvecs);
  free(tmp_vec1);
  free(tmp_vec2);
}

void JA_inv(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int pos;
  unsigned int dimension = (int)(vecSize / varsCount);
  double * eigenvals = (double*)malloc(2 * varsCount * sizeof(double));
  double ** eigenvecs = (double**)malloc(2 * varsCount * sizeof(double*));
  for(size_t i = 0; i < 2 * varsCount; ++i)
    eigenvecs[i] = (double*)calloc(dimension, sizeof(double));

  double *tmp_vec1 = (double*)malloc(dimension * sizeof(double));
  double *tmp_vec2 = (double*)malloc(dimension * sizeof(double));
  double sqrt_eigenval1, sqrt_eigenval2;

  JA_eigenvals(vec, vecSize, varsCount, eigenvals);
  JA_eigenvecs(vec, vecSize, varsCount, eigenvecs);

  for(size_t i = 0; i < 2 * varsCount; i += 2)
  {
    sqrt_eigenval1 = 1. / eigenvals[i];
    sqrt_eigenval2 = 1. / eigenvals[i + 1];
    pos = (i / 2) * dimension;
    NV_const_add(eigenvecs[i], dimension, sqrt_eigenval1, 0, tmp_vec1);
    NV_const_add(eigenvecs[i + 1], dimension, sqrt_eigenval2, 0, tmp_vec2);
    NV_add(tmp_vec1, tmp_vec2, dimension, out + pos);
  }

  free(eigenvals);
  for(size_t i = 0; i < 2 * varsCount; ++i)
    free(eigenvecs[i]);
  free(eigenvecs);
  free(tmp_vec1);
  free(tmp_vec2);
}

void JA_det(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  for(size_t i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;
    out[i] = vec[pos] * vec[pos] - cblas_ddot(dimension - 1, vec + pos + 1, 1, vec + pos + 1, 1);
  }
}

/* Returns the 2-norm of a vector - uses long double - based on blas_dnrm2 */
long double dnrm2l(const unsigned int n, const double * x)
{
  long double norm, scale, ssq, absxi, quo;

  if (n < 1)
    norm = 0.0;
  else if (n == 1)
    norm = fabs(x[0]);
  else
  {
    scale = 0.0;
    ssq = 1.0;
    for (size_t i = 0; i < n; i++)
    {
      if (x[i] != 0)
      {
        absxi = fabsl(x[i]);
        if (scale < absxi)
        {
          quo = scale/absxi;
          ssq = 1.0 + ssq * (quo * quo);
          scale = absxi;
        }
        else
        {
          quo = absxi/scale;
          ssq = ssq + (quo * quo);
        }
      }
      norm = scale * sqrtl(ssq);
    }
  }
  return norm;
}

/* Returns the step length for variables update in IPM */
double getNewtonStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
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
    //   bi = gamma*bi;

    NV_power2(xi, dimension, xi2);
    ci = xi2[0] - NV_reduce((xi2 + 1), dimension - 1);
    //    ci = gamma*gamma*ci;

    di = bi * bi - ai * ci;

    if(ai < 0 || (bi < 0 && ai < (bi * bi) / ci))
      alpha = ((-bi - sqrt(di)) / ai);
    else if((fabs(ai) < DBL_EPSILON) && (bi < 0))
      alpha = (-ci / (2 * bi));
    else
      alpha = DBL_MAX;
    //NV_display(xi2, dimension);
    //printf("**************** %3i ai = %9.2e b = %9.2e ci = %9.2e alpha = %9.2e\n",i, ai, bi, ci, alpha);

    if(fabs(alpha) < DBL_EPSILON)
      alpha = 0.0;

    alpha_list[i] = alpha;
  }

  min_alpha = NV_min(alpha_list, varsCount);

  free(xi2);
  free(dxi2);
  free(xi_dxi);
  free(alpha_list);

  return fmin(1.0, gamma * min_alpha);
  //return fmin(1.0, min_alpha);
}

/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double. */
double getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                            const unsigned int varsCount, const double gamma)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  long double aL, bL, cL, dL, alphaL;
  double min_alpha;

  min_alpha = 1e20; //1.0;

  for(unsigned int i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;
    aL = dnrm2l(dimension-1, dx+pos+1);
    aL = (dx[pos] - aL)*(dx[pos] + aL);
    bL = x[pos]*dx[pos];
    for (int k = 1; k < dimension; bL -= x[pos+k]*dx[pos+k], k++);
    cL = dnrm2l(dimension-1, x+pos+1);
    cL = (x[pos] - cL)*(x[pos] + cL);
    dL = bL*bL - aL*cL;
    if(aL < 0 || (bL < 0 && dL > 0 ))
      if (bL>0)
        alphaL = -(bL+sqrtl(dL))/aL;
      else
        alphaL = cL/(-bL+sqrtl(dL));
    else if((fabsl(aL) == 0.0) && (bL < 0))
      alphaL = -cL/bL/2;
    else
      alphaL = DBL_MAX;
    min_alpha = ((alphaL < min_alpha) ? alphaL : min_alpha);
  }
  min_alpha = gamma*min_alpha;
  min_alpha = ((min_alpha < 1.0) ? min_alpha : 1.0);
  return min_alpha;
}

/* Returns the primal constraint vector for global fricprob: out = velocity - H x globalVelocity - w */
/* and the relative 2-norm of this vector: |out|/max{|velocity|, |H x globalVelocity|, |w|} */
void primalResidual(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                           double * out, double * rnorm)
{
  size_t nd = H->size0;
  double rn;


  /* The memory for the result vectors should be allocated using calloc
   * since H is a sparse matrix. In other case the behaviour will be undefined.*/
  //  double *Hv = (double*)calloc(nd, sizeof(double));
  //double *u_minus_Hv = (double*)calloc(nd, sizeof(double));

  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  /* rn = cblas_dnrm2(nd, out, 1); */
  cblas_daxpy(nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(nd, -1.0, w, 1, out, 1);
  /* rn = fmax(rn, cblas_dnrm2(nd, velocity, 1)); */
  /* rn = fmax(rn, cblas_dnrm2(nd, w, 1)); */
  /* *rnorm = (rn > DBL_EPSILON ? cblas_dnrm2(nd, out, 1)/rn : cblas_dnrm2(nd, out, 1)); */
  *rnorm = cblas_dnrm2(nd, out, 1);
}

/* Returns the dual constraint vector for global fricprob ( M*globalVelocity - f - H'*reaction ) */
void dualResidual(NumericsMatrix * M, const double * globalVelocity, NumericsMatrix * H, const double * reaction, const double * f,
                         double * out, double * rnorm )
{
  double m = H->size1;
  double *HTr = (double*)calloc(m, sizeof(double));
  double rn;

  NM_gemv(1.0, M, globalVelocity, 0.0, out);
  /* rn = cblas_dnrm2(m, out, 1); */
  cblas_daxpy(m, -1.0, f, 1, out, 1);
  NM_tgemv(1.0, H, reaction, 0.0, HTr);
  cblas_daxpy(m, -1.0, HTr, 1, out, 1);
  /* rn = fmax(rn, cblas_dnrm2(m, f, 1)); */
  /* rn = fmax(rn, cblas_dnrm2(m, HTr, 1)); */
  /* *rnorm = (rn >DBL_EPSILON ? cblas_dnrm2(m, out, 1)/rn : cblas_dnrm2(m, out, 1)); */
  *rnorm = cblas_dnrm2(m, out, 1);
  free(HTr);
}

/* Returns the 2-norm of primal residual vector = | H * globalVelocity + w - velocity |_2 / (1 + |w|_inf) */
/* double primalResidualNorm(const double * velocity, NumericsMatrix * H, */
/*                                  const double * globalVelocity, const double * w) */
/* { */
/*   double * resid = (double*)calloc(H->size0, sizeof(double)); */
/*   primalResidualVector(velocity, H, globalVelocity, w, resid); */
/*   double norm_2 = cblas_dnrm2(H->size0, resid, 1); */
/*   free(resid); */
/*   return norm_2 / (1 + NV_norm_inf(w, H->size0)); */
/* } */

/* Returns the 2-norm of the dual residual vector  = | M * globalVelocity - H * reaction + f |_2 / (1 + |f|_inf)  */
/* double dualResidualNorm(NumericsMatrix * M, const double * globalVelocity, */
/*                                NumericsMatrix * H, const double * reaction, const double * f) */
/* { */
/*   double * resid = (double*)calloc(H->size1, sizeof(double)); */
/*   dualResidualVector(M, globalVelocity, H, reaction, f, resid); */
/*   double norm_2 = cblas_dnrm2(H->size1, resid, 1); */
/*   free(resid); */
/*   return norm_2 / (1 + NV_norm_inf(f, H->size1)); */
/* } */

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product velocity o reaction  */
double complemResidualNorm(const double * const velocity, const double * const reaction,
                                  const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  JA_prod(velocity, reaction, vecSize, varsCount, resid);
  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  free(resid);
  return norm2;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
double complemResidualNorm_p(const double * const velocity, const double * const reaction,
                                    const unsigned int vecSize, const unsigned int varsCount)
{

  double * resid = (double*)calloc(vecSize, sizeof(double));
  double * u_p = (double*)calloc(vecSize, sizeof(double));
  double * r_p = (double*)calloc(vecSize, sizeof(double));
  //double * p_inv = (double*)calloc(vecSize, sizeof(double));
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(velocity, reaction, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(velocity, b, vecSize, varsCount, a);

  Qx50y(a, velocity, vecSize, varsCount, u_p);
  Qx05y(a, reaction, vecSize, varsCount, r_p);
  JA_prod(u_p, r_p, vecSize, varsCount, resid);

  double norm2 = cblas_dnrm2(vecSize, resid, 1);

  /* Qxy(p, velocity, vecSize, varsCount, u_p); */
  /* JA_inv(p, vecSize, varsCount, p_inv); */
  /* Qxy(p_inv, reaction, vecSize, varsCount, r_p); */
  /* JA_prod(u_p, r_p, vecSize, varsCount, resid); */

  /* norm2 = cblas_dnrm2(vecSize, resid, 1); */
  /* printf("complemnt-2 = %.15e\n",norm2); */
  //norm2 = norm2/ cblas_dnrm2(vecSize, velocity, 1);
  //norm2 = norm2 / cblas_dnrm2(vecSize, reaction, 1);
  free(resid);
  free(u_p);
  free(r_p);
  //free(p_inv);
  free(a);
  free(b);

  return norm2;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
/* This computation is done with the formula "F" */
double complemResidualNorm_p_F(NumericsMatrix * Qp, NumericsMatrix * Qpinv,
              const double * const velocity, const double * const reaction,
                                      const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  double * u_p = (double*)calloc(vecSize, sizeof(double));
  double * r_p = (double*)calloc(vecSize, sizeof(double));

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
double dualGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m)
{
  double * Mv = (double*)calloc(m, sizeof(double));
  double vMv, pval, dval;

  NM_gemv(0.5, M, globalVelocity, 0.0, Mv);
  vMv = cblas_ddot(m, globalVelocity, 1, Mv, 1);
  free(Mv);
  pval = vMv - cblas_ddot(m, f, 1, globalVelocity, 1);
  dval = -vMv - cblas_ddot(nd, w, 1, reaction, 1);
  return (pval - dval)/ (1 + fabs(pval) + fabs(dval));
}

void setErrorArray(double * error, const double pinfeas, const double dinfeas,
                          const double dualgap, const double complem, const double complem_p)
{
  error[0] = pinfeas;
  error[1] = dinfeas;
  error[2] = dualgap;
  error[3] = complem_p;
  error[4] = complem;
}

/* Return the 2-norm of the difference between two vectors */
double norm2VecDiff (const double * vec1, const double * vec2, const unsigned int vecSize)
{
  double *vecDiff;
  double nvd;
  vecDiff = (double*)calloc(vecSize,sizeof(double));
  cblas_dcopy(vecSize, vec1, 1, vecDiff, 1);
  cblas_daxpy(vecSize, -1.0, vec2, 1, vecDiff, 1);
  nvd = cblas_dnrm2(vecSize, vecDiff, 1);
  free(vecDiff);
  return nvd;
}

/* Returns the product Q_sqrt(x)*y */
void Qx05y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, c1y, c2y, nxb, fx1, fx2, dx;
  size_t j;
  long double *xb = (long double*)calloc(dimension-1, sizeof(long double));

  for (int i = 0; i < dimension - 1; xb[i] = 1/sqrtl(dimension-1), i++);

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    nxb = dnrm2l(dimension-1, x+j+1);
    if (nxb > 0)
      for (int k = 0; k < dimension-1; xb[k] = x[j+1+k]/nxb, k++);
    l1 = x[j]+nxb;
    l2 = x[j]-nxb;
    dx = sqrtl(l1*l2);
    c1y = y[j];
    for (int k = 0; k < dimension-1; c1y += xb[k]*y[j+1+k], k++);
    c2y = 2*y[j] - c1y;
    fx1 = (l1*c1y + dx*c2y)/2;
    fx2 = (dx*c1y + l2*c2y)/2;
    out[j] = fx1 + fx2 - dx*y[j];
    for (int k = 0; k < dimension-1; out[j+k+1] = fx1*xb[k] - fx2*xb[k] + dx*y[j+k+1], k++);
  }
  free(xb);
}

/* Returns the product Q_inv_sqrt(x)*y */
void Qx50y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, c1y, c2y, nxb, fx1, fx2, dx;
  size_t j;
  long double *xb = (long double*)calloc(dimension-1, sizeof(long double));

  for (int i = 0; i < dimension - 1; xb[i] = 1/sqrtl(dimension-1), i++);

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    nxb = dnrm2l(dimension-1, x+j+1);
    if (nxb > 0)
      for (int k = 0; k < dimension-1; xb[k] = x[j+1+k]/nxb, k++);

    l1 = x[j]+nxb;
    l2 = x[j]-nxb;
    dx = 1/sqrtl(l1*l2);
    c1y = y[j];
    for (int k = 0; k < dimension-1; c1y += xb[k]*y[j+1+k], k++);
    c2y = 2*y[j] - c1y;
    fx1 = (c1y/l1 + dx*c2y)/2;
    fx2 = (dx*c1y + c2y/l2)/2;
    out[j] = fx1 + fx2 - dx*y[j];
    for (int k = 0; k < dimension-1; out[j+k+1] = fx1*xb[k] - fx2*xb[k] + dx*y[j+k+1], k++);
  }
  free(xb);
}
/* void Qx50y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out) */
/* { */
/*   unsigned int dimension = (int)(vecSize / varsCount); */
/*   double l1, l2, c1y, c2y, normx, fx1, fx2, dx;  */
/*   size_t j; */
/*   double *xn = (double*)calloc(dimension-1, sizeof(double)); */

/*   for (int i = 0; i < dimension - 1; xn[i] = 1/sqrt(dimension-1), i++); */

/*   for(size_t i = 0; i < varsCount; i++) */
/*     { */
/*       j = i*dimension; */
/*       normx = cblas_dnrm2(dimension-1, x+j+1, 1); */
/*       if (normx > 0) */
/*  { */
/*    cblas_dcopy(dimension-1, x+j+1, 1, xn, 1); */
/*    cblas_dscal(dimension-1, 1.0/normx, xn, 1); */
/*  } */
/*       cblas_dcopy(dimension-1, x+j+1, 1, xn, 1); */
/*       cblas_dscal(dimension-1, 1.0/normx, xn, 1); */
/*       l1 = 1/(x[j]+normx); */
/*       l2 = 1/(x[j]-normx); */
/*       dx = sqrt(l1*l2); */
/*       c1y = y[j] + cblas_ddot(dimension-1, xn, 1, y+j+1, 1); */
/*       c2y = 2*y[j] - c1y; */
/*       fx1 = (l1*c1y + dx*c2y)/2; */
/*       fx2 = (dx*c1y + l2*c2y)/2; */
/*       out[j] = fx1 + fx2 - dx*y[j]; */
/*       for (int k = 0; k < dimension-1; out[j+k+1] = fx1*xn[k] - fx2*xn[k] + dx*y[j+k+1], k++); */
/*     } */
/*   free(xn); */
/* } */

/* PA: Jordan algebra, returns inv(x) */
void Jinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, normx;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    normx = dnrm2l(dimension-1, x+j+1);
    l1 = 1/(x[j]+normx)/2;
    l2 = 1/(x[j]-normx)/2;
    out[j] = l1+l2;
    for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
  }
}

/* /\* PA: Return J_sqrt(x) *\/ */
/* void Jsqrt(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out) */
/* { */
/*   unsigned int dimension = (int)(vecSize / varsCount); */
/*   double l1, l2, normx; */
/*   size_t j; */

/*   for(size_t i = 0; i < varsCount; i++) */
/*     { */
/*       j = i*dimension; */
/*       normx = cblas_dnrm2(dimension-1, x+j+1, 1); */
/*       l1 = sqrt(x[j]+normx)/2; */
/*       l2 = sqrt(x[j]-normx)/2; */
/*       out[j] = l1+l2; */
/*       for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++); */
/*     } */
/* } */
/* /\* PA: Return J_sqrtinv(x) *\/ */
/* void Jsqrtinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out) */
/* { */
/*   unsigned int dimension = (int)(vecSize / varsCount); */
/*   double l1, l2, normx; */
/*   size_t j; */

/*   for(size_t i = 0; i < varsCount; i++) */
/*     { */
/*       j = i*dimension; */
/*       normx = cblas_dnrm2(dimension-1, x+j+1, 1); */
/*       l1 = 1/sqrt(x[j]+normx)/2; */
/*       l2 = 1/sqrt(x[j]-normx)/2; */
/*       out[j] = l1+l2; */
/*       for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++); */
/*     } */
/* } */

/* Returns J_sqrt(x) */
void Jsqrt(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, normx;

  for(size_t j = 0; j < vecSize; j+=dimension)
  {
    normx = dnrm2l(dimension-1, x+j+1);
    l1 = sqrtl(x[j]+normx)/2;
    l2 = sqrtl(x[j]-normx)/2;
    out[j] = l1+l2;
    for(int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
  }
}
/* Returns J_sqrtinv(x) */
void Jsqrtinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, normx;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    normx = dnrm2l(dimension-1, x+j+1);
    l1 = 1/sqrtl(x[j]+normx)/2;
    l2 = 1/sqrtl(x[j]-normx)/2;
    out[j] = l1+l2;
    for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
  }
}

/* PA: Return the Nesterov-Todd vector */
void Nesterov_Todd_vector(short T, const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * p)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);
  switch(T)
  {
  case 0:
  {
    Jsqrtinv(a, vecSize, varsCount, p); // NT-vector p
    break;
  }
  case 1:
  {
    Jsqrt(a, vecSize, varsCount, p);   // NT-vector p-inverse
    break;
  }
  case 2:
  {
    Jinv(a, vecSize, varsCount, p);   // NT-vector p-square
    break;
  }
  default:
  {
    printf("error in Nesterov_Todd_vector: T must be 0, 1 or 2\n");
    break;
  }
  }
  free(a);
  free(b);
}

/* PA: Return the Nesterov-Todd vector by means of the second formula */
void Nesterov_Todd_vector_b(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * p)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(y, x, vecSize, varsCount,a);
  Jsqrt(a, vecSize, varsCount, b);
  Qx50y(y, b, vecSize, varsCount, a);
  Jsqrtinv(a, vecSize, varsCount, p);

  free(a);
  free(b);
}

/* Computation of Qx*y by means of the formula 2*(x'*y)*x - det(x)*R*y */
void Qxy(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * z)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  size_t j;
  double xy;
  long double dx, nxb;

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    xy = 2*cblas_ddot(dimension, x+j, 1, y+j, 1);
    /* dx = x[j]*x[j]; */
    /* for (int k = 1; k < dimension; dx -= x[j+k]*x[j+k], k++); */
    nxb = dnrm2l(dimension-1, x+j+1);
    dx = (x[j] + nxb)*(x[j] - nxb);
    z[j] = xy*x[j] - dx*y[j];
    for (int k = 1; k < dimension; z[j+k] = xy*x[j+k] + dx*y[j+k], k++);
  }
}

/* Returns the product Q_{p}*z where p is the NT vector related to the pair (x,y) */
void QNTpz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);
  Qx50y(a, z, vecSize, varsCount, out);

  free(a);
  free(b);
}

/* Returns the product Q_{p^{-1}}*z where p is the NT vector related to the pair (x,y) */
void QNTpinvz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);
  Qx05y(a, z, vecSize, varsCount, out);

  free(a);
  free(b);
}

/* returns the Jordan product x^{-1} o y by using the formula x^{-1} = R*x/det(x), where R is the reflection matrix */
void Jxinvprody(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double nxb, detx, tmp;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;

    nxb = dnrm2l(dimension-1, x+j+1);
    detx = (x[j] + nxb) * (x[j] - nxb);

    tmp = x[j]*y[j];
    for (int k = 1; k < dimension; tmp -= x[j+k]*y[j+k], k++);
    out[j] = tmp/detx;

    for (int k = 1; k < dimension; out[j+k] = (x[j]*y[j+k] - y[j]*x[j+k])/detx, k++);
  }
}

/* returns the quadratic representation of a vector vec */
NumericsMatrix* QRmat(const double* const vec, const unsigned int vecSize, const size_t varsCount)
{
  size_t dimension = (size_t)(vecSize / varsCount);
  NumericsMatrix* out = NM_create(NM_SPARSE, vecSize, vecSize);
  NM_triplet_alloc(out, dimension * dimension * varsCount);
  //NM_fill(out, NM_SPARSE, vecSize, vecSize, out->matrix2);

  NumericsMatrix* quad_tmp = NM_create(NM_DENSE, dimension, dimension);

  long double nvec, nvecb, dvec;

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    for(unsigned int j = 0; j < dimension; ++j)
    {
      for(unsigned int k = 0; k < dimension; ++k)
        quad_tmp->matrix0[j+k*dimension] = 2.0*vec[i+j]*vec[i+k];
    }
    nvec = dnrm2l(dimension, vec+i);
    nvecb = dnrm2l(dimension-1, vec+i+1);
    dvec = (vec[i] + nvecb) * (vec[i] - nvecb);
    quad_tmp->matrix0[0] = nvec*nvec;
    for(unsigned int j = 1; j < dimension; ++j)
      quad_tmp->matrix0[j*(1+dimension)] += dvec;
    NM_insert(out, quad_tmp, i, i);
  }
  NM_clear(quad_tmp);
  free(quad_tmp);
  return out;
}

/* Returns a long double as the square root of determinant of a vector related to the Jordan product */
long double gammal(const double * const x, const size_t dimension)
{
  long double nxb, detx;
  nxb = dnrm2l(dimension-1, x+1);
  detx = (x[0] + nxb) * (x[0] - nxb);
  return(sqrtl(detx));
}

/*
   Returns the NT matrix by performing computation as in
   Solving semidefinite-linear programs using SDPT3
   by Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196
*/
NumericsMatrix* NTmat(const double* const x, const double* const z, const unsigned int vecSize, const size_t varsCount)
{
  size_t dimension = (size_t)(vecSize / varsCount);
  NumericsMatrix* out = NM_create(NM_SPARSE, vecSize, vecSize);
  NM_triplet_alloc(out, dimension * dimension * varsCount);
  //NM_fill(out, NM_SPARSE, vecSize, vecSize, out->matrix2);

  NumericsMatrix* G = NM_create(NM_DENSE, dimension, dimension);

  long double nvec, nvecb, dvec;

  long double gamx, gamz, w, gamt;
  double * t = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    gamx = gammal(x+i, dimension);
    gamz = gammal(z+i, dimension);
    w = sqrtl(gamz/gamx);
    t[0] = z[i]/w + w*x[i];
    for(unsigned int j = 1; j < dimension; ++j)
    {
      t[j] = z[i+j]/w - w*x[i+j];
    }
    gamt = gammal(t, dimension);
    for(unsigned int j = 0; j < dimension; ++j)
    {
      t[j] = t[j]/gamt;
    }
    G->matrix0[0] = t[0]*w;
    for(unsigned int j = 1; j < dimension; ++j)
    {
      G->matrix0[j] = t[j]*w;
      G->matrix0[j*dimension] = G->matrix0[j];
    }
    for(unsigned int j = 1; j < dimension; ++j)
    {
      for(unsigned int k = 1; k < dimension; ++k)
  G->matrix0[j+dimension*k] = (j==k) ? (1+t[j]*t[k]/(1+t[0]))*w : t[j]*t[k]/(1+t[0])*w;
    }
    NM_insert(out, G, i, i);
  }
  NM_clear(G);
  free(G);
  return out;
}

/*
   Returns the inverse of NT matrix by performing computation as in
   Solving semidefinite-linear programs using SDPT3
   by Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196
*/
NumericsMatrix* NTmatinv(const double* const x, const double* const z, const unsigned int vecSize, const size_t varsCount)
{
  size_t dimension = (size_t)(vecSize / varsCount);
  NumericsMatrix* out = NM_create(NM_SPARSE, vecSize, vecSize);
  NM_triplet_alloc(out, dimension * dimension * varsCount);
  //NM_fill(out, NM_SPARSE, vecSize, vecSize, out->matrix2);

  NumericsMatrix* G = NM_create(NM_DENSE, dimension, dimension);

  long double nvec, nvecb, dvec;

  long double gamx, gamz, w, gamt;
  double * t = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    gamx = gammal(x+i, dimension);
    gamz = gammal(z+i, dimension);
    w = sqrtl(gamz/gamx);
    t[0] = z[i]/w + w*x[i];
    for(unsigned int j = 1; j < dimension; ++j)
    {
      t[j] = z[i+j]/w - w*x[i+j];
    }
    gamt = gammal(t, dimension);
    for(unsigned int j = 0; j < dimension; ++j)
    {
      t[j] = t[j]/gamt;
    }
    G->matrix0[0] = t[0]/w;
    for(unsigned int j = 1; j < dimension; ++j)
    {
      G->matrix0[j] = -t[j]/w;
      G->matrix0[j*dimension] = G->matrix0[j];
    }
    for(unsigned int j = 1; j < dimension; ++j)
    {
      for(unsigned int k = 1; k < dimension; ++k)
  G->matrix0[j+dimension*k] = (j==k) ? (1+t[j]*t[k]/(1+t[0]))/w : t[j]*t[k]/(1+t[0])/w;
    }
    NM_insert(out, G, i, i);
  }
  NM_clear(G);
  free(G);
  return out;
}

/*
   Returns the square of NT matrix by performing computation as in
   Solving semidefinite-linear programs using SDPT3
   by Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196
*/
NumericsMatrix* NTmatsqr(const double* const x, const double* const z, const unsigned int vecSize, const size_t varsCount)
{
  size_t dimension = (size_t)(vecSize / varsCount);
  NumericsMatrix* out = NM_create(NM_SPARSE, vecSize, vecSize);
  NM_triplet_alloc(out, dimension * dimension * varsCount);
  //NM_fill(out, NM_SPARSE, vecSize, vecSize, out->matrix2);

  NumericsMatrix* G = NM_create(NM_DENSE, dimension, dimension);

  long double nvec, nvecb, dvec;

  long double gamx, gamz, w, w2, gamt;
  double * t = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    gamx = gammal(x+i, dimension);
    gamz = gammal(z+i, dimension);
    w2 = gamz/gamx;
    w = sqrtl(w2);
    t[0] = z[i]/w + w*x[i];
    for(unsigned int j = 1; j < dimension; ++j)
    {
      t[j] = z[i+j]/w - w*x[i+j];
    }
    gamt = gammal(t, dimension);
    for(unsigned int j = 0; j < dimension; ++j)
    {
      t[j] = t[j]/gamt;
    }
    G->matrix0[0] = w2*cblas_ddot(dimension, t, 1, t, 1);
    for(unsigned int j = 1; j < dimension; ++j)
    {
      G->matrix0[j] = 2*t[0]*t[j]*w2;
      G->matrix0[j*dimension] = G->matrix0[j];
    }
    for(unsigned int j = 1; j < dimension; ++j)
    {
      for(unsigned int k = 1; k < dimension; ++k)
  G->matrix0[j+dimension*k] = (j==k) ? (1+2*t[j]*t[k])*w2 : 2*t[j]*t[k]*w2;
    }
    NM_insert(out, G, i, i);
  }
  NM_clear(G);
  free(G);
  return out;
}




