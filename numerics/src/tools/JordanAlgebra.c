/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

//#define DEBUG_MESSAGES
#include "debug.h"


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
  unsigned int nzmax = (dimension * 3 - 2) * varsCount;
  NM_triplet_alloc(Arw_mat, nzmax);
  NM_fill(Arw_mat, NM_SPARSE, vecSize, vecSize, Arw_mat->matrix2);

  /* Arrow matirx filling */
  unsigned int pos;
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

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    NV_dott(vec + i, vec + i, dimension, quad_tmp);

    for(unsigned int j = 0; j < dimension; ++j)
    {
      for(unsigned int k = 0; k < dimension; ++k)
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

  assert(!NV_isnan(out,vecSize));


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
    for(unsigned int j = 1; j < dimension; ++j)
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
    assert(!isnan(vec[pos]));
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
  const double EPS = 1e-12;
  unsigned int dimension = (int)(vecSize / varsCount);
  register unsigned int pos;
  double xi_bar_norm;

  for(size_t i = 0; i < 2*varsCount; i += 2)
  {
    pos = (i / 2.) * dimension;
    xi_bar_norm = NV_norm_2(vec + pos + 1, dimension - 1);
    out[i][0] = 0.5;
    NV_const_add(vec + pos + 1, dimension - 1, 1. / (2 * xi_bar_norm + EPS), 0, out[i] + 1);
  }

  for(size_t i = 1; i < 2*varsCount; i += 2)
  {
    pos = ((i - 1) / 2.) * dimension;
    xi_bar_norm = NV_norm_2(vec + pos + 1, dimension - 1);
    out[i][0] = 0.5;
    NV_const_add(vec + pos + 1, dimension - 1, - 1. / (2 * xi_bar_norm + EPS), 0, out[i] + 1);
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
    assert(eigenvals[i]>0);
    assert(eigenvals[i+1]>0);
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
