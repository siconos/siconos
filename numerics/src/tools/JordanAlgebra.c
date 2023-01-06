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

#include "JordanAlgebra.h"
#include "cblas.h"
#include "math.h"
#include "NumericsVector.h"
#include <float.h>               // for DBL_EPSILON, DBL_MAX

//#define DEBUG_MESSAGES
#include "siconos_debug.h"
typedef long double float_type;
/* typedef double float_type; */
#define EPS 1e-40

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
    // out[i] = vec[pos] + NV_norm_2(vec + pos + 1, dimension - 1);
    out[i] = vec[pos] + dnrm2l(dimension - 1, vec + pos + 1);
  }

  for(size_t i = 1; i < 2*varsCount; i += 2)
  {
    pos = ((i - 1) / 2.) * dimension;
    // out[i] = vec[pos] - NV_norm_2(vec + pos + 1, dimension - 1);
    out[i] = vec[pos] - dnrm2l(dimension - 1, vec + pos + 1);
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
  // double xi_bar_norm;
  float_type xi_bar_norm;

  for(size_t i = 0; i < 2*varsCount; i += 2)
  {
    pos = (i / 2.) * dimension;
    // xi_bar_norm = NV_norm_2(vec + pos + 1, dimension - 1);
    xi_bar_norm = dnrm2l(dimension - 1, vec + pos + 1);
    out[i][0] = 0.5;
    //NV_const_add(vec + pos + 1, dimension - 1, 1. / (2 * xi_bar_norm + EPS), 0, out[i] + 1);
    NV_const_add(vec + pos + 1, dimension - 1, 1. / (2 * xi_bar_norm), 0, out[i] + 1);
  }

  for(size_t i = 1; i < 2*varsCount; i += 2)
  {
    pos = ((i - 1) / 2.) * dimension;
    // xi_bar_norm = NV_norm_2(vec + pos + 1, dimension - 1);
    xi_bar_norm = dnrm2l(dimension - 1, vec + pos + 1);
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
    // if (eigenvals[i] == 0.) sqrt_eigenval1 = 1. / 1e-20; else sqrt_eigenval1 = 1. / eigenvals[i];
    // if (eigenvals[i+1] == 0.) sqrt_eigenval2 = 1. / 1e-20; else sqrt_eigenval2 = 1. / eigenvals[i+1];
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
float_type dnrm2l(const unsigned int n, const double * x)
{
  float_type norm, scale, ssq, absxi, quo;

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
      if(ssq < 0.) printf("\n dnrm2l. NaN\n");
      norm = scale * sqrtl(ssq);
    }
  }
  return norm;
}


/* Returns the square of 2-norm of a vector - uses long double - based on blas_dnrm2 */
float_type dnrm2sqrl(const unsigned int n, const double * x)
{
  float_type norm, scale, ssq, absxi, quo;

  if (n < 1)
    norm = 0.0;
  else if (n == 1)
    norm = fabsl(x[0]);
  else
  {
    scale = 0.0;
    ssq = 1.0;
    for (size_t i = 0; i < n; i++)
    {
      if (x[i] != 0)
      {
        absxi = fabs(x[i]);
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
      norm = scale * scale * ssq;
    }
  }
  return norm;
}

/* Returns the product Q_sqrt(x)*y */
void Qx05y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  float_type l1, l2, c1y, c2y, nxb, fx1, fx2, dx;
  size_t j;
  float_type *xb = (float_type*)calloc(dimension-1, sizeof(float_type));

  for (int i = 0; i < dimension - 1; xb[i] = 1/sqrtl(dimension-1), i++);

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    nxb = dnrm2l(dimension-1, x+j+1);
    if (nxb > 0)
      for (int k = 0; k < dimension-1; xb[k] = x[j+1+k]/nxb, k++);
    l1 = x[j]+nxb;
    l2 = x[j]-nxb;
    // if (l2 <= 0.) {printf("Qx05y. i = %zu: l2 = %3.50Le => l2 = EPS\n", i, l2); l2 = EPS;} // to avoid negative number b/c of different data types
    // if (l2 <= 0.) l2 = fabsl(l2); // to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Qx05y. i = %zu: l2 = %3.50Le => l2 = fabsl(l2)\n", i, l2); l2 = fabsl(l2);} // to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Qx05y. i = %zu: l2 = %3.50Le => l2 = fabsl(l2)/10\n", i, l2); l2 = fabsl(l2)/10.;} // to avoid negative number b/c of different data types
    dx = sqrtl(l1*l2);
    c1y = y[j];
    for (int k = 0; k < dimension-1; c1y += xb[k]*y[j+1+k], k++);
    c2y = 2*y[j] - c1y;
    fx1 = (l1*c1y + dx*c2y)/2;
    fx2 = (dx*c1y + l2*c2y)/2;
    out[j] = fx1 + fx2 - dx*y[j];
    //    out[j] = cblas_ddot(dimension, x+j, 1, y+j, 1);
    for (int k = 0; k < dimension-1; out[j+k+1] = fx1*xb[k] - fx2*xb[k] + dx*y[j+k+1], k++);
  }
  free(xb);
}

/* Returns the product Q_inv_sqrt(x)*y */
void Qx50y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  float_type l1, l2, c1y, c2y, nxb, fx1, fx2, dx;
  size_t j;
  float_type *xb = (float_type*)calloc(dimension-1, sizeof(float_type));

  for (int i = 0; i < dimension - 1; xb[i] = 1/sqrtl(dimension-1), i++); // useful ?

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    nxb = dnrm2l(dimension-1, x+j+1);
    /* if (isnan(nxb)) */
    /*   { */
    /* 	printf("%Le\n",nxb); */
    /* 	//getchar(); */
    /*   } */

    if (nxb > 0)
      for (int k = 0; k < dimension-1; xb[k] = x[j+1+k]/nxb, k++);

    l1 = x[j]+nxb;
    l2 = x[j]-nxb;
    //if (l2 <= 0.) {printf("Qx50y. i = %zu: l2 = %3.50Le => l2 = EPS\n", i, l2); l2 = EPS;} // to avoid negative number b/c of different data types
    // if (l2 <= 0.) l2 = fabsl(l2); // to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Qx50y. i = %zu: l2 = %3.50Le => l2 = fabsl(l2)\n", i, l2); l2 = fabsl(l2);} // to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Qx50y. i = %zu: l2 = %3.50Le => l2 = fabsl(l2)/10\n", i, l2); l2 = fabsl(l2)/10.;} // to avoid negative number b/c of different data types
    dx = 1/sqrtl(l1*l2);
    /* if (isnan(dx)) */
    /*   { */
    /* 	printf("%Le %Le\n",l1,l2); */
    /* 	getchar(); */
    /*   } */
	      
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

/* PA: Jordan algebra, returns inv(x) */
void Jinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  float_type l1, l2, normx;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    normx = dnrm2l(dimension-1, x+j+1);
    l1 = 1/(x[j]+normx)/2;
    l2 = 1/(x[j]-normx)/2;
    // l2 = x[j]-normx;
    // if (l2 == 0.) l2 = 1./1e-20; else l2 = 1/l2/2;// to avoid divide by 0 b/c of different data types
    // if (l2 == 0.) {printf("Jinv. i = %zu: l2 = %3.50Le => l2 = 1./1e-20\n", i, l2); l2 = 1./1e-20;} else l2 = 1/(x[j]-normx)/2;// to avoid divide by 0 b/c of different data types
    out[j] = l1+l2;
    for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
  }
}

/* Returns J_sqrt(x) */
void Jsqrt(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  float_type l1, l2, normx;

  for(size_t j = 0; j < vecSize; j+=dimension)
  {
    normx = dnrm2l(dimension-1, x+j+1);
    l1 = sqrtl(x[j]+normx)/2;
    l2 = sqrtl(x[j]-normx)/2;
    // l2 = x[j]-normx;
    // if (l2 <= 0.) {printf("Jsqrt. i = %zu: l2 = %3.50Le => l2 = EPS\n", j, l2); l2 = EPS;} else l2 = sqrtl(x[j]-normx)/2;// to avoid negative number b/c of different data types
    // if (l2 <= 0.) l2 = sqrtl(fabsl(l2))/2; else l2 = sqrtl(l2)/2;// to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Jsqrt. i = %zu: l2 = %3.50Le => l2 = sqrtl(fabsl(x[j]-normx))/2\n", j, l2); l2 = sqrtl(fabsl(x[j]-normx))/2;} else l2 = sqrtl(x[j]-normx)/2;// to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Jsqrt. i = %zu: l2 = %3.50Le => l2 = sqrtl(fabsl(x[j]-normx)/10)/2\n", j, l2); l2 = sqrtl(fabsl(x[j]-normx)/10.)/2;} else l2 = sqrtl(x[j]-normx)/2;// to avoid negative number b/c of different data types
    out[j] = l1+l2;
    for(int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
  }
}
/* Returns J_sqrtinv(x) */
void Jsqrtinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  float_type l1, l2, normx;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;
    normx = dnrm2l(dimension-1, x+j+1);
    l1 = 1/sqrtl(x[j]+normx)/2;
    l2 = 1/sqrtl(x[j]-normx)/2;
    // l2 = x[j]-normx;
    // if (l2 <= 0.) {printf("Jsqrtinv. i = %zu: l2 = %3.50Le => l2 = 1./sqrtl(EPS)\n", i, l2); l2 = 1./sqrtl(EPS);} else l2 = 1./sqrtl(x[j]-normx)/2;// to avoid negative number b/c of different data types
    // if (l2 <= 0.) l2 = 1./sqrtl(fabsl(l2))/2; else l2 = 1./sqrtl(l2)/2;// to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Jsqrtinv. i = %zu: l2 = %3.50Le => l2 = 1./sqrtl(fabsl(x[j]-normx))/2\n", i, l2); l2 = 1./sqrtl(fabsl(x[j]-normx))/2;} else l2 = 1./sqrtl(x[j]-normx)/2;// to avoid negative number b/c of different data types
    // if (l2 <= 0.) {printf("Jsqrtinv. i = %zu: l2 = %3.50Le => l2 = 1./sqrtl(fabsl(x[j]-normx)/10)/2\n", i, l2); l2 = 1./sqrtl(fabsl(x[j]-normx)/10.)/2;} else l2 = 1./sqrtl(x[j]-normx)/2;// to avoid negative number b/c of different data types

    /* if (x[j]-normx<1e-14) */
    /*   { */
    /* 	printf("Jsqrtinv: %e %Le %Le\n",x[j], x[j]+normx, x[j]-normx); */
    /* 	getchar(); */
    /*   } */
    out[j] = l1+l2;
    for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
  }
}

/* Computation of the Nesterov-Todd vector 
   T = 0 -> p 
   T = 1 -> p^{-1}
   T = 2 -> p^2
   T = 3 -> p^{-2}
*/
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
  case 3:
  {
    NV_copy(a, vecSize, p); 
    break;
  }
  default:
  {
    printf("error in Nesterov_Todd_vector: T must be 0, 1, 2 or 3\n");
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
  float_type dx, nxb;

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

  /* alternative */
  /* Qx05y(y, x, vecSize, varsCount,a); */
  /* Jsqrt(a, vecSize, varsCount, b); */
  /* Qx50y(y, b, vecSize, varsCount, a); */
  /* Qx05y(a, z, vecSize, varsCount, out); */

  free(a);
  free(b);
}

/* Returns the product Q_{p^{-2}}*z where p is the NT vector related to the pair (x,y) */
void QNTpinv2z(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);
  Qxy(a, z, vecSize, varsCount, out);

  free(a);
  free(b);
}

/* returns the Jordan product x^{-1} o y by using the formula x^{-1} = R*x/det(x), where R is the reflection matrix */
void Jxinvprody(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  float_type nxb, detx, tmp;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
  {
    j = i*dimension;

    nxb = dnrm2l(dimension-1, x+j+1);
    detx = (x[j] + nxb) * (x[j] - nxb);
    // detx = x[j]-nxb;
    // if (detx <= 0.) {printf("Jxinvprody. i = %zu: l2 = %3.50Le => l2 = EPS\n", i, detx); detx = (x[j] + nxb) * EPS;}
    // if (detx <= 0.) detx = fabsl(detx);
    // if (detx <= 0.) {printf("Jxinvprody. i = %zu: l2 = %3.50Le => det = (x[j] + nxb) * fabsl(x[j]-nxb)\n", i, detx); detx = (x[j] + nxb) * fabsl(x[j]-nxb);}
    // if (detx <= 0.) {printf("Jxinvprody. i = %zu: l2 = %3.50Le => det = (x[j] + nxb) * fabsl(x[j]-nxb)/10\n", i, detx); detx = (x[j] + nxb) * fabsl(x[j]-nxb)/10.;}
    // else detx = (x[j] + nxb) * (x[j] - nxb);// to avoid divide by 0 b/c of different data types

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

  float_type nvec, nvecb, dvec;

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
float_type ld_gammal(const double * const x, const size_t dimension)
{
  float_type nxb, detx;
  nxb = dnrm2l(dimension-1, x+1);
  detx = (x[0] + nxb) * (x[0] - nxb);
  // detx = x[0]-nxb;
  // if (detx <= 0.) {printf("ld_gammal. l2 = %3.50Le => l2 = EPS\n", detx); detx = (x[0] + nxb) * EPS;}
  // if (detx <= 0.) detx = fabsl(detx);
  // if (detx <= 0.) {printf("ld_gammal. l2 = %3.50Le => detx = (x[0] + nxb) * fabsl(x[0]-nxb)\n", detx); detx = (x[0] + nxb) * fabsl(x[0]-nxb);}
  // if (detx <= 0.) {printf("ld_gammal. l2 = %3.50Le => detx = (x[0] + nxb) * fabsl(x[0]-nxb)/10\n", detx); detx = (x[0] + nxb) * fabsl(x[0]-nxb)/10.;}
  // else detx = (x[0] + nxb) * (x[0] - nxb);// to avoid negative number b/c of different data types

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

  float_type nvec, nvecb, dvec;

  float_type gamx, gamz, w, gamt;
  double * t = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    gamx = ld_gammal(x+i, dimension);
    gamz = ld_gammal(z+i, dimension);
    w = sqrtl(gamz/gamx);
    t[0] = z[i]/w + w*x[i];
    for(unsigned int j = 1; j < dimension; ++j)
    {
      t[j] = z[i+j]/w - w*x[i+j];
    }
    gamt = ld_gammal(t, dimension);
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

  float_type nvec, nvecb, dvec;

  float_type gamx, gamz, w, gamt;
  double * t = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    gamx = ld_gammal(x+i, dimension);
    gamz = ld_gammal(z+i, dimension);
    w = sqrtl(gamz/gamx);
    t[0] = z[i]/w + w*x[i];
    for(unsigned int j = 1; j < dimension; ++j)
    {
      t[j] = z[i+j]/w - w*x[i+j];
    }
    gamt = ld_gammal(t, dimension);
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

  float_type nvec, nvecb, dvec;

  float_type gamx, gamz, w, w2, gamt;
  double * t = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    gamx = ld_gammal(x+i, dimension);
    gamz = ld_gammal(z+i, dimension);
    w2 = gamz/gamx;
    w = sqrtl(w2);
    t[0] = z[i]/w + w*x[i];
    for(unsigned int j = 1; j < dimension; ++j)
    {
      t[j] = z[i+j]/w - w*x[i+j];
    }
    gamt = ld_gammal(t, dimension);
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
