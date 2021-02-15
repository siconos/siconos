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
#include "NumericsVector.h"
#include <math.h>    // for fabs
#include <stdio.h>   // for fprintf, printf, FILE, stderr
#include <stdlib.h>  // for exit, EXIT_FAILURE
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"   // for DEBUG_PRINTF
#include "SiconosBlas.h"
#include "NumericsMatrix.h"
#include "float.h"

void NV_display(const double * const m, int nRow)
{
  int lin;
  printf("vector of size\t%d\t =\n[", nRow);
  if(nRow == 0)
  {
    printf("]\n");
  }
  for(lin = 0; lin < nRow; lin++)
  {
    printf(" %.15e", m[lin]);
    if(lin != nRow - 1)
      printf(", ");
    else
      printf("]\n");
  }

}

void NV_copy(const double * const vec, unsigned int vecSize, double * out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = vec[i];
}

void NV_write_in_file_python(double * m,  int nRow, FILE* file)
{
  if(! m)
  {
    fprintf(stderr, "Numerics, NV_write_in_file_python  failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "size = %d; \n", nRow);
  fprintf(file, "data= [");
  for(int i = 0; i < nRow; i++)
  {
    fprintf(file, "%32.24e,\t ", m[i]);
  }
  fprintf(file, "]");
}

bool NV_equal(double * x, double * y, int n, double tol)
{
  for(int i =0; i< n ; i++)
  {
    if(fabs(x[i] - y[i]) >= tol)
    {
      DEBUG_PRINTF("error %i = %e\n",i, fabs(x[i]) - y[i]);
      return false;
    }
  }
  return true;
}

void NV_insert(double * x, const unsigned int xSize,
               const double * const y, const unsigned int ySize,
               unsigned int i)
{
  if(xSize < ySize)
  {
    fprintf(stderr, "NV_insert ::  the vector to be inserted is greater than the given vector: size_x < size_y - %d < %d\n", xSize, ySize);
    exit(EXIT_FAILURE);
  }
  if(i + ySize > xSize)
  {
    fprintf(stderr, "NV_insert ::  the vector to be inserted is too big for insertion from position %d\n", i);
    exit(EXIT_FAILURE);
  }
  for(unsigned int j = i; j < i + ySize; ++j)
    x[j] = y[j - i];
}

void NV_power2(const double * const vec, const unsigned int vecSize, double * out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = vec[i] * vec[i];
}


double NV_reduce(const double * const vec, const unsigned int vecSize)
{
  register double sum = 0.0;
  for(unsigned int i = 0; i < vecSize; ++i)
    sum += vec[i];
  return sum;
}

void NV_prod(const double * const vec1, const double * const vec2, const unsigned int vecSize, double * out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = vec1[i] * vec2[i];
}

double* NV_div(const double * const x, const double * const y, const unsigned int vecSize)
{
  double * out = (double*)malloc(vecSize * sizeof(double));
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = x[i] / (y[i] + 1e-12);
  return out;
}

double NV_min(const double * const vec, const unsigned int vecSize)
{
  double min_elem = DBL_MAX;
  for(unsigned int i = 0; i < vecSize; ++i)
    if(vec[i] < min_elem)
      min_elem = vec[i];
  return min_elem;
}

double NV_max(const double * const vec, const unsigned int vecSize)
{
  double max_elem = DBL_MIN;
  for(unsigned int i = 0; i < vecSize; ++i)
    if(vec[i] > max_elem)
      max_elem = vec[i];
  return max_elem;
}

double * NV_abs(const double * const vec, const unsigned int vecSize)
{
  double * out = (double*)malloc(vecSize * sizeof(double));
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = fabs(vec[i]);
  return out;
}

void NV_add(const double * const vec1, const double * const vec2, const unsigned int vecSize, double * out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = vec1[i] + vec2[i];
}

void NV_const_add(const double * const vec, const unsigned int vecSize, const double alpha, const double beta, double * out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = alpha * vec[i] + beta;
}

void NV_sub(const double * const vec1, const double * const vec2, const unsigned int vecSize, double * out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = vec1[i] - vec2[i];
}

double NV_norm_inf(const double * const vec, const unsigned int vecSize)
{
  /* double * abs_vec = NV_abs(vec, vecSize); */
  /* return NV_max(abs_vec, vecSize); */
  double norm = DBL_MIN;
  for(unsigned int i = 0; i < vecSize; ++i)
  {
    norm = fmax(norm, fabs(vec[i]));
  }
  return norm;
}

double NV_norm_2(const double * const vec, const unsigned int vecSize)
{
  /* double * vec2 = (double*)calloc(vecSize, sizeof(double)); */
  /* NV_power2(vec, vecSize, vec2); */
  /* double sum = NV_reduce(vec2, vecSize); */
  /* free(vec2); */
  /* return sqrt(sum); */
  double norm = cblas_dnrm2(vecSize, vec, 1);
  assert(!isnan(norm));
  return norm;
}

void NV_sqrt(const double * const vec, const unsigned int vecSize, double * out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    out[i] = sqrt(vec[i]);
}

void NV_dott(const double * const vec1, const double * const vec2, const unsigned int vecSize, NumericsMatrix* out)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    for(unsigned int j = 0; j < vecSize; ++j)
      NM_entry(out, i, j, vec1[i] * vec2[j]);
}

int NV_isnan(const double * const vec,  const unsigned int vecSize)
{
  for(unsigned int i = 0; i < vecSize; ++i)
    if(isnan(vec[i]))
      return 1;
  return 0;
}

