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

#include "JordanAlgebra.h"
#include "cblas.h"


RawNumericsMatrix* Arrow_repr(const unsigned int vecSize, const double* const vec, const size_t varsCount)
{
    /* validation */
    if (vecSize % varsCount != 0)
    {
        fprintf(stderr, "Arrow_repr: %d variables can not be extracted from vector of size %d.\n", varsCount, vecSize);
        exit(EXIT_FAILURE);
    }

    size_t dimension = (size_t)(vecSize / varsCount);
    if (dimension < 2)
    {
        fprintf(stderr, "Arrow_repr: The dimension of variables can not be less than 2 but given %d.\n", dimension);
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
        NM_zentry(Arw_mat, pos, pos, vec[pos]);

        for(size_t j = 1; j < dimension; ++j)
        {
            NM_zentry(Arw_mat, pos, pos + j, vec[pos + j]);
            NM_zentry(Arw_mat, pos + j, pos, vec[pos + j]);
            NM_zentry(Arw_mat, pos + j, pos + j, vec[pos]);
        }
    }
    return Arw_mat;
}


/* Jordan product of two vectors */
double* JA_prod(const double * const x, const double * const y, const unsigned int vecSize, const int varsCount)
{
    assert(x);
    assert(y);
    double * out = (double *)calloc(vecSize, sizeof(double));
    unsigned int dimension = (int)(vecSize / varsCount);

    unsigned int pos;
    for(size_t i = 0; i < varsCount; ++i)
    {
        pos = i * dimension;
        out[pos] = cblas_ddot(dimension, x + pos, 1, y + pos, 1);
        for(size_t j = 1; j < dimension; ++j)
            out[pos + j] = x[pos] * y[pos + j] + y[pos] * x[pos + j];
    }
    return out;
}
