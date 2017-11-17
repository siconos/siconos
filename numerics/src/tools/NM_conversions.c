/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <stdio.h>

#include "NM_conversions.h"
#include "SiconosConfig.h"
#include "csparse.h"

#ifdef WITH_MKL_SPBLAS
#include "tlsdef.h"
#include "MKL_common.h"
typedef void (*mkl_dcsrcoo_t) (const __INT_T *job , const __INT_T *n , double *acsr , __INT_T *ja , __INT_T *ia , __INT_T *nnz , double *acoo , __INT_T *rowind , __INT_T *colind , __INT_T *info );
typedef void (*mkl_dcsrcsc_t) (const __INT_T *job , const __INT_T *n , double *acsr , __INT_T *ja , __INT_T *ia , double *acsc , __INT_T *ja1 , __INT_T *ia1 , __INT_T *info );

tlsvar mkl_dcsrcoo_t mkl_dcsrcoo_p = NULL;
tlsvar mkl_dcsrcsc_t mkl_dcsrcsc_p = NULL;

#endif

CSparseMatrix* NM_csc_to_triplet(CSparseMatrix* csc)
{
  assert(csc);
  CSparseMatrix* triplet = cs_spalloc(csc->m, csc->n, csc->p[csc->n], 1, 1);

  csi* Ap = csc->p;
  csi* Ai = csc->i;
  double* val = csc->x;
  for (csi j = 0; j < csc->n; ++j)
  {
    for (csi i = Ap[j]; i < Ap[j+1]; ++i)
    {
      cs_zentry(triplet, Ai[i], j, val[i]);
    }
  }
  return triplet;
}

CSparseMatrix* NM_triplet_to_csr(CSparseMatrix* triplet)
{
#ifdef WITH_MKL_SPBLAS
  assert(triplet);
  CHECK_MKL(load_mkl_function("mkl_dcsrcoo", (void**)&mkl_dcsrcoo_p));
  if (triplet->m != triplet->n)
  {
    fprintf(stderr, "NM_triplet_to_csr :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  CSparseMatrix* csr = cs_spalloc(NS_NROW_CSR(triplet), NS_NCOL_CSR(triplet), triplet->nz, 1, 0);
  assert(csr);
  csr->nz = -2;

  csi n = csr->n;
  csi job[6] = {0};
  csi info = 0;
  job[0] = 2;
  (*mkl_dcsrcoo_p)(job, &n, csr->x, csr->i, csr->p, &(triplet->nz), triplet->x, triplet->i, triplet->p, &info);

  return csr;
#else
  fprintf(stderr, "NM_triplet_to_csr :: MKL not enabled\n");
  exit(EXIT_FAILURE);
#endif
}

CSparseMatrix* NM_csr_to_triplet(CSparseMatrix* csr)
{
#ifdef WITH_MKL_SPBLAS
  assert(csr);
  CHECK_MKL(load_mkl_function("mkl_dcsrcoo", (void**)&mkl_dcsrcoo_p));
  if (csr->m != csr->n)
  {
    fprintf(stderr, "NM_csr_to_triplet :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  CSparseMatrix* triplet = cs_spalloc(csr->m, csr->n, csr->p[csr->m], 1, 1);
  assert(triplet);

  csi n = csr->n;
  csi job[6] = {0};
  job[4] = csr->p[csr->m];
  job[5] = 3;
  csi info = 0;
  (*mkl_dcsrcoo_p)(job, &n, csr->x, csr->i, csr->p, &(csr->p[csr->m]), triplet->x, triplet->i, triplet->p, &info);
  triplet->nz = csr->p[csr->m];

  return triplet;
#else
  fprintf(stderr, "NM_csr_to_triplet :: MKL not enabled\n");
  exit(EXIT_FAILURE);
#endif
}

CSparseMatrix* NM_csc_to_csr(CSparseMatrix* csc)
{
#ifdef WITH_MKL_SPBLAS
  assert(csc);
  CHECK_MKL(load_mkl_function("mkl_dcsrcsc", (void**)&mkl_dcsrcsc_p));
  if (csc->m != csc->n)
  {
    fprintf(stderr, "NM_csc_to_csr :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  CSparseMatrix* csr = cs_spalloc(NS_NROW_CSR(csc), NS_NCOL_CSR(csc), csc->p[csc->n], 1, 0);
  assert(csr);
  csr->nz = -2;

  csi n = csr->n;
  csi job[6] = {0};
  csi info = 0;
  job[0] = 1;
  job[5] = 1;
  (*mkl_dcsrcsc_p)(job, &n, csr->x, csr->i, csr->p, csc->x, csc->i, csc->p, &info);

  return csr;
#else
  fprintf(stderr, "NM_csc_to_csr :: MKL not enabled\n");
  exit(EXIT_FAILURE);
#endif
}

CSparseMatrix* NM_csr_to_csc(CSparseMatrix* csr)
{
#ifdef WITH_MKL_SPBLAS
  assert(csr);
  CHECK_MKL(load_mkl_function("mkl_dcsrcsc", (void**)&mkl_dcsrcsc_p));
  if (csr->m != csr->n)
  {
    fprintf(stderr, "NM_csr_to_csc :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  assert(csr);
  CSparseMatrix* csc = cs_spalloc(csr->m, csr->n, csr->nzmax, 1, 0);
  assert(csc);

  csi n = csr->n;
  csi job[6] = {0};
  job[5] = 1;
  csi info = 0;
  (*mkl_dcsrcsc_p)(job, &n, csr->x, csr->i, csr->p, csc->x, csc->i, csc->p, &info);

  return csc;
#else
  fprintf(stderr, "NM_csr_to_csc :: MKL not enabled\n");
  exit(EXIT_FAILURE);
#endif
}
