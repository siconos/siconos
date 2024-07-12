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



#include <assert.h>         // for assert
#include <stdio.h>          // for fprintf, stderr
#include <stdlib.h>         // for exit, EXIT_FAILURE
#include "CSparseMatrix_internal.h"  // for CSparseMatrix, CS_INT
#include "NM_conversions.h"
#include "SiconosConfig.h"  // for WITH_MKL_SPBLAS  // IWYU pragma: keep

#ifdef WITH_MKL_SPBLAS
#include "tlsdef.h"
#include "MKL_common.h"
typedef void (*mkl_dcsrcoo_t)(const __INT_T *job, const __INT_T *n, double *acsr, __INT_T *ja, __INT_T *ia, __INT_T *nnz, double *acoo, __INT_T *rowind, __INT_T *colind, __INT_T *info);
typedef void (*mkl_dcsrcsc_t)(const __INT_T *job, const __INT_T *n, double *acsr, __INT_T *ja, __INT_T *ia, double *acsc, __INT_T *ja1, __INT_T *ia1, __INT_T *info);

tlsvar mkl_dcsrcoo_t mkl_dcsrcoo_p = NULL;
tlsvar mkl_dcsrcsc_t mkl_dcsrcsc_p = NULL;

#endif

CSparseMatrix* NM_csc_to_triplet(CSparseMatrix* csc)
{
  assert(csc);
  CSparseMatrix* triplet = cs_spalloc(csc->m, csc->n, csc->p[csc->n], 1, 1);

  CS_INT* Ap = csc->p;
  CS_INT* Ai = csc->i;
  double* val = csc->x;
  for(CS_INT j = 0; j < csc->n; ++j)
  {
    for(CS_INT i = Ap[j]; i < Ap[j+1]; ++i)
    {
      CSparseMatrix_entry(triplet, Ai[i], j, val[i]);
    }
  }
  return triplet;
}

CSparseMatrix* NM_csc_to_half_triplet(CSparseMatrix* csc)
{
  assert(csc);
  CSparseMatrix* triplet = cs_spalloc(csc->m, csc->n, csc->p[csc->n], 1, 1);
  if(!triplet) return (cs_done(triplet, NULL, NULL, 0)) ;
  CS_INT* Ap = csc->p;
  CS_INT* Ai = csc->i;
  double* val = csc->x;
  for(CS_INT j = 0; j < csc->n; ++j)
  {
    for(CS_INT i = Ap[j]; i < Ap[j+1]; ++i)
    {
      CSparseMatrix_symmetric_entry(triplet, Ai[i], j, val[i]);
    }
  }
  return triplet;
}

CSparseMatrix* NM_triplet_to_csr(CSparseMatrix* triplet)
{
#ifdef WITH_MKL_SPBLAS
  assert(triplet);
  CHECK_MKL(load_mkl_function("mkl_dcsrcoo", (void**)&mkl_dcsrcoo_p));
  if(triplet->m != triplet->n)
  {
    fprintf(stderr, "NM_triplet_to_csr :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  CSparseMatrix* csr = cs_spalloc(NSM_NROW_CSR(triplet), NSM_NCOL_CSR(triplet), triplet->nz, 1, 0);
  assert(csr);
  csr->nz = -2;

  CS_INT n = csr->n;
  CS_INT job[6] = {0};
  CS_INT info = 0;
  job[0] = 2;
  (*mkl_dcsrcoo_p)(job, &n, csr->x, csr->i, csr->p, &(triplet->nz), triplet->x, triplet->i, triplet->p, &info);

  return csr;
#else

    CS_INT m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    CS_ENTRY *Cx, *Tx ;
    if (!CS_TRIPLET (triplet)) return (NULL) ;                /* check inputs */

    m = triplet->m ; n = triplet->n ; Ti = triplet->i ; Tj = triplet->p ; Tx = triplet->x ; nz = triplet->nz ;

    CSparseMatrix* csr = cs_spalloc (m, n, nz, Tx != NULL, 0) ;          /* allocate result */

    w = cs_calloc (n, sizeof (CS_INT)) ;

    if (!csr || !w) return (cs_done (csr, w, NULL, 0)) ;    /* out of memory */
    for (k = 0 ; k < nz ; k++) w [Ti [k]]++ ;           /* rows counts */

    Cp = csr->p ; Ci = csr->i ; Cx = csr->x ;

    cs_cumsum (Cp, w, n) ;    /* row pointers */
    for (k = 0 ; k < nz ; k++)
    {
        Ci [p = w [Ti [k]]++] = Tj [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
    }

    csr->nz=-2;

    return (cs_done (csr, w, NULL, 1)) ;      /* success; free w and return C */


  fprintf(stderr, "NM_triplet_to_csr :: MKL not enabled\n");
  exit(EXIT_FAILURE);
#endif
}

CSparseMatrix* NM_csr_to_triplet(CSparseMatrix* csr)
{
#ifdef WITH_MKL_SPBLAS
  assert(csr);
  CHECK_MKL(load_mkl_function("mkl_dcsrcoo", (void**)&mkl_dcsrcoo_p));
  if(csr->m != csr->n)
  {
    fprintf(stderr, "NM_csr_to_triplet :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  CSparseMatrix* triplet = cs_spalloc(csr->m, csr->n, csr->p[csr->m], 1, 1);
  assert(triplet);

  CS_INT n = csr->n;
  CS_INT job[6] = {0};
  job[4] = csr->p[csr->m];
  job[5] = 3;
  CS_INT info = 0;
  (*mkl_dcsrcoo_p)(job, &n, csr->x, csr->i, csr->p, &(csr->p[csr->m]), triplet->x, triplet->i, triplet->p, &info);
  triplet->nz = csr->p[csr->m];

  return triplet;
#else

  // Ugly
  CSparseMatrix* csc = NM_csr_to_csc(csr);
  CSparseMatrix* triplet= NM_csc_to_triplet(csc);

  cs_spfree(csc);
  return triplet;

  fprintf(stderr, "NM_csr_to_triplet :: MKL not enabled\n");
  exit(EXIT_FAILURE);
#endif
}

CSparseMatrix* NM_csc_to_csr(CSparseMatrix* csc)
{
#ifdef WITH_MKL_SPBLAS
  assert(csc);
  CHECK_MKL(load_mkl_function("mkl_dcsrcsc", (void**)&mkl_dcsrcsc_p));
  if(csc->m != csc->n)
  {
    fprintf(stderr, "NM_csc_to_csr :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  CSparseMatrix* csr = cs_spalloc(NSM_NROW_CSR(csc), NSM_NCOL_CSR(csc), csc->p[csc->n], 1, 0);
  assert(csr);
  csr->nz = -2;

  CS_INT n = csr->n;
  CS_INT job[6] = {0};
  CS_INT info = 0;
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
  if(csr->m != csr->n)
  {
    fprintf(stderr, "NM_csr_to_csc :: the matrix has to be square\n");
    exit(EXIT_FAILURE);
  }
  assert(csr);
  CSparseMatrix* csc = cs_spalloc(csr->m, csr->n, csr->nzmax, 1, 0);
  assert(csc);

  CS_INT n = csr->n;
  CS_INT job[6] = {0};
  job[5] = 1;
  CS_INT info = 0;
  (*mkl_dcsrcsc_p)(job, &n, csr->x, csr->i, csr->p, csc->x, csc->i, csc->p, &info);

  return csc;
#else
  if(csr->m != csr->n)
    {
      fprintf(stderr, "NM_csr_to_csc :: the matrix has to be square\n");
      exit(EXIT_FAILURE);
    }
  assert(csr);

  CSparseMatrix* csc = cs_spalloc(csr->m, csr->n, csr->nzmax, 1, 0);
  assert(csc);
  CS_INT* Ap = csr->p;
  CS_INT* Ai = csr->i;
  CS_ENTRY* Ax = csr->x;

  CS_INT* Bp = csc->p;
  CS_INT* Bi = csc->i;
  CS_ENTRY* Bx = csc->x;


  CS_INT n = csr->n;
  CS_INT nnz = csr->p[n];


  for (CS_INT col =0; col < csr->n+1 ; col++)
    {
      Bp[col] =0.0;
    }


  //compute number of non-zero entries per column of csr
  for (CS_INT i =0; i < nnz ; i++)
    {

      /* printf("Ai[%i] = %i ", i, Ai[i]); */
      Bp[Ai[i]]++;
    }

    // cumsum the nnz per column to get Bp[]
    CS_INT cumsum = 0;
    for (CS_INT col = 0; col < csr->n; col++) {
      CS_INT temp = Bp[col];
      Bp[col] = cumsum;
      cumsum += temp;
    }
    Bp[csr->n] = csr->nzmax;

    /* for (CS_INT col = 0; col < csr->n +1; col++) */
    /*   { */
    /* 	printf("Bp[%i] = %i ", col, Bp[col]); */
    /*   } */


    for(CS_INT row = 0; row < csr->m; row++){
        for(CS_INT jj = Ap[row]; jj < Ap[row+1]; jj++){
            CS_INT col  = Ai[jj];
            CS_INT dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }
    CS_INT last = 0;
    for(CS_INT col = 0; col <= csr->n; col++){
      CS_INT temp  = Bp[col];
      Bp[col] = last;
      last    = temp;
    }


    /* cs_print(csc,1); */
    return csc;
  /* fprintf(stderr, "NM_csr_to_csc :: MKL not enabled\n"); */
  /* exit(EXIT_FAILURE); */

#endif
}
