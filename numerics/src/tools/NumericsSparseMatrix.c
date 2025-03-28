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
#include "NumericsSparseMatrix.h"

#include <assert.h>  // for assert

#include "CSparseMatrix.h"  // for CSparseMatrix, CS_INT, cs_dl_spfree
#ifndef __cplusplus
#include <stdbool.h>  // for bool, false, true
#endif
#include <stdio.h>   // for NULL, size_t, printf, fprintf, fscanf
#include <stdlib.h>  // for free, exit, malloc, realloc, EXIT_FAILURE

#include "NumericsMatrix.h"  // for NumericsMatrix, NM_csc, numericsSparse...
#include "SiconosConfig.h"   // for HAVE_SORT
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "numerics_verbose.h"  // for numerics_error_nonfatal, CHECK_IO
#include "siconos_debug.h"     // for DEBUG_BEGIN, DEBUG_END, DEBUG_EXPR
#include "string.h"            // for memcpy, memset

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#elif !(__INTEL_COMPILER || __APPLE__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

typedef struct {
  CS_INT i;
  size_t indx;
} sort_indices_struct;

#ifdef HAVE_SORT
#define SORT_NAME sorter
#define SORT_TYPE sort_indices_struct
#define SORT_CMP(x, y) ((x).i - (y).i)
#include "sort.h"  // for sorter_tim_sort
#else
#include "stdlib.h"  // qsort
static int sort_indices_struct_cmp(const void* a, const void* b) {
  const sort_indices_struct* sa = (const sort_indices_struct*)a;
  const sort_indices_struct* sb = (const sort_indices_struct*)b;
  return (sa->i > sb->i) - (sa->i < sb->i);
}
#endif

#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__)
#pragma GCC diagnostic pop
#endif

void NSM_set_version(NumericsSparseMatrix* M, NSM_t type, version_t value) {
  NDV_set_value(&(M->versions[type]), value);
}

void NSM_inc_version(NumericsSparseMatrix* M, NSM_t type) {
  version_t new_version = NSM_max_version(M) + 1;
  NSM_set_version(M, type, new_version);
}

/* internal compare function */

CSparseMatrix* NSM_latest(const NumericsSparseMatrix* M) {
  assert(M);

  switch (NSM_latest_id(M)) {
    case NSM_TRIPLET:
      return M->triplet;
    case NSM_HALF_TRIPLET:
      return M->half_triplet;
    case NSM_CSR:
      return M->csr;
    case NSM_CSC:
      return M->csc;
    default:
      numerics_error("NSM_latest", "unknown matrix type");
  }
  return 0;
}

void NSM_reset_version(NumericsSparseMatrix* M, NSM_t id) { NDV_reset(&(M->versions[id])); }

void NSM_reset_versions(NumericsSparseMatrix* M) {
  NSM_reset_version(M, NSM_TRIPLET);
  NSM_reset_version(M, NSM_HALF_TRIPLET);
  NSM_reset_version(M, NSM_CSC);
  NSM_reset_version(M, NSM_CSR);
}

void NSM_null(NumericsSparseMatrix* A) {
  A->linearSolverParams = NULL;
  A->triplet = NULL;
  A->half_triplet = NULL;
  A->csc = NULL;
  A->trans_csc = NULL;
  A->csr = NULL;
  A->diag_indx = NULL;
  A->origin = NSM_UNKNOWN;

  NSM_reset_versions(A);
}

double* NSM_data(NumericsSparseMatrix* A) {
  assert(NSM_version(A, NSM_latest_id(A)) == NSM_version(A, A->origin));

  switch (A->origin) {
    case NSM_CSC: {
      assert(A->csc);
      return A->csc->x;
      break;
    }
    case NSM_CSR: {
      assert(A->csr);
      return A->csr->x;
      break;
    }
    case NSM_TRIPLET: {
      assert(A->triplet);
      return A->triplet->x;
      break;
    }
    case NSM_HALF_TRIPLET: {
      assert(A->half_triplet);
      return A->half_triplet->x;
      break;
    }
    default:
      printf("NSM_data :: unknown matrix origin %d", A->origin);
      exit(EXIT_FAILURE);
  }
}

NumericsSparseMatrix* NSM_new(void) {
  NumericsSparseMatrix* p = (NumericsSparseMatrix*)malloc(sizeof(NumericsSparseMatrix));

  NSM_null(p);
  p->linearSolverParams = NSM_linearSolverParams_new();

  return p;
}

NumericsSparseMatrix* NSM_clear(NumericsSparseMatrix* A) {
  if (A->linearSolverParams) {
    NSM_linearSolverParams_free(A->linearSolverParams);
    A->linearSolverParams = NULL;
  }
  if (A->triplet) {
    cs_spfree(A->triplet);
    A->triplet = NULL;
  }
  if (A->half_triplet) {
    cs_spfree(A->half_triplet);
    A->half_triplet = NULL;
  }
  if (A->csc) {
    cs_spfree(A->csc);
    A->csc = NULL;
  }
  if (A->trans_csc) {
    cs_spfree(A->trans_csc);
    A->trans_csc = NULL;
  }
  if (A->csr) {
    cs_spfree(A->csr);
    A->csr = NULL;
  }
  if (A->diag_indx) {
    free(A->diag_indx);
    A->diag_indx = NULL;
  }

  NSM_reset_versions(A);

  return NULL;
}

void NSM_version_copy(const NumericsSparseMatrix* const A, NumericsSparseMatrix* B) {
  assert(A);
  assert(B);
  switch (A->origin) {
    case NSM_TRIPLET: {
      NSM_set_version(B, NSM_TRIPLET, NSM_version(A, NSM_TRIPLET));
      break;
    }
    case NSM_HALF_TRIPLET: {
      NSM_set_version(B, NSM_HALF_TRIPLET, NSM_version(A, NSM_HALF_TRIPLET));
      break;
    }
    case NSM_CSR: {
      NSM_set_version(B, NSM_CSR, NSM_version(A, NSM_CSR));
      break;
    }
    case NSM_CSC: {
      NSM_set_version(B, NSM_CSC, NSM_version(B, NSM_CSC));
      break;
    }
    default: {
      numerics_error("NSM_version_copy", "unknown id");
    }
      assert(false);
  }
}

void NSM_copy(NumericsSparseMatrix* A, NumericsSparseMatrix* B) {
  assert(A);
  assert(B);

  B->origin = A->origin;

  CSparseMatrix* A_;
  CSparseMatrix* B_;

  switch (A->origin) {
    case NSM_TRIPLET: {
      A_ = A->triplet;

      if (!B->triplet) {
        B->triplet = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 1);
      }

      B_ = B->triplet;
      NSM_set_version(B, NSM_TRIPLET, NSM_version(A, NSM_TRIPLET));
      break;
    }
    case NSM_HALF_TRIPLET: {
      A_ = A->half_triplet;

      if (!B->half_triplet) {
        B->half_triplet = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 1);
      }

      B_ = B->half_triplet;
      NSM_set_version(B, NSM_HALF_TRIPLET, NSM_version(A, NSM_HALF_TRIPLET));
      break;
    }
    case NSM_CSC: {
      assert(A->csc);

      A_ = A->csc;

      if (!B->csc) {
        B->csc = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 0);
      }

      B_ = B->csc;
      NSM_set_version(B, NSM_CSC, NSM_version(A, NSM_CSC));
      break;
    }
    case NSM_CSR: {
      assert(A->csr);

      A_ = A->csr;

      if (!B->csr) {
        B->csr = cs_spalloc(A_->m, A_->n, A_->nzmax, 1, 0);
        B->csr->nz = -2;
        B->csr->m = A_->m;
        B->csr->n = A_->n;
      }

      B_ = B->csr;
      NSM_set_version(B, NSM_CSR, NSM_version(A, NSM_CSR));
      break;
    }
    default: {
      fprintf(stderr, "NSM_copy :: error unknown origin %u for sparse matrix\n", A->origin);
      exit(EXIT_FAILURE);
    }
  }
  CSparseMatrix_copy(A_, B_);

  /* diag_indx   */
  /* We remove diag_indx from B and  we copy it from A if it exists */
  if (B->diag_indx) {
    free(B->diag_indx);
    B->diag_indx = NULL;
  }
  if (A->diag_indx) {
    B->diag_indx = (CS_INT*)malloc((size_t)A_->m * sizeof(CS_INT));
    memcpy(B->diag_indx, A->diag_indx, (size_t)A_->m * sizeof(CS_INT));
  }

  /*  NSM_linear_solver_params */
  if (B->linearSolverParams) {
    NSM_linearSolverParams_free(B->linearSolverParams);
    B->linearSolverParams = NULL;
  }
  /* We copy only the id of the solver to keep track of if NM_preserve is used before
   * factorization. The remaining part is not copied since we have blind pointer.
   */
  if (A->linearSolverParams) {
    B->linearSolverParams = NSM_linearSolverParams_new();
    B->linearSolverParams->solver = A->linearSolverParams->solver;
    B->linearSolverParams->LDLT_solver = A->linearSolverParams->LDLT_solver;
  }
}

NSM_linear_solver_params* NSM_linearSolverParams_new(void) {
  NSM_linear_solver_params* p =
      (NSM_linear_solver_params*)malloc(sizeof(NSM_linear_solver_params));

#if defined(WITH_MUMPS)
  p->solver = NSM_MUMPS;
#elif defined(WITH_UMFPACK)
  p->solver = NSM_UMFPACK;
#elif defined(WITH_SUPERLU)
  p->solver = NSM_SUPERLU;
#elif defined(WITH_SUPERLU_MT)
  p->solver = NSM_SUPERLU_MT;
#elif defined(WITH_MKL_PARDISO)
  p->solver = NSM_MKL_PARDISO;
#else
  p->solver = NSM_CSPARSE;  // default solver
#endif

  p->LDLT_solver = p->solver;
#if defined(WITH_MA57)
  p->LDLT_solver = NSM_HSL;
#endif

  p->parent_matrix = NULL;

  p->linear_solver_data = NULL;
  p->solver_free_hook = NULL;

  p->iWorkSize = 0;
  p->dWorkSize = 0;

  p->iWork = NULL;
  p->dWork = NULL;
  p->linalg_data = NULL;

  return p;
}

NSM_linear_solver_params* NSM_linearSolverParams(NumericsMatrix* A) {
  if (!numericsSparseMatrix(A)->linearSolverParams) {
    numericsSparseMatrix(A)->linearSolverParams = NSM_linearSolverParams_new();
  }
  return numericsSparseMatrix(A)->linearSolverParams;
}

NSM_linear_solver_params* NSM_linearSolverParams_free(NSM_linear_solver_params* p) {
  /* First free linear_solver_data if some additional information has been given  */
  if (p->solver_free_hook) {
    (*p->solver_free_hook)(p);
    p->solver_free_hook = NULL;
  }

  if (p->iWork) {
    // assert(p->iWorkSize>0); This assertion is commented since it may happen that we
    //                         allocate a chunk of memory of size 0 (for an empty matrix !!)
    free(p->iWork);
    p->iWork = NULL;
  }

  if (p->dWork) {
    // assert(p->dWorkSize>0); This assertion is commented since it may happen that we
    //                         allocate a chunk of memory of size 0 (for an empty matrix !!)
    free(p->dWork);
    p->dWork = NULL;
  }

  if (p->linear_solver_data) {
    free(p->linear_solver_data);
    p->linear_solver_data = NULL;
  }

  if (p->linalg_data) {
    p->linalg_data->free_fn(p->linalg_data);
    free(p->linalg_data);
    p->linalg_data = NULL;
  }

  free(p);
  return NULL;
}

void NSM_clear_p(void* p) {
  assert(p);
  NSM_linear_solver_params* ptr = (NSM_linear_solver_params*)p;
  CSparseMatrix_factors* cs_lu_A = (CSparseMatrix_factors*)NSM_linear_solver_data(ptr);

  CSparseMatrix_free_lu_factors(cs_lu_A);

  ptr->linear_solver_data = NULL;
}

size_t NSM_nnz(const CSparseMatrix* const A) {
  if (A->nz >= 0) {
    return (size_t)A->nz;
  } else if (A->nz == NSM_CS_CSC) {
    return (size_t)A->p[A->n];
  } else if (A->nz == NSM_CS_CSR) {
    return (size_t)A->p[A->m];
  } else {
    fprintf(stderr, "NSM_nnz :: unsupported nz number " CS_ID, A->nz);
    exit(EXIT_FAILURE);
  }
}

void NSM_fix_csc(CSparseMatrix* A) {
  CS_INT* Ap = A->p;
  CS_INT* Ai = A->i;
  double* xbck = NULL;
  sort_indices_struct* s = NULL;
  for (size_t j = 0; j < (size_t)A->n; ++j) {
    bool need_sorting = false;
    CS_INT max_indx = -1;
    CS_INT p = Ap[j];
    for (; p < Ap[j + 1]; ++p) {
      if (Ai[p] <= max_indx) {
        need_sorting = true;
        break;
      } else {
        max_indx = Ai[p];
      }
    }
    if (need_sorting) {
      double* Ax = A->x;
      CS_INT min_indx = Ai[p];
      CS_INT ps = p - 1;
      for (; ps > Ap[j]; --ps) {
        if (Ai[ps] < min_indx) {
          break;
        }
      }
      size_t len = Ap[j + 1] - ps;
      s = (sort_indices_struct*)realloc(s, len * sizeof(sort_indices_struct));
      xbck = (double*)realloc(xbck, len * sizeof(double));
      memcpy(xbck, &Ax[ps], len * sizeof(double));
      for (size_t i = 0, pp = ps; i < len; ++i, ++pp) {
        s[i].i = Ai[pp];
        s[i].indx = i;
      }

#ifdef HAVE_SORT
      sorter_tim_sort(s, len);
#else
      qsort(s, len, sizeof(sort_indices_struct), sort_indices_struct_cmp);
#endif

      for (size_t i = 0, pp = ps; i < len; ++i, ++pp) {
        Ai[pp] = s[i].i;
        Ax[pp] = xbck[s[i].indx];
      }
    }
  }

  if (xbck) {
    free(xbck);
    xbck = NULL;
  }
  if (s) {
    free(s);
    s = NULL;
  }
}

void NSM_sort_csc(CSparseMatrix* A) {
  CSparseMatrix* AT = cs_transpose(A, 1);
  CSparseMatrix* ATT = cs_transpose(AT, 1);
  cs_spfree(AT);

  free(A->i);
  free(A->x);
  A->i = ATT->i;
  A->x = ATT->x;

  ATT->i = NULL;
  ATT->x = NULL;

  cs_spfree(ATT);
}

int NSM_to_dense(const NumericsSparseMatrix* const A, double* B) {
  if (!A) {
    printf("NSM_to_dense :: A = null\n");
    return (0);
  }
  return (int)CSparseMatrix_to_dense(NSM_get_origin(A), B);
}

unsigned NSM_origin(const NumericsSparseMatrix* M) {
  assert(M);
  if (!M) return -1;
  assert(NSM_version(M, NSM_latest_id(M)) == NSM_version(M, M->origin));
  return M->origin;
}

CSparseMatrix* NSM_get_origin(const NumericsSparseMatrix* M) {
  assert(M);
  switch (M->origin) {
    case NSM_CSC:
      return M->csc;
    case NSM_TRIPLET:
      return M->triplet;
    case NSM_HALF_TRIPLET:
      return M->half_triplet;
    case NSM_CSR:
      return M->csr;
    default:
      numerics_error_nonfatal("NSM_get_origin", "Unknown matrix origin %d", M->origin);
      return NULL;
  }
}

void NSM_write_in_file(const NumericsSparseMatrix* m, FILE* file) {
  assert(m);
  fprintf(file, "%d\n", NSM_origin(m));
  CSparseMatrix_print_in_file(NSM_get_origin(m), 0, file);
}

NumericsSparseMatrix* NSM_new_from_file(FILE* file) {
  int info;
  NSM_t _origin = 0;
  CHECK_IO(fscanf(file, "%d", &_origin), &info);
  NumericsSparseMatrix* out = NSM_new();
  out->origin = _origin;

  CSparseMatrix* C = CSparseMatrix_new_from_file(file);

  if (C->nz >= 0) {
    assert(out->origin == NSM_TRIPLET);
    out->triplet = C;
    out->origin = NSM_TRIPLET;
    NSM_inc_version(out, NSM_TRIPLET);
  } else {
    if (out->origin == NSM_CSC) {
      out->csc = C;
      out->origin = NSM_CSC;
      NSM_inc_version(out, NSM_CSC);
    } else if (out->origin == NSM_CSR) {
      out->csr = C;
      out->origin = NSM_CSR;
      NSM_inc_version(out, NSM_CSR);
    } else if (out->origin == NSM_HALF_TRIPLET) {
      out->half_triplet = C;
      NSM_inc_version(out, NSM_HALF_TRIPLET);
    }
  }
  return out;
}

NumericsSparseMatrix* NSM_triplet_eye(unsigned int size) {
  NSM_t _origin = NSM_TRIPLET;
  NumericsSparseMatrix* out = NSM_new();
  out->origin = _origin;

  CSparseMatrix* C = cs_spalloc(size, size, size, 1, 1);

  for (unsigned int k = 0; k < size; k++) {
    C->nz++;
    C->i[k] = k;
    C->p[k] = k;
    C->x[k] = 1.0;
  }
  assert(out->origin == NSM_TRIPLET);
  out->triplet = C;
  out->origin = NSM_TRIPLET;
  NSM_inc_version(out, NSM_TRIPLET);
  return out;
}

NumericsSparseMatrix* NSM_triplet_scalar(unsigned int size, double s) {
  int _origin = NSM_TRIPLET;
  NumericsSparseMatrix* out = NSM_new();
  out->origin = _origin;

  CSparseMatrix* C = cs_spalloc(size, size, size, 1, 1);

  for (unsigned int k = 0; k < size; k++) {
    C->nz++;
    C->i[k] = k;
    C->p[k] = k;
    C->x[k] = s;
  }
  assert(out->origin == NSM_TRIPLET);
  out->triplet = C;
  out->origin = NSM_TRIPLET;
  return out;
}

static CS_INT* NSM_diag_indices_trivial(NumericsMatrix* M) {
  NumericsSparseMatrix* A = M->matrix2;
  assert(A);
  if (A->diag_indx) return A->diag_indx;

  CS_INT* indices = (CS_INT*)malloc((size_t)M->size0 * sizeof(CS_INT));
  A->diag_indx = indices;
  /* XXX hack --xhub  */
  if (A->origin == NSM_TRIPLET) {
    NM_csc(M);
    A->origin = NSM_CSC;
  }
  switch (A->origin) {
    case NSM_CSC: {
      assert(A->csc);

      CS_INT* Ai = A->csc->i;
      CS_INT* Ap = A->csc->p;

      for (CS_INT j = 0; j < (CS_INT)M->size0; ++j) {
        int is_diag_index_found = 0;
        for (CS_INT p = Ap[j]; p < Ap[j + 1]; ++p) {
          if (Ai[p] == j) {
            indices[j] = p;
            is_diag_index_found = 1;
            break;
          }
        }
        if (!is_diag_index_found) {
          free(indices);
          A->diag_indx = NULL;
          return NULL;
        }
      }
      break;
    }
    case NSM_CSR: {
      assert(A->csr);

      CS_INT* Ai = A->csr->i;
      CS_INT* Ap = A->csr->p;

      for (CS_INT j = 0; j < (CS_INT)M->size0; ++j) {
        int is_diag_index_found = 0;
        for (CS_INT p = Ap[j]; p < Ap[j + 1]; ++p) {
          if (Ai[p] == j) {
            indices[j] = p;
            is_diag_index_found = 1;
            break;
          }
        }
        if (!is_diag_index_found) {
          free(indices);
          A->diag_indx = NULL;
          return NULL;
        }
      }
      break;
    }
    case NSM_TRIPLET:
    default:
      printf("NSM_diag_indices :: unknown matrix origin %d", A->origin);
      exit(EXIT_FAILURE);
  }

  return indices;
}

CS_INT* NSM_diag_indices(NumericsMatrix* M) {
  DEBUG_BEGIN("NSM_diag_indices(NumericsMatrix* M)\n");
  NumericsSparseMatrix* A = M->matrix2;
  assert(A);
  if (A->diag_indx) return A->diag_indx;

  /* 1-  we assume that all diagonal elements exist, and we search it in a trivial way */
  CS_INT* indices = NSM_diag_indices_trivial(M);
  if (indices) {
    return indices;
  }

  /* 2- if not, we compute the diagonal indices in the general case without any assumption
   * on the existence  of the diagonal elements. This makes sure that the matrix is replaced
   * by one that has all diagonal elements. The copy of the matrix could be optimized. */

  indices = (CS_INT*)malloc((size_t)M->size0 * sizeof(CS_INT));
  A->diag_indx = indices;
  /* XXX hack --xhub  */
  if (A->origin == NSM_TRIPLET) {
    NM_csc(M);
    A->origin = NSM_CSC;
  }
  switch (A->origin) {
    case NSM_CSC: {
      assert(A->csc);
      indices[0] = 0;
      CSparseMatrix* newMat =
          cs_spalloc(M->size0, M->size1, A->csc->p[M->size0] + M->size0, 1, 0);
      CS_INT* Ai = A->csc->i;
      CS_INT* Ap = A->csc->p;
      double* Ax = A->csc->x;
      CS_INT* Ni = newMat->i;
      CS_INT* Np = newMat->p;
      double* Nx = newMat->x;
      CS_INT end = Ap[1];
      CS_INT inc = 0;
      Np[0] = 0;
      if (Ai[0] == 0) {
        memcpy(Ni, Ai, (size_t)end * sizeof(CS_INT));
        Np[1] = Ap[1];
        memcpy(Nx, Ax, (size_t)end * sizeof(double));
      } else {
        Ni[0] = 0;
        Np[1] = Ap[1] + 1;
        Nx[0] = 0.;
        memcpy(&Ni[1], Ai, (size_t)end * sizeof(CS_INT));
        memcpy(&Nx[1], Ax, (size_t)end * sizeof(double));
        ++inc;
      }

      /* Could optimize further and copy everything using memcpy */
      for (CS_INT j = 1; j < M->size0; ++j) {
        CS_INT rem = 0;
        for (CS_INT p = Ap[j]; (rem == 0) && (p < Ap[j + 1]); ++p) {
          if (Ai[p] < j) {
            Ni[p + inc] = Ai[p];
            Nx[p + inc] = Ax[p];
          } else {
            if (Ai[p] > j) {
              Ni[p + inc] = j;
              Nx[p + inc] = 0.;
              indices[j] = p + inc;
              ++inc;
            } else {
              indices[j] = p + inc;
            }
            rem = p;
            Np[j] = Ap[j] + inc;
          }
          end = Ap[j + 1] - rem;
          memcpy(&Ni[rem + inc], &Ai[rem], (size_t)end * sizeof(CS_INT));
          memcpy(&Nx[rem + inc], &Ax[rem], (size_t)end * sizeof(double));
          assert(inc <= M->size0);
        }
      }
      Np[M->size0] = Ap[M->size0] + inc;
      NM_clearSparseStorage(M);
      A->origin = NSM_CSC;
      A->csc = newMat;
      DEBUG_EXPR(NM_display(M););
      break;
    }
    case NSM_TRIPLET:
    case NSM_CSR:
    default:

      exit(EXIT_FAILURE);
  }
  DEBUG_END("NSM_diag_indices(NumericsMatrix* M)\n");
  return indices;
}

void NSM_extract_block(NumericsMatrix* M, double* blockM, size_t pos_row, size_t pos_col,
                       size_t block_row_size, size_t block_col_size) {
  assert(M);
  assert(M->storageType == NM_SPARSE);
  assert(blockM);
  assert(pos_row < (size_t)M->size0);
  assert(pos_col < (size_t)M->size1);
  assert(block_row_size > 0 && block_row_size + pos_row <= (unsigned long int)M->size0);
  assert(block_col_size > 0 && block_col_size + pos_col <= (unsigned long int)M->size1);

  assert(M->matrix2);

  /* Clear memory */
  memset(blockM, 0, block_row_size * block_col_size * sizeof(double));

  switch (M->matrix2->origin) {
    case NSM_TRIPLET:
    case NSM_CSC: {
      CSparseMatrix* Mcsc = NM_csc(M);
      assert(Mcsc);
      CS_INT* Mp = Mcsc->p;
      CS_INT* Mi = Mcsc->i;
      double* Mx = Mcsc->x;
      for (size_t j = pos_col; j < pos_col + block_col_size; ++j) {
        for (CS_INT p = Mp[j]; p < Mp[j + 1]; ++p) {
          CS_INT row_nb = Mi[p];
          /* Warning : the following strategy work only if the csc storage is
             correclty ordered. Use NSM_fix_csc to be sure*/
          if (row_nb >= (CS_INT)pos_row) {
            if (row_nb >= (CS_INT)(pos_row + block_row_size)) {
              break;
            } else {
              blockM[(j - pos_col) * block_col_size + row_nb - pos_row] = Mx[p];
            }
          }
        }
      }
      break;
    }
    case NSM_CSR: {
      CSparseMatrix* Mcsr = M->matrix2->csr;
      assert(Mcsr);
      CS_INT* Mp = Mcsr->p;
      CS_INT* Mi = Mcsr->i;
      double* Mx = Mcsr->x;

      for (size_t i = pos_row; i < pos_row + block_row_size; ++i) {
        for (CS_INT p = Mp[i]; p < Mp[i + 1]; ++p) {
          CS_INT col_nb = Mi[p];
          /* Warning : the following strategy work only if the csc storage is
             correclty ordered. Use NSM_fix_csc to be sure*/
          if (col_nb >= (CS_INT)pos_col) {
            if (col_nb >= (CS_INT)(pos_col + block_col_size)) {
              break;
            } else {
              blockM[(col_nb - pos_col) * block_col_size + i - pos_row] = Mx[p];
            }
          }
        }
      }
      break;
    }
    default:
      printf("NSM_extract_block :: unsupported matrix type %d\n", M->matrix2->origin);
      exit(EXIT_FAILURE);
  }
}

double** NSM_extract_diagonal_blocks(NumericsMatrix* M, size_t block_size) {
  assert(M);
  assert(M->storageType == NM_SPARSE);
  assert(M->matrix2);

  double** diag_blocks = NULL;
  CS_INT number_of_diagonal_blocks = 0;
  switch (M->matrix2->origin) {
    case NSM_TRIPLET:
    case NSM_CSC: {
      CSparseMatrix* Mcsc = NM_csc(M);
      assert(Mcsc);

      CS_INT n = Mcsc->n;
      CS_INT m = Mcsc->m;

      if (n != m) return diag_blocks;

      if (n % block_size != 0) {
        return diag_blocks;
      } else {
        number_of_diagonal_blocks = n / block_size;
      }

      break;
    }
    case NSM_CSR: {
      CSparseMatrix* Mcsr = M->matrix2->csr;
      assert(Mcsr);
      CS_INT n = Mcsr->n;
      CS_INT m = Mcsr->m;

      if (n != m) return diag_blocks;

      if (n % block_size != 0) {
        return diag_blocks;
      } else {
        number_of_diagonal_blocks = n / block_size;
      }

      break;
    }
    default:
      printf("NSM_extract_block :: unsupported matrix type %d\n", M->matrix2->origin);
      exit(EXIT_FAILURE);
  }

  diag_blocks = (double**)malloc(number_of_diagonal_blocks * sizeof(double*));

  for (CS_INT b = 0; b < number_of_diagonal_blocks; b++) {
    diag_blocks[b] = (double*)malloc(block_size * block_size * sizeof(double));

    size_t pos_row = b * block_size;
    size_t pos_col = b * block_size;

    NSM_extract_block(M, diag_blocks[b], pos_row, pos_col, block_size, block_size);
  }

  return diag_blocks;
}

CSparseMatrix* NSM_remove_diagonal_blocks(NumericsMatrix* M, size_t block_size) {
  assert(M);
  assert(M->storageType == NM_SPARSE);
  assert(M->matrix2);

  CSparseMatrix* out = cs_spalloc(M->size0, M->size1, 0, 1, 1);

  CSparseMatrix* M_triplet = NM_triplet(M);

  CS_INT* Mp = M_triplet->p;
  CS_INT* Mi = M_triplet->i;
  CS_ENTRY* Mx = M_triplet->x;

  CS_INT nz = M_triplet->nz;

  for (CS_INT e = 0; e < nz; e++) {
    CS_INT i = Mi[e];
    CS_INT j = Mp[e];

    int is_block_diagonal_element = 0;

    int index = i % block_size;

    if ((j >= i - index) && (j < i - index + block_size)) {
      is_block_diagonal_element = 1;
    }

    /* if (is_block_diagonal_element) */
    /* 	printf("%i %i is diagonal \n", i, j); */
    /* else */
    /* 	printf("%i %i is NOT diagonal \n", i, j); */

    if (!is_block_diagonal_element) {
      CSparseMatrix_entry(out, i, j, Mx[e]);  // to be improved
    }
  }

  return out;
}
