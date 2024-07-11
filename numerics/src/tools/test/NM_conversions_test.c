#include <stdio.h>  // for printf, fclose, fopen, FILE, NULL
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

#include <stdlib.h>  // for malloc

#include "CSparseMatrix.h"
#include "CSparseMatrix_internal.h"  // for cs_dl_entry, CS_INT, cs_dl_print
#include "NM_conversions.h"          // for NV_display
#include "NumericsFwd.h"             // for NumericsMatrix, NumericsSparseMatrix
#include "NumericsMatrix.h"          // for NM_entry, NM_display, NM_create
#include "NumericsSparseMatrix.h"    // for NumericsSparseMatrix, NSM_TRIPLET
#include "NumericsVector.h"          // for NV_display

static int compare_csc(CSparseMatrix *A, CSparseMatrix *B) {
  CS_INT m, n, nz, p, k, j, *Ap, *Ai, *Bi, *Bp;
  CS_ENTRY *Ax, *Bx;

  Ai = A->i;
  Ap = A->p;
  Ax = A->x;

  m = B->m;
  n = B->n;
  Bi = B->i;
  Bp = B->p;
  Bx = B->x;
  nz = B->nz;

  if (m != A->m) return 1;
  if (n != A->n) return 1;
  if (nz != A->nz) return 1;

  for (j = 0; j < n; j++) {
    if (Ap[j] != Bp[j]) {
      printf("p is not equal\n");
      return 1;
    }
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      if (Ai[p] != Bi[p]) {
        printf("i is not equal\n");
        return 1;
      }
      if (Ax[p] != Bx[p]) {
        printf("x is not equal\n");
        return 1;
      }
    }
  }
  return 0;
}

static int test_NM_conversion(NumericsMatrix *M) {
  CSparseMatrix *triplet = NM_triplet(M);

  printf("\n ###### conversion triplet to csc\n");
  CSparseMatrix *csc = NM_csc(M);
  /* CSparseMatrix_print(csc, 1); */

  printf("\n ###### conversion csc to triplet\n");
  CSparseMatrix *triplet_2 = NM_csc_to_triplet(csc);
  /* CSparseMatrix_print(triplet_2, 1); */

  int is_equal = CSparseMatrix_is_equal(triplet, triplet_2, 1e-14);

  if (!is_equal) {
    printf("triplet and triplet_2 are not equal\n");
    return 1;
  } else
    printf("triplet and triplet_2 are equal\n");

  printf("\n ###### conversion triplet to csr\n");
  CSparseMatrix *csr = NM_triplet_to_csr(NM_triplet(M));
  printf("csr from triplet \n:");
  /* CSparseMatrix_print(csr, 1); */

  is_equal = CSparseMatrix_is_equal(triplet, csr, 1e-14);

  if (!is_equal) {
    printf("triplet and csr are not equal\n");
    return 1;
  } else
    printf("triplet and csr are equal\n");

  printf("\n ###### conversion csr to csc\n");
  CSparseMatrix *csc_2 = NM_csr_to_csc(NM_csr(M));
  is_equal = CSparseMatrix_is_equal(csc, csc_2, 1e-14);
  if (!is_equal){
    printf("csc and csc_2 are not equal\n\n");
    return 1;
  } else
    printf("csc and csc_2 are equal\n\n");

  printf("\n ###### conversion csc to triplet\n");
  CSparseMatrix *triplet_3 = NM_csc_to_triplet(csc_2);

  /* CSparseMatrix_print(triplet_3, 1); */

  is_equal = CSparseMatrix_is_equal(triplet, triplet_3, 1e-14);

  if (!is_equal){
    printf("triplet and triplet_3 are not equal\n");
    return 1;}
  else
    printf("triplet and triplet_3 are equal\n");

  return 0;
};


static int test_read_write_sparse(NumericsMatrix *M_orig) {


  NumericsMatrix * M = NM_new();

  NM_copy(M_orig, M);


  FILE * foutput  =  fopen("matrix_sparse_triplet.dat", "w");
  NM_write_in_file(M, foutput);
  fclose(foutput);

  FILE * finput  =  fopen("matrix_sparse_triplet.dat", "r");
  NumericsMatrix * M_new = NM_new_from_file(finput);
  fclose(finput);

  int is_equal =  NM_equal(M, M_new);

  if (is_equal)
    printf("equal triplet\n");
  else
    printf("equal not triplet\n");

  CSparseMatrix * csc = NM_csc(M);

  NM_clearSparseStorage(M_new);
  M_new->matrix2->origin = NSM_CSC;
  M_new->matrix2->csc = csc;
  //NM_display(M_new);


  foutput  =  fopen("matrix_sparse_csc.dat", "w");
  NM_write_in_file(M_new, foutput);
  fclose(foutput);

  finput  =  fopen("matrix_sparse_csc.dat", "r");
  NumericsMatrix * M_csc = NM_new_from_file(finput);
  fclose(finput);

  is_equal =  NM_equal(M_new, M_csc);

  if (is_equal)
    printf("equal csc\n");
  else
    printf("equal not csc\n");

  CSparseMatrix * csr = NM_csr(M);

  NM_clearSparseStorage(M_new);
  M_new->matrix2->origin = NSM_CSR;
  M_new->matrix2->csr = csr;
  //NM_display(M_new);


  foutput  =  fopen("matrix_sparse_csr.dat", "w");
  NM_write_in_file(M_new, foutput);
  fclose(foutput);

  finput  =  fopen("matrix_sparse_csr.dat", "r");
  NumericsMatrix * M_csr = NM_new_from_file(finput);
  //NM_display(M_csr);
  fclose(finput);

  is_equal =  NM_equal(M_new, M_csr);

  if (is_equal)
    printf("equal csr\n");
  else
    printf("equal not csr\n");


  return 0;
};


int main() {
  int n = 5;
  int m = 5;

  NumericsMatrix *M = NM_create(NM_SPARSE, n, m);
  NM_triplet_alloc(M, 0);
  M->matrix2->origin = NSM_TRIPLET;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if ((i != j) && (i != j + 2)) NM_entry(M, i, j, i + j);
    }
  }

  test_read_write_sparse(M);




  NM_display(M);

  int info = test_NM_conversion(M);

  /* const char * filename =  "./data/NSM_csc_162x162.dat"; */
  /* NumericsMatrix *A = NM_new_from_filename(filename); */


  /* info = test_NM_conversion(A); */

  return info;
}
