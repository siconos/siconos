/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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

/*!\file test_globalFrictionContact_reformulation.c
 * \brief Example/test for globalFrictionContact_reduce_to_FrictionContact
 *        with various NumericsMatrix storage types in 2D and 3D.
 *
 * For each data file we build a reference reduced FC problem with dense
 * storage, then check that the reformulation succeeds for homogeneous
 * storage-type combinations (dense/dense, sparse/sparse, original/original)
 * and that mixed-type combinations fall back to NM_SPARSE for both M and H.
 */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FrictionContactProblem.h"
#include "GlobalFrictionContactProblem.h"
#include "NM_types.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"

typedef struct {
  const char* name;
  const char* filename;
  int dimension;
} TestProblem;

static const char* storage_type_name(NM_types type) {
  switch (type) {
    case NM_DENSE:
      return "NM_DENSE";
    case NM_SPARSE_BLOCK:
      return "NM_SPARSE_BLOCK";
    case NM_SPARSE:
      return "NM_SPARSE";
    default:
      return "NM_UNKNOWN";
  }
}

static NumericsMatrix* convert_to_type(NumericsMatrix* src, NM_types target) {
  if (src->storageType == target) {
    return NM_copy(src);
  }
  if (target == NM_DENSE) {
    NumericsMatrix* dst = NM_create(NM_DENSE, src->size0, src->size1);
    NM_to_dense(src, dst);
    return dst;
  }
  if (target == NM_SPARSE) {
    return NM_copy_to_sparse(src, DBL_EPSILON);
  }
  fprintf(stderr, "convert_to_type: conversion to %s is not supported\n",
          storage_type_name(target));
  return NULL;
}

static GlobalFrictionContactProblem* problem_with_storage(GlobalFrictionContactProblem* source,
                                                          NM_types M_type, NM_types H_type) {
  GlobalFrictionContactProblem* problem = globalFrictionContact_copy(source);
  assert(problem);

  NumericsMatrix* new_M = convert_to_type(problem->M, M_type);
  NumericsMatrix* new_H = convert_to_type(problem->H, H_type);
  assert(new_M);
  assert(new_H);

  NM_clear(problem->M);
  free(problem->M);
  NM_clear(problem->H);
  free(problem->H);

  problem->M = new_M;
  problem->H = new_H;
  return problem;
}

static int compare_vectors(int n, double* a, double* b, double tol) {
  for (int i = 0; i < n; ++i) {
    if (fabs(a[i] - b[i]) > tol) {
      return 0;
    }
  }
  return 1;
}

static int compare_matrices(NumericsMatrix* A, NumericsMatrix* B, double tol) {
  assert(A->size0 == B->size0);
  assert(A->size1 == B->size1);
  if (A->storageType == NM_DENSE && B->storageType == NM_DENSE) {
    return compare_vectors(A->size0 * A->size1, A->dense, B->dense, tol);
  }

  NumericsMatrix* A_dense = NULL;
  NumericsMatrix* B_dense = NULL;
  int own_A = 0;
  int own_B = 0;

  if (A->storageType == NM_DENSE) {
    A_dense = A;
  } else {
    A_dense = NM_create(NM_DENSE, A->size0, A->size1);
    NM_to_dense(A, A_dense);
    own_A = 1;
  }

  if (B->storageType == NM_DENSE) {
    B_dense = B;
  } else {
    B_dense = NM_create(NM_DENSE, B->size0, B->size1);
    NM_to_dense(B, B_dense);
    own_B = 1;
  }

  int same = compare_vectors(A->size0 * A->size1, A_dense->dense, B_dense->dense, tol);

  if (own_A) {
    NM_clear(A_dense);
    free(A_dense);
  }
  if (own_B) {
    NM_clear(B_dense);
    free(B_dense);
  }
  return same;
}

static void test_combination(GlobalFrictionContactProblem* source, NM_types M_type,
                             NM_types H_type, FrictionContactProblem* reference, double tol) {
  printf("  Testing M=%s, H=%s\n", storage_type_name(M_type), storage_type_name(H_type));
  GlobalFrictionContactProblem* problem = problem_with_storage(source, M_type, H_type);

  FrictionContactProblem* reduced = globalFrictionContact_reduce_to_FrictionContact(problem);

  assert(reduced);
  assert(reduced->numberOfContacts == problem->numberOfContacts);
  assert(reduced->dimension == problem->dimension);
  assert(reduced->M);
  assert(reduced->M->size0 == problem->dimension * problem->numberOfContacts);
  assert(reduced->M->size1 == problem->dimension * problem->numberOfContacts);
  assert(reduced->q);

  /* The reformulation function may fall back to NM_SPARSE internally, but it
   * must restore the original matrices in the caller's problem. */
  assert(problem->M->storageType == M_type);
  assert(problem->H->storageType == H_type);

  assert(compare_matrices(reduced->M, reference->M, tol));
  assert(compare_vectors(reduced->M->size0, reduced->q, reference->q, tol));

  frictionContactProblem_free(reduced);
  globalFrictionContact_free(problem);
}

int main() {
  verbose = 0;

  TestProblem problems[] = {
      {"GFC2D_Example0", "./data/GFC3D_Example0_2d.dat", 2},
      {"GFC2D_Example1", "./data/GFC3D_Example1_2d.dat", 2},
      {"GFC2D_OneContact", "./data/GFC3D_OneContact_2d.dat", 2},
      {"GFC3D_Example0", "./data/GFC3D_Example0.dat", 3},
      {"GFC3D_Example1", "./data/GFC3D_Example1.dat", 3},
      {"GFC3D_OneContact", "./data/GFC3D_OneContact.dat", 3},
      {"GFC3D_Example0_SBM", "./data/GFC3D_Example0_SBM.dat", 3},
  };
  const int n_problems = (int)(sizeof(problems) / sizeof(problems[0]));
  const double tol = 1e-6;

  printf("Running globalFrictionContact reformulation tests (%d cases)\n", n_problems);
  for (int p = 0; p < n_problems; ++p) {
    printf("[%d/%d] %s (%dD) from %s\n", p + 1, n_problems, problems[p].name,
           problems[p].dimension, problems[p].filename);

    GlobalFrictionContactProblem* source =
        globalFrictionContact_new_from_filename(problems[p].filename);
    assert(source);
    assert(source->dimension == problems[p].dimension);
    assert(source->numberOfContacts > 0);

    /* Reference reduced problem computed with dense storage. */
    GlobalFrictionContactProblem* ref_problem =
        problem_with_storage(source, NM_DENSE, NM_DENSE);
    FrictionContactProblem* reference =
        globalFrictionContact_reduce_to_FrictionContact(ref_problem);
    assert(reference);
    globalFrictionContact_free(ref_problem);

    /* Original storage types from the data file. */
    printf("  Original M=%s, H=%s\n", storage_type_name(source->M->storageType),
           storage_type_name(source->H->storageType));
    test_combination(source, source->M->storageType, source->H->storageType, reference, tol);

    /* Homogeneous dense and sparse combinations. */
    test_combination(source, NM_DENSE, NM_DENSE, reference, tol);
    test_combination(source, NM_SPARSE, NM_SPARSE, reference, tol);

    /* Mixed combinations: must fall back to NM_SPARSE. */
    test_combination(source, NM_DENSE, NM_SPARSE, reference, tol);
    test_combination(source, NM_SPARSE, NM_DENSE, reference, tol);

    frictionContactProblem_free(reference);
    globalFrictionContact_free(source);
    printf("  OK\n");
  }

  printf("\nSUCCESS: all %d reformulation tests passed.\n", n_problems);
  return 0;
}
