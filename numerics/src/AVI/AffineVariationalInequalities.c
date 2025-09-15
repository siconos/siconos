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
#ifndef AVI_PROBLEM_C
#define AVI_PROBLEM_C

#include "AffineVariationalInequalities.h"

#include <assert.h>  // for assert
#include <stdio.h>   // for printf, NULL, fprintf, fscanf, FILE
#include <stdlib.h>  // for free, malloc, exit, EXIT_FAILURE

#include "NumericsMatrix.h"    // for NM_display, NM_clear, NM_new_from_file
#include "SiconosSets.h"       // for polyhedron_set, free_polyhedron, free_...
#include "numerics_verbose.h"  // for CHECK_IO

void AVI_display(AffineVariationalInequalities* avi) {
  assert(avi);
  size_t i;
  size_t n = avi->size;
  printf("AffineVariationalInequalities Display :\n-------------\n");
  printf("size :%zu \n", avi->size);
  if (avi->M) {
    printf("M matrix:\n");
    NM_display(avi->M);
  } else
    printf("No M matrix.\n");

  if (avi->q) {
    printf("q vector:\n");
    for (i = 0; i < n; i++) printf("q[ %zd ] = %12.8e\n", i, avi->q[i]);
  } else
    printf("No q vector.\n");

  if (avi->d) {
    printf("d vector:\n");
    for (i = 0; i < n; i++) printf("d[ %zd ] = %12.8e\n", i, avi->d[i]);
  } else
    printf("No d vector.\n");
}

int AVI_printInFile(AffineVariationalInequalities* avi, FILE* file) {
  if (!avi) {
    fprintf(stderr,
            "Numerics, AffineVariationalInequalities printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;
  size_t n = avi->size;
  fprintf(file, "%zu\n", n);
  NM_write_in_file(avi->M, file);
  for (i = 0; i < avi->M->size1; i++) {
    fprintf(file, "%32.24e ", avi->q[i]);
  }
  return 1;
}

int AVI_newFromFile(AffineVariationalInequalities* avi, FILE* file) {
  size_t n = 0;
  int i;

  CHECK_IO(fscanf(file, "%zu\n", &n));
  avi->size = n;
  avi->M = NM_new_from_file(file);

  avi->q = (double*)malloc(avi->M->size1 * sizeof(double));
  for (i = 0; i < avi->M->size1; i++) {
    CHECK_IO(fscanf(file, "%lf ", &(avi->q[i])));
  }
  return 1;
}
int AVI_newFromFilename(AffineVariationalInequalities* avi, char* filename) {
  int info = 0;
  FILE* file = fopen(filename, "r");

  info = AVI_newFromFile(avi, file);

  fclose(file);
  return info;
}

void freeAVI(AffineVariationalInequalities* avi) {
  assert(avi);
  if (avi->M) {
    NM_clear(avi->M);
    free(avi->M);
    avi->M = NULL;
  }

  if (avi->poly.set->id == SICONOS_SET_POLYHEDRON) {
    free_polyhedron(avi->poly.split);
    avi->poly.split = NULL;
  } else if (avi->poly.set->id == SICONOS_SET_POLYHEDRON_UNIFIED) {
    free_polyhedron_unified(avi->poly.unif);
    avi->poly.unif = NULL;
  }

  if (avi->q) {
    free(avi->q);
    avi->q = NULL;
  }
  if (avi->d) {
    free(avi->d);
    avi->d = NULL;
  }
  if (avi->lb) {
    free(avi->lb);
    avi->lb = NULL;
  }
  if (avi->ub) {
    free(avi->ub);
    avi->ub = NULL;
  }

  free(avi);
}

AffineVariationalInequalities* newAVI(void) {
  AffineVariationalInequalities* avi;
  avi = (AffineVariationalInequalities*)malloc(sizeof(AffineVariationalInequalities));
  avi->size = 0;
  avi->M = NULL;
  avi->q = NULL;
  avi->d = NULL;
  avi->poly.set = NULL;
  avi->lb = NULL;
  avi->ub = NULL;
  avi->cones = NULL;

  return avi;
}

#endif
