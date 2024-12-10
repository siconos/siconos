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
#include "MohrCoulomb2DProblem.h"

#include <assert.h>     // for assert
#include <math.h>       // for fabs
#include <stdio.h>      // for printf, fprintf, fscanf, NULL, fclose
#include <stdlib.h>     // for malloc, free, exit, EXIT_FAILURE
#include <string.h>     // for memcpy
#include <sys/errno.h>  // for errno

#include "NumericsMatrix.h"     // for NumericsMatrix, NM_create, RawNumeric...
#include "SiconosBlas.h"        // for cblas_dscal
#include "SparseBlockMatrix.h"  // for SBM_extract_component_3x3
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "io_tools.h"
#include "numerics_verbose.h"  // for CHECK_IO, numerics_error, numerics_pr...
#include "siconos_debug.h"     // for DEBUG_PRINT, DEBUG_PRINTF
#if defined(WITH_FCLIB)
#include "fclib_interface.h"
#endif

void mohrCoulomb2D_display(MohrCoulomb2DProblem* problem) {
  assert(problem);
  int n = problem->dimension * problem->numberOfCones;
  printf("MohrCoulomb2D Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfCones:%d \n", problem->numberOfCones);

  if (problem->M) {
    printf("M matrix:\n");
    NM_display(problem->M);
  } else
    printf("No M matrix:\n");

  if (problem->q) {
    printf("q vector:\n");
    NM_vector_display(problem->q, n);
  } else
    printf("No q vector:\n");

  if (problem->eta) {
    printf("mu vector:\n");
    NM_vector_display(problem->eta, problem->numberOfCones);
  } else
    printf("No eta vector:\n");

  if (problem->mu) {
    printf("mu vector:\n");
    NM_vector_display(problem->mu, problem->numberOfCones);
  } else
    printf("No mu vector:\n");
}

int mohrCoulomb2D_printInFile(MohrCoulomb2DProblem* problem, FILE* file) {
  if (!problem) {
    fprintf(stderr, "Numerics, MohrCoulomb2DProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;

  int d = problem->dimension;
  fprintf(file, "%d\n", d);
  int nc = problem->numberOfCones;
  fprintf(file, "%d\n", nc);
  NM_write_in_file(problem->M, file);
  for (i = 0; i < problem->M->size1; i++) {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < nc; i++) {
    fprintf(file, "%32.24e ", problem->eta[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < nc; i++) {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  return 0;
}

int mohrCoulomb2D_printInFilename(MohrCoulomb2DProblem* problem, char* filename) {
  int info = 0;
  FILE* file = fopen(filename, "w");

  if (!file) {
    return errno;
  }

  info = mohrCoulomb2D_printInFile(problem, file);

  fclose(file);
  return info;
}

MohrCoulomb2DProblem* mohrCoulomb2D_newFromFile(FILE* file) {
  MohrCoulomb2DProblem* problem = mohrCoulomb2DProblem_new();
  assert(file);
  DEBUG_PRINT(
      "Start -- int mohrCoulomb2D_newFromFile(MohrCoulomb2DProblem* problem, FILE* "
      "file)\n");
  int nc = 0, d = 0;
  int i;
  CHECK_IO(fscanf(file, "%d\n", &d));
  problem->dimension = d;
  DEBUG_PRINTF("problem->dimension = %i \n", problem->dimension);
  CHECK_IO(fscanf(file, "%d\n", &nc));
  problem->numberOfCones = nc;
  problem->M = NM_new_from_file(file);

  problem->q = (double*)malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++) {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }
  problem->eta = (double*)malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++) {
    CHECK_IO(fscanf(file, "%lf ", &(problem->eta[i])));
  }

  problem->mu = (double*)malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++) {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu[i])));
  }
  DEBUG_PRINT(
      "End --  int mohrCoulomb2D_newFromFile(MohrCoulomb2DProblem* problem, FILE* "
      "file)\n");

  return problem;
}

MohrCoulomb2DProblem* mohrCoulomb2D_new_from_filename(const char* filename) {
  MohrCoulomb2DProblem* problem = NULL;
  int is_hdf5 = check_hdf5_file(filename);
  // if the input file is an hdf5 file, we try to read it with fclib interface function.
  if (is_hdf5) {
    /* #if defined(WITH_FCLIB) */
    /* #else */
    numerics_error("MohrCoulomb2DProblem",
                   "Try to read an hdf5 file, while fclib interface is not active. Recompile "
                   "Siconos with fclib.",
                   filename);
    //#endif
  } else {
    FILE* file = fopen(filename, "r");
    if (!file) numerics_error("MohrCoulomb2DProblem", "Can not open file ", filename);

    problem = mohrCoulomb2D_newFromFile(file);
    fclose(file);
  }
  return problem;
}

void mohrCoulomb2DProblem_free(MohrCoulomb2DProblem* problem) {
  assert(problem);
  if (problem->M) {
    NM_clear(problem->M);
    free(problem->M);
    problem->M = NULL;
  }

  if (problem->mu) {
    free(problem->mu);
    problem->mu = NULL;
  }
  if (problem->eta) {
    free(problem->eta);
    problem->eta = NULL;
  }

  if (problem->q) {
    free(problem->q);
    problem->q = NULL;
  }

  free(problem);
}

MohrCoulomb2DProblem* mohrCoulomb2DProblem_new(void) {
  MohrCoulomb2DProblem* mc2d = (MohrCoulomb2DProblem*)malloc(sizeof(MohrCoulomb2DProblem));
  mc2d->dimension = 0;
  mc2d->numberOfCones = 0;
  mc2d->M = NULL;
  mc2d->q = NULL;
  mc2d->eta = NULL;
  mc2d->mu = NULL;

  return mc2d;
}

MohrCoulomb2DProblem* mohrCoulomb2DProblem_new_with_data(int dim, int nc, NumericsMatrix* M,
                                                         double* q, double* eta, double* mu) {
  MohrCoulomb2DProblem* mc2d = (MohrCoulomb2DProblem*)malloc(sizeof(MohrCoulomb2DProblem));

  mc2d->dimension = dim;
  mc2d->numberOfCones = nc;
  mc2d->M = M;
  mc2d->q = q;
  mc2d->eta = eta;
  mc2d->mu = mu;

  return mc2d;
}

MohrCoulomb2DProblem* mohrCoulomb2D_copy(MohrCoulomb2DProblem* problem) {
  assert(problem);

  int nc = problem->numberOfCones;
  int n = problem->M->size0;
  MohrCoulomb2DProblem* new = (MohrCoulomb2DProblem*)malloc(sizeof(MohrCoulomb2DProblem));
  new->dimension = problem->dimension;
  new->numberOfCones = problem->numberOfCones;
  new->M = NM_new();
  NM_copy(problem->M, new->M);
  new->q = (double*)malloc(n * sizeof(double));
  memcpy(new->q, problem->q, n * sizeof(double));
  new->eta = (double*)malloc(nc * sizeof(double));
  memcpy(new->eta, problem->eta, nc * sizeof(double));
  new->mu = (double*)malloc(nc * sizeof(double));
  memcpy(new->mu, problem->mu, nc * sizeof(double));
  return new;
}

void mohrCoulomb2D_rescaling(MohrCoulomb2DProblem* problem, double alpha, double gamma) {
  int n = problem->M->size0;
  /* scaling of M */
  NM_scal(alpha * gamma * gamma, problem->M);
  /* scaling of q */
  cblas_dscal(n, alpha * gamma, problem->q, 1);
}
