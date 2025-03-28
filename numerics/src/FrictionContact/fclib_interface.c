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
#include "CSparseMatrix.h"
#include "SiconosConfig.h"  // for WITH_FCLIB // IWYU pragma: keep

#ifdef WITH_FCLIB
#define DEBUG_NOCOLOR
#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include <assert.h>  // for assert
#include <fclib.h>   // for fclib_matrix, fclib_global, fclib_global_rolling, fclib_read...
#include <stdio.h>   // for NULL, fprintf, stderr
#include <stdlib.h>  // for malloc, free, exit, EXIT_F...

#include "CSparseMatrix.h"                        // for CSparseMatrix, CS_INT, cs_...
#include "FrictionContactProblem.h"               // for FrictionContactProblem
#include "GlobalFrictionContactProblem.h"         // for GlobalFrictionContactProblem
#include "GlobalRollingFrictionContactProblem.h"  // for GlobalRollingFrictionContactProblem
#include "NumericsMatrix.h"                       // for NumericsMatrix, RawNumeric...
#include "NumericsSparseMatrix.h"                 // for NumericsSparseMatrix, NSM_CSC
#include "SiconosConfig.h"                        // for WITH_FCLIB
#include "SparseBlockMatrix.h"                    // for SBM_from_csparse, SBM_to_s...
#include "fclib_interface.h"
#include "siconos_debug.h"  // for DEBUG_PRINT, DEBUG_PRINTF
#include "timers_interf.h"  // for MAYBE_UNUSED

// avoid a conflict with old csparse.h in case fclib includes it
#define _CS_H

static void int_to_csi(int* o, CS_INT* d, unsigned int n) {
  for (unsigned int i = 0; i < n; ++i) {
    d[i] = (CS_INT)o[i];
  }
}

static void csi_to_int(CS_INT* o, int* d, unsigned int n) {
  for (unsigned int i = 0; i < n; ++i) {
    d[i] = (int)o[i];
  }
}

FrictionContactProblem* from_fclib_local(const fclib_local* fclib_problem) {
  FrictionContactProblem* problem;

  problem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));

  problem->dimension = fclib_problem->spacedim;
  problem->mu = fclib_problem->mu;
  problem->q = fclib_problem->q;

  problem->numberOfContacts =
      fclib_problem->W->m / fclib_problem->spacedim; /* cf fclib spec */

  problem->M = NM_create(NM_SPARSE_BLOCK, fclib_problem->W->m, fclib_problem->W->n);

  problem->M->matrix1->block = NULL;
  problem->M->matrix1->index1_data = NULL;
  problem->M->matrix1->index2_data = NULL;

  CSparseMatrix W;

  W.nzmax = (CS_INT)fclib_problem->W->nzmax;
  W.m = (CS_INT)fclib_problem->W->m;
  W.n = (CS_INT)fclib_problem->W->n;

  if (fclib_problem->W->nz == -1) {
    /* compressed colums */
    W.p = (CS_INT*)malloc(sizeof(CS_INT) * (W.n + 1));
    int_to_csi(fclib_problem->W->p, W.p, (unsigned)(W.n + 1));
  } else if (fclib_problem->W->nz == -2) {
    /* compressed rows */
    W.p = (CS_INT*)malloc(sizeof(CS_INT) * (W.m + 1));
    int_to_csi(fclib_problem->W->p, W.p, (unsigned)(W.m + 1));
  } else {
    /* triplet */
    W.p = (CS_INT*)malloc(sizeof(CS_INT) * W.nzmax);
    int_to_csi(fclib_problem->W->p, W.p, (unsigned)W.nzmax);
  }

  W.i = (CS_INT*)malloc(sizeof(CS_INT) * W.nzmax);
  int_to_csi(fclib_problem->W->i, W.i, (unsigned)W.nzmax);

  W.x = fclib_problem->W->x;

  W.nz = fclib_problem->W->nz;

  SBM_from_csparse(problem->dimension, &W, problem->M->matrix1);

  free(W.p);
  free(W.i);

  NM_reset_versions(problem->M);

  return problem;
}
FrictionContactProblem* from_fclib_local_sparse(const fclib_local* fclib_problem) {
  FrictionContactProblem* problem;

  problem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));

  problem->dimension = fclib_problem->spacedim;
  problem->mu = fclib_problem->mu;
  problem->q = fclib_problem->q;

  problem->numberOfContacts =
      fclib_problem->W->m / fclib_problem->spacedim; /* cf fclib spec */

  problem->M = NM_create(NM_SPARSE, fclib_problem->W->m, fclib_problem->W->n);

  CSparseMatrix* W = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));

  W->nzmax = (CS_INT)fclib_problem->W->nzmax;
  W->m = (CS_INT)fclib_problem->W->m;
  W->n = (CS_INT)fclib_problem->W->n;

  if (fclib_problem->W->nz == -1) {
    /* compressed colums */
    W->p = (CS_INT*)malloc(sizeof(CS_INT) * (W->n + 1));
    int_to_csi(fclib_problem->W->p, W->p, (unsigned)(W->n + 1));
    problem->M->matrix2->csc = W;
    problem->M->matrix2->origin = NSM_CSC;
  } else if (fclib_problem->W->nz == -2) {
    /* compressed rows */
    W->p = (CS_INT*)malloc(sizeof(CS_INT) * (W->m + 1));
    int_to_csi(fclib_problem->W->p, W->p, (unsigned)(W->m + 1));
    problem->M->matrix2->csr = W;
    problem->M->matrix2->origin = NSM_CSR;
  } else {
    /* triplet */
    W->p = (CS_INT*)malloc(sizeof(CS_INT) * W->nzmax);
    int_to_csi(fclib_problem->W->p, W->p, (unsigned)W->nzmax);
    problem->M->matrix2->triplet = W;
    problem->M->matrix2->origin = NSM_TRIPLET;
  }

  W->i = (CS_INT*)malloc(sizeof(CS_INT) * W->nzmax);
  int_to_csi(fclib_problem->W->i, W->i, (unsigned)W->nzmax);

  W->x = fclib_problem->W->x;
  W->nz = fclib_problem->W->nz;

  NM_reset_versions(problem->M);

  return problem;
}

FrictionContactProblem* frictionContact_fclib_read(const char* path) {
  fclib_local* fclib_problem;

  fclib_problem = fclib_read_local(path);

  if (!fclib_problem) {
    return NULL;
  }

  return from_fclib_local(fclib_problem);
  // return from_fclib_local_sparse(fclib_problem);
}

int frictionContact_fclib_write_csr(FrictionContactProblem* problem, char* title,
                                    char* description, char* mathInfo, const char* path,
                                    int ndof) {
  int info = 0;

  fclib_local* fclib_problem;

  fclib_problem = (fclib_local*)malloc(sizeof(fclib_local));

  fclib_problem->spacedim = problem->dimension;
  fclib_problem->mu = problem->mu;
  fclib_problem->q = problem->q;

  fclib_problem->s = NULL;

  fclib_problem->info = (struct fclib_info*)malloc(sizeof(struct fclib_info));
  fclib_problem->info->title = title;
  fclib_problem->info->description = description;
  fclib_problem->info->math_info = mathInfo;

  fclib_problem->W = (struct fclib_matrix*)malloc(sizeof(struct fclib_matrix));
  fclib_problem->R = NULL;
  fclib_problem->V = NULL;

  fclib_problem->W->m = problem->M->size0;
  fclib_problem->W->n = problem->M->size1;

  CSparseMatrix* spmat = NULL;

  if (problem->M->storageType == NM_DENSE) /* Dense Matrix */
  {
    /* DEBUG_PRINT("NM_DENSE case\n"); */
    fclib_problem->W->nzmax = problem->M->size0 * problem->M->size1;
    fclib_problem->W->p = (int*)malloc((fclib_problem->W->m + 1) * sizeof(int));
    fclib_problem->W->i = (int*)malloc((fclib_problem->W->nzmax) * sizeof(int));
    fclib_problem->W->x = (double*)malloc((fclib_problem->W->nzmax) * sizeof(double));
    fclib_problem->W->nz = -2;
    fclib_problem->W->info = NULL;
    for (int i = 0; i < problem->M->size0; i++) {
      fclib_problem->W->p[i] = i * problem->M->size1;
      for (int j = 0; j < problem->M->size1; j++) {
        fclib_problem->W->x[i * problem->M->size1 + j] =
            problem->M->matrix0[j * problem->M->size0 + i];
        fclib_problem->W->i[i * problem->M->size1 + j] = j;
      }
    }
    fclib_problem->W->p[fclib_problem->W->m] = (fclib_problem->W->m) * problem->M->size1;

  } else if (problem->M->storageType == NM_SPARSE_BLOCK) /* Sparse block storage */
  {
    /* DEBUG_PRINT("NM_SPARSE_BLOCK case\n"); */
    spmat = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));
    int MAYBE_UNUSED res = SBM_to_sparse_init_memory(problem->M->matrix1, spmat);
    res = SBM_to_sparse(problem->M->matrix1, spmat);
    fclib_problem->W->nzmax = (int)spmat->nzmax;
    fclib_problem->W->m = (int)spmat->m;
    fclib_problem->W->n = (int)spmat->n;
    fclib_problem->W->x = spmat->x;
    fclib_problem->W->nz = (int)spmat->nz;

    if (spmat->nz == -1) {
      fclib_problem->W->p = (int*)malloc(sizeof(int) * (spmat->n + 1));
      csi_to_int(spmat->p, fclib_problem->W->p, (unsigned)(spmat->n + 1));
    } else if (spmat->nz == -2) {
      fclib_problem->W->p = (int*)malloc(sizeof(int) * (spmat->m + 1));
      csi_to_int(spmat->p, fclib_problem->W->p, (unsigned)(spmat->m + 1));
    } else {
      fclib_problem->W->p = (int*)malloc(sizeof(int) * spmat->nzmax);
      csi_to_int(spmat->p, fclib_problem->W->p, (unsigned)spmat->nzmax);
    }

    fclib_problem->W->i = (int*)malloc(sizeof(int) * spmat->nzmax);
    csi_to_int(spmat->i, fclib_problem->W->i, (unsigned)spmat->nzmax);

    fclib_problem->W->info = NULL;
  } else {
    fprintf(stderr, "frictionContact_fclib_write, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
    ;
  }

  info = fclib_write_local(fclib_problem, path);

  info = fclib_create_int_attributes_in_info(path, "numberOfDegreeOfFreedom", ndof);

  /*   fclib_delete_local (fclib_problem); */

  if (problem->M->storageType == NM_DENSE) /* Dense Matrix */
  {
    free(fclib_problem->W->x);
  } else if (problem->M->storageType == NM_SPARSE_BLOCK) {
    cs_spfree(spmat);
  }
  free(fclib_problem->W->p);
  free(fclib_problem->W->i);
  free(fclib_problem->W);
  free(fclib_problem->info);

  free(fclib_problem);

  return info;
}
int frictionContact_fclib_write(FrictionContactProblem* problem, char* title,
                                char* description, char* mathInfo, const char* path,
                                int ndof) {
  int info = 0;

  fclib_local* fclib_problem;

  fclib_problem = (fclib_local*)malloc(sizeof(fclib_local));

  fclib_problem->spacedim = problem->dimension;
  fclib_problem->mu = problem->mu;
  fclib_problem->q = problem->q;

  fclib_problem->s = NULL;

  fclib_problem->info = (struct fclib_info*)malloc(sizeof(struct fclib_info));
  fclib_problem->info->title = title;
  fclib_problem->info->description = description;
  fclib_problem->info->math_info = mathInfo;

  fclib_problem->R = NULL;
  fclib_problem->V = NULL;

  CSparseMatrix* spmat = NM_triplet(problem->M);

  fclib_problem->W = (struct fclib_matrix*)malloc(sizeof(struct fclib_matrix));
  fclib_problem->W->n = (int)spmat->n;
  fclib_problem->W->m = (int)spmat->m;

  /* We output only values up to nz and we set nzmax to nz */
  /* There is no interest to save in fclib file values from nz to nzmax */
  fclib_problem->W->nzmax = (int)spmat->nz;
  fclib_problem->W->p = (int*)malloc(sizeof(int) * (spmat->nz));
  csi_to_int(spmat->p, fclib_problem->W->p, (unsigned)spmat->nz);
  fclib_problem->W->i = (int*)malloc(sizeof(int) * (spmat->nz));
  csi_to_int(spmat->i, fclib_problem->W->i, (unsigned)spmat->nz);
  fclib_problem->W->x = spmat->x;
  fclib_problem->W->nz = (int)spmat->nz;
  fclib_problem->W->info = NULL;

  info = fclib_write_local(fclib_problem, path);

  info = fclib_create_int_attributes_in_info(path, "numberOfDegreeOfFreedom", ndof);

  free(fclib_problem->W->p);
  free(fclib_problem->W->i);
  free(fclib_problem->W);
  free(fclib_problem->info);

  free(fclib_problem);

  return info;
}
int frictionContact_fclib_write_guess(double* reaction, double* velocity, const char* path) {
  int info = 0;
  int number_of_guesses = 1;
  fclib_solution* guesses =
      (fclib_solution*)malloc(number_of_guesses * sizeof(fclib_solution));
  guesses->v = NULL;
  guesses->l = NULL;
  guesses->u = velocity;
  guesses->r = reaction;

  info = fclib_write_guesses(number_of_guesses, guesses, path);
  return info;
}

GlobalFrictionContactProblem* from_fclib_global(const fclib_global* fclib_problem) {
  GlobalFrictionContactProblem* problem = globalFrictionContactProblem_new();

  problem->dimension = fclib_problem->spacedim;
  problem->mu = fclib_problem->mu;
  problem->q = fclib_problem->f;
  problem->b = fclib_problem->w;
  problem->env = NULL;

  problem->numberOfContacts =
      fclib_problem->H->n / fclib_problem->spacedim; /* cf fclib spec */

  problem->M = NM_create(NM_SPARSE, fclib_problem->M->m, fclib_problem->M->n);

  CSparseMatrix* M = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));
  M->nzmax = (CS_INT)fclib_problem->M->nzmax;
  M->m = (CS_INT)fclib_problem->M->m;
  M->n = (CS_INT)fclib_problem->M->n;

  M->x = fclib_problem->M->x;

  if (fclib_problem->M->nz == -1) {
    /* compressed colums */
    problem->M->matrix2->csc = M;
    problem->M->matrix2->origin = NSM_CSC;
    M->nz = (CS_INT)fclib_problem->M->nz;
    M->p = (CS_INT*)malloc(sizeof(CS_INT) * (M->n + 1));
    int_to_csi(fclib_problem->M->p, M->p, (unsigned)(M->n + 1));
  } else if (fclib_problem->M->nz == -2) {
    /* compressed rows */
    M->nz = (CS_INT)fclib_problem->M->nz;
    M->p = (CS_INT*)malloc(sizeof(CS_INT) * (M->m + 1));
    int_to_csi(fclib_problem->M->p, M->p, (unsigned)(M->m + 1));
    /* since  problem->M->matrix2->csr does not exist, we need
       to fill transform M into a triplet or csc before returning
     */

    fprintf(stderr, "from_fclib_global not implemented for csr matrices.\n");
    exit(EXIT_FAILURE);
    ;
  } else {
    /* triplet */
    problem->M->matrix2->triplet = M;
    problem->M->matrix2->origin = NSM_TRIPLET;
    M->nz = (CS_INT)fclib_problem->M->nz;
    M->p = (CS_INT*)malloc(sizeof(CS_INT) * M->nzmax);
    int_to_csi(fclib_problem->M->p, M->p, (unsigned)M->nz);
  }
  M->i = (CS_INT*)malloc(sizeof(CS_INT) * M->nzmax);
  int_to_csi(fclib_problem->M->i, M->i, (unsigned)M->nz);

  problem->H = NM_create(NM_SPARSE, fclib_problem->H->m, fclib_problem->H->n);

  CSparseMatrix* H = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));

  H->nzmax = (CS_INT)fclib_problem->H->nzmax;
  H->m = (CS_INT)fclib_problem->H->m;
  H->n = (CS_INT)fclib_problem->H->n;
  H->nz = (CS_INT)fclib_problem->H->nz;
  H->x = fclib_problem->H->x;

  if (fclib_problem->H->nz == -1) {
    /* compressed colums */
    problem->H->matrix2->csc = H;
    problem->H->matrix2->origin = NSM_CSC;
    H->p = (CS_INT*)malloc(sizeof(CS_INT) * (H->n + 1));
    int_to_csi(fclib_problem->H->p, H->p, (unsigned)(H->n + 1));
  } else if (fclib_problem->H->nz == -2) {
    /* compressed rows */
    fprintf(stderr, "from_fclib_global not implemented for csr matrices.\n");
    exit(EXIT_FAILURE);
    ;
  } else {
    /* triplet */
    problem->H->matrix2->triplet = H;
    problem->H->matrix2->origin = NSM_TRIPLET;
    H->p = (CS_INT*)malloc(sizeof(CS_INT) * H->nzmax);
    int_to_csi(fclib_problem->H->p, H->p, (unsigned)H->nz);
  }

  H->i = (CS_INT*)malloc(sizeof(CS_INT) * H->nzmax);
  int_to_csi(fclib_problem->H->i, H->i, (unsigned)H->nz);

  NM_reset_versions(problem->M);
  return problem;
}

GlobalFrictionContactProblem* globalFrictionContact_fclib_read(const char* path) {
  fclib_global* fclib_problem;

  fclib_problem = fclib_read_global(path);

  if (!fclib_problem) {
    return NULL;
  }

  GlobalFrictionContactProblem* gfc3d = from_fclib_global(fclib_problem);

  // gfc3d problem points to the following chunks of memory
  // allocated in fclib_read_global
  // we set them to NULL before calling fclib_delete_global to keep them
  fclib_problem->f = NULL;
  fclib_problem->w = NULL;
  fclib_problem->mu = NULL;
  fclib_problem->M->x = NULL;
  fclib_problem->H->x = NULL;
  if (fclib_problem->G) fclib_problem->G->x = NULL;

  fclib_delete_global(fclib_problem);
  free(fclib_problem);

  return gfc3d;
}

int globalFrictionContact_fclib_write(GlobalFrictionContactProblem* problem, char* title,
                                      char* description, char* mathInfo, const char* path) {
  int rinfo = 0;

  DEBUG_PRINTF("construction of fclib_problem in %s with title = %s and description = %s\n",
               path, title, description);
  if (problem->numberOfContacts == 0) {
    DEBUG_PRINT("zero contacts");
    return rinfo;
  }
  /* globalFrictionContact_display(problem); */
  /* FILE * file  =  fopen("toto.dat", "w"); */
  /* globalFrictionContact_printInFile(problem, file); */

  fclib_global* fclib_problem;
  fclib_problem = (fclib_global*)malloc(sizeof(fclib_global));

  fclib_problem->info = (struct fclib_info*)malloc(sizeof(struct fclib_info));

  fclib_problem->info->title = title;
  fclib_problem->info->description = description;
  fclib_problem->info->math_info = mathInfo;

  fclib_problem->spacedim = problem->dimension;
  fclib_problem->mu = problem->mu;
  fclib_problem->w = problem->b;
  fclib_problem->f = problem->q;

  /* only sparse storage */
  assert(problem->M->matrix2);
  assert(problem->H->matrix2);

  CSparseMatrix* spmat = NM_triplet(problem->M);

  fclib_problem->M = (struct fclib_matrix*)malloc(sizeof(struct fclib_matrix));
  fclib_problem->M->n = (int)spmat->n;
  fclib_problem->M->m = (int)spmat->m;

  /* We output only values up to nz and we set nzmax to nz */
  /* There is no interest to save in fclib file values from nz to nzmax */
  fclib_problem->M->nzmax = (int)spmat->nz;
  fclib_problem->M->p = (int*)malloc(sizeof(int) * (spmat->nz));
  csi_to_int(spmat->p, fclib_problem->M->p, (unsigned)spmat->nz);
  fclib_problem->M->i = (int*)malloc(sizeof(int) * (spmat->nz));
  csi_to_int(spmat->i, fclib_problem->M->i, (unsigned)spmat->nz);
  fclib_problem->M->x = spmat->x;
  fclib_problem->M->nz = (int)spmat->nz;
  fclib_problem->M->info = NULL;

  spmat = NM_triplet(problem->H);

  fclib_problem->H = (struct fclib_matrix*)malloc(sizeof(struct fclib_matrix));
  fclib_problem->H->n = (int)spmat->n;
  fclib_problem->H->m = (int)spmat->m;
  fclib_problem->H->nzmax = (int)spmat->nz;
  fclib_problem->H->p = (int*)malloc(sizeof(int) * spmat->nz);
  csi_to_int(spmat->p, fclib_problem->H->p, (unsigned)spmat->nz);
  fclib_problem->H->i = (int*)malloc(sizeof(int) * spmat->nz);
  csi_to_int(spmat->i, fclib_problem->H->i, (unsigned)spmat->nz);
  fclib_problem->H->x = spmat->x;
  fclib_problem->H->nz = (int)spmat->nz;
  fclib_problem->H->info = NULL;

  fclib_problem->G = NULL;
  fclib_problem->b = NULL;

  DEBUG_PRINT("write in fclib of fclib_problem\n");

  rinfo = fclib_write_global(fclib_problem, path);
  DEBUG_PRINT("end of write in fclib of fclib_problem\n");

  free(fclib_problem->M->p);
  free(fclib_problem->M->i);
  free(fclib_problem->H->p);
  free(fclib_problem->H->i);
  free(fclib_problem->H);
  free(fclib_problem->M);
  free(fclib_problem->info);
  free(fclib_problem);

  return rinfo;
}

GlobalRollingFrictionContactProblem* from_fclib_global_rolling(
    const fclib_global_rolling* fclib_problem) {
  GlobalRollingFrictionContactProblem* problem = globalRollingFrictionContactProblem_new();

  problem->dimension = fclib_problem->spacedim;
  problem->mu = fclib_problem->mu;
  problem->mu_r = fclib_problem->mu_r;
  problem->q = fclib_problem->f;
  problem->b = fclib_problem->w;

  problem->numberOfContacts =
      fclib_problem->H->n / fclib_problem->spacedim; /* cf fclib spec */

  problem->M = NM_create(NM_SPARSE, fclib_problem->M->m, fclib_problem->M->n);

  CSparseMatrix* M = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));
  M->nzmax = (CS_INT)fclib_problem->M->nzmax;
  M->m = (CS_INT)fclib_problem->M->m;
  M->n = (CS_INT)fclib_problem->M->n;

  M->x = fclib_problem->M->x;

  if (fclib_problem->M->nz == -1) {
    /* compressed colums */
    problem->M->matrix2->csc = M;
    problem->M->matrix2->origin = NSM_CSC;
    M->nz = (CS_INT)fclib_problem->M->nz;
    M->p = (CS_INT*)malloc(sizeof(CS_INT) * (M->n + 1));
    int_to_csi(fclib_problem->M->p, M->p, (unsigned)(M->n + 1));
  } else if (fclib_problem->M->nz == -2) {
    /* compressed rows */
    M->nz = (CS_INT)fclib_problem->M->nz;
    M->p = (CS_INT*)malloc(sizeof(CS_INT) * (M->m + 1));
    int_to_csi(fclib_problem->M->p, M->p, (unsigned)(M->m + 1));
    /* since  problem->M->matrix2->csr does not exist, we need
       to fill transform M into a triplet or csc before returning
     */

    fprintf(stderr, "from_fclib_local not implemented for csr matrices.\n");
    exit(EXIT_FAILURE);
    ;
  } else {
    /* triplet */
    problem->M->matrix2->triplet = M;
    problem->M->matrix2->origin = NSM_TRIPLET;
    M->nz = (CS_INT)fclib_problem->M->nz;
    M->p = (CS_INT*)malloc(sizeof(CS_INT) * M->nzmax);
    int_to_csi(fclib_problem->M->p, M->p, (unsigned)M->nz);
  }
  M->i = (CS_INT*)malloc(sizeof(CS_INT) * M->nzmax);
  int_to_csi(fclib_problem->M->i, M->i, (unsigned)M->nz);

  problem->H = NM_create(NM_SPARSE, fclib_problem->H->m, fclib_problem->H->n);

  CSparseMatrix* H = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));

  H->nzmax = (CS_INT)fclib_problem->H->nzmax;
  H->m = (CS_INT)fclib_problem->H->m;
  H->n = (CS_INT)fclib_problem->H->n;
  H->nz = (CS_INT)fclib_problem->H->nz;
  H->x = fclib_problem->H->x;

  if (fclib_problem->H->nz == -1) {
    /* compressed colums */
    problem->H->matrix2->csc = H;
    problem->H->matrix2->origin = NSM_CSC;
    H->p = (CS_INT*)malloc(sizeof(CS_INT) * (H->n + 1));
    int_to_csi(fclib_problem->H->p, H->p, (unsigned)(H->n + 1));
  } else if (fclib_problem->H->nz == -2) {
    /* compressed rows */
    fprintf(stderr, "from_fclib_local not implemented for csr matrices.\n");
    exit(EXIT_FAILURE);
    ;
  } else {
    /* triplet */
    problem->H->matrix2->triplet = H;
    problem->H->matrix2->origin = NSM_TRIPLET;
    H->p = (CS_INT*)malloc(sizeof(CS_INT) * H->nzmax);
    int_to_csi(fclib_problem->H->p, H->p, (unsigned)H->nz);
  }

  H->i = (CS_INT*)malloc(sizeof(CS_INT) * H->nzmax);
  int_to_csi(fclib_problem->H->i, H->i, (unsigned)H->nz);

  NM_reset_versions(problem->M);
  return problem;
}

GlobalRollingFrictionContactProblem* globalRollingFrictionContact_fclib_read(
    const char* path) {
  fclib_global_rolling* fclib_problem;
  fclib_problem = fclib_read_global_rolling(path);

  if (!fclib_problem) {
    return NULL;
  }

  return from_fclib_global_rolling(fclib_problem);
}

int globalRollingFrictionContact_fclib_write(GlobalRollingFrictionContactProblem* problem,
                                             char* title, char* description, char* mathInfo,
                                             const char* path) {
  int rinfo = 0;

  DEBUG_PRINTF("construction of fclib_problem in %s with title = %s and description = %s\n",
               path, title, description);
  if (problem->numberOfContacts == 0) {
    DEBUG_PRINT("zero contacts");
    return rinfo;
  }
  /* globalFrictionContact_display(problem); */
  /* FILE * file  =  fopen("toto.dat", "w"); */
  /* globalFrictionContact_printInFile(problem, file); */

  fclib_global_rolling* fclib_problem;
  fclib_problem = (fclib_global_rolling*)malloc(sizeof(fclib_global_rolling));

  fclib_problem->info = (struct fclib_info*)malloc(sizeof(struct fclib_info));

  fclib_problem->info->title = title;
  fclib_problem->info->description = description;
  fclib_problem->info->math_info = mathInfo;

  fclib_problem->spacedim = problem->dimension;
  fclib_problem->mu = problem->mu;
  fclib_problem->mu_r = problem->mu_r;
  fclib_problem->w = problem->b;
  fclib_problem->f = problem->q;

  /* only sparse storage */
  assert(problem->M->matrix2);
  assert(problem->H->matrix2);

  CSparseMatrix* spmat = NM_triplet(problem->M);

  fclib_problem->M = (struct fclib_matrix*)malloc(sizeof(struct fclib_matrix));
  fclib_problem->M->n = (int)spmat->n;
  fclib_problem->M->m = (int)spmat->m;

  /* we output only values up to nz and we set nzmax to nz */
  /* There is no interest to save in fclib file values from nz to nzmax */
  fclib_problem->M->nzmax = (int)spmat->nz;

  fclib_problem->M->p = (int*)malloc(sizeof(int) * (spmat->nz));
  csi_to_int(spmat->p, fclib_problem->M->p, (unsigned)spmat->nz);
  fclib_problem->M->i = (int*)malloc(sizeof(int) * (spmat->nz));
  csi_to_int(spmat->i, fclib_problem->M->i, (unsigned)spmat->nz);
  fclib_problem->M->x = spmat->x;
  fclib_problem->M->nz = (int)spmat->nz;
  fclib_problem->M->info = NULL;

  spmat = NM_triplet(problem->H);

  fclib_problem->H = (struct fclib_matrix*)malloc(sizeof(struct fclib_matrix));
  fclib_problem->H->n = (int)spmat->n;
  fclib_problem->H->m = (int)spmat->m;
  fclib_problem->H->nzmax = (int)spmat->nz;
  fclib_problem->H->p = (int*)malloc(sizeof(int) * spmat->nz);
  csi_to_int(spmat->p, fclib_problem->H->p, (unsigned)spmat->nz);
  fclib_problem->H->i = (int*)malloc(sizeof(int) * spmat->nz);
  csi_to_int(spmat->i, fclib_problem->H->i, (unsigned)spmat->nz);
  fclib_problem->H->x = spmat->x;
  fclib_problem->H->nz = (int)spmat->nz;
  fclib_problem->H->info = NULL;

  fclib_problem->G = NULL;
  fclib_problem->b = NULL;

  DEBUG_PRINT("write in fclib of fclib_problem\n");

  rinfo = fclib_write_global_rolling(fclib_problem, path);
  DEBUG_PRINT("end of write in fclib of fclib_problem\n");

  free(fclib_problem->M->p);
  free(fclib_problem->M->i);
  free(fclib_problem->H->p);
  free(fclib_problem->H->i);
  free(fclib_problem->H);
  free(fclib_problem->M);
  free(fclib_problem->info);
  free(fclib_problem);

  return rinfo;
}

#endif
