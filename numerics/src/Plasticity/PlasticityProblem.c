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
#include "PlasticityProblem.h"

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
#include "numerics_verbose.h"  // for check_io, numerics_error, numerics_pr...
#include "numerics_errors.h"
#include "siconos_debug.h"     // for DEBUG_PRINT, DEBUG_PRINTF
#if defined(WITH_FCLIB)
#include "fclib_interface.h"
#endif

/* Helper macros for accessing model parameters */
#define GET_ETA(problem) ((problem)->model ? (problem)->model->eta : NULL)
#define GET_THETA(problem) ((problem)->model ? (problem)->model->theta : NULL)

Plasticity_DruckerPrager_model *plasticity_DruckerPrager_model_new(double *eta, double *theta) {
  Plasticity_DruckerPrager_model *model = (Plasticity_DruckerPrager_model *)malloc(
      sizeof(Plasticity_DruckerPrager_model));
  model->eta = eta;
  model->theta = theta;
  return model;
}

void plasticity_DruckerPrager_model_free(Plasticity_DruckerPrager_model *model) {
  if (model) {
    free(model);
  }
}

void plasticity_display(PlasticityProblem *problem) {
  assert(problem);
  int n = problem->dimension * problem->numberOfCones;
  printf("PlasticityProblem Display :\n-------------\n");
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

  double *eta = GET_ETA(problem);
  if (eta) {
    printf("eta vector:\n");
    NM_vector_display(eta, problem->numberOfCones);
  } else
    printf("No eta vector:\n");

  double *theta = GET_THETA(problem);
  if (theta) {
    printf("theta vector:\n");
    NM_vector_display(theta, problem->numberOfCones);
  } else
    printf("No theta vector:\n");
}

int plasticity_printInFile(PlasticityProblem *problem, FILE *file) {
  if (!problem) {
    CHECK_ARG(0, "Numerics, PlasticityProblem printInFile failed, NULL input.\n");
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
  
  double *eta = GET_ETA(problem);
  if (eta) {
    for (i = 0; i < nc; i++) {
      fprintf(file, "%32.24e ", eta[i]);
    }
    fprintf(file, "\n");
  }
  
  double *theta = GET_THETA(problem);
  if (theta) {
    for (i = 0; i < nc; i++) {
      fprintf(file, "%32.24e ", theta[i]);
    }
    fprintf(file, "\n");
  }
  return 0;
}

int plasticity_printInFilename(PlasticityProblem *problem, char *filename) {
  int info = 0;
  FILE *file = fopen(filename, "w");

  if (!file) {
    return errno;
  }

  info = plasticity_printInFile(problem, file);

  fclose(file);
  return info;
}

PlasticityProblem *plasticity_newFromFile(FILE *file) {
  PlasticityProblem *problem = plasticityProblem_new();
  if (!file) return NULL;
  DEBUG_PRINT(
      "Start -- int plasticity_newFromFile(PlasticityProblem* problem, FILE* "
      "file)\n");
  int nc = 0, d = 0;
  int i;
  check_io(fscanf(file, "%d\n", &d));
  problem->dimension = d;
  DEBUG_PRINTF("problem->dimension = %i \n", problem->dimension);
  check_io(fscanf(file, "%d\n", &nc));
  problem->numberOfCones = nc;
  problem->M = NM_new_from_file(file);

  problem->q = (double *)malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++) {
    check_io(fscanf(file, "%lf ", &(problem->q[i])));
  }
  
  /* Allocate and read model parameters */
  Plasticity_DruckerPrager_model *model = (Plasticity_DruckerPrager_model *)malloc(
      sizeof(Plasticity_DruckerPrager_model));
  model->eta = (double *)malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++) {
    check_io(fscanf(file, "%lf ", &(model->eta[i])));
  }

  model->theta = (double *)malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++) {
    check_io(fscanf(file, "%lf ", &(model->theta[i])));
  }
  problem->model = model;
  
  DEBUG_PRINT(
      "End --  int plasticity_newFromFile(PlasticityProblem* problem, FILE* "
      "file)\n");

  return problem;
}

PlasticityProblem *plasticity_new_from_filename(const char *filename) {
  PlasticityProblem *problem = NULL;
  int is_hdf5 = check_hdf5_file(filename);
  // if the input file is an hdf5 file, we try to read it with fclib interface function.
  if (is_hdf5) {
    /* #if defined(WITH_FCLIB) */
    /* #else */
    numerics_error("PlasticityProblem",
                   "Try to read an hdf5 file, while fclib interface is not active. Recompile "
                   "Siconos with fclib.",
                   filename);
    //#endif
  } else {
    FILE *file = fopen(filename, "r");
    if (!file) numerics_error("PlasticityProblem", "Can not open file ", filename);

    problem = plasticity_newFromFile(file);
    fclose(file);
  }
  return problem;
}

void plasticityProblem_free(PlasticityProblem *problem) {
  assert(problem);
  if (problem->M) {
    NM_clear(problem->M);
    free(problem->M);
    problem->M = NULL;
  }

  if (problem->model) {
    if (problem->model->theta) {
      free(problem->model->theta);
      problem->model->theta = NULL;
    }
    if (problem->model->eta) {
      free(problem->model->eta);
      problem->model->eta = NULL;
    }
    free(problem->model);
    problem->model = NULL;
  }

  if (problem->q) {
    free(problem->q);
    problem->q = NULL;
  }

  free(problem);
}

PlasticityProblem *plasticityProblem_new(void) {
  PlasticityProblem *problem = (PlasticityProblem *)malloc(sizeof(PlasticityProblem));
  problem->dimension = 0;
  problem->numberOfCones = 0;
  problem->M = NULL;
  problem->q = NULL;
  problem->model = NULL;

  return problem;
}

PlasticityProblem *plasticityProblem_new_with_data(int dim, int nc, NumericsMatrix *M,
                                                    double *q, Plasticity_DruckerPrager_model *model) {
  PlasticityProblem *problem = (PlasticityProblem *)malloc(sizeof(PlasticityProblem));

  problem->dimension = dim;
  problem->numberOfCones = nc;
  problem->M = M;
  problem->q = q;
  problem->model = model;

  return problem;
}

PlasticityProblem *plasticity_copy(PlasticityProblem *problem) {
  if (!problem) return NULL;

  int nc = problem->numberOfCones;
  int n = problem->M->size0;
  PlasticityProblem *new_problem = (PlasticityProblem *)malloc(sizeof(PlasticityProblem));
  new_problem->dimension = problem->dimension;
  new_problem->numberOfCones = problem->numberOfCones;
  new_problem->M = NM_new();
  NM_copy(problem->M, new_problem->M);
  new_problem->q = (double *)malloc(n * sizeof(double));
  memcpy(new_problem->q, problem->q, n * sizeof(double));
  
  if (problem->model) {
    new_problem->model = (Plasticity_DruckerPrager_model *)malloc(
        sizeof(Plasticity_DruckerPrager_model));
    new_problem->model->eta = (double *)malloc(nc * sizeof(double));
    memcpy(new_problem->model->eta, problem->model->eta, nc * sizeof(double));
    new_problem->model->theta = (double *)malloc(nc * sizeof(double));
    memcpy(new_problem->model->theta, problem->model->theta, nc * sizeof(double));
  } else {
    new_problem->model = NULL;
  }
  return new_problem;
}

void plasticity_rescaling(PlasticityProblem *problem, double alpha, double gamma) {
  int n = problem->M->size0;
  /* scaling of M */
  NM_scal(alpha * gamma * gamma, problem->M);
  /* scaling of q */
  cblas_dscal(n, alpha * gamma, problem->q, 1);
}

/* ============================================================================
 * Backward compatibility functions
 * ============================================================================ */

void plasticity2D_display(Plasticity2DProblem *problem) {
  plasticity_display(problem);
}

int plasticity2D_printInFile(Plasticity2DProblem *problem, FILE *file) {
  return plasticity_printInFile(problem, file);
}

int plasticity2D_printInFilename(Plasticity2DProblem *problem, char *filename) {
  return plasticity_printInFilename(problem, filename);
}

Plasticity2DProblem *plasticity2D_newFromFile(FILE *file) {
  return plasticity_newFromFile(file);
}

Plasticity2DProblem *plasticity2D_new_from_filename(const char *filename) {
  return plasticity_new_from_filename(filename);
}

void plasticity2DProblem_free(Plasticity2DProblem *problem) {
  plasticityProblem_free(problem);
}

Plasticity2DProblem *plasticity2DProblem_new(void) {
  return plasticityProblem_new();
}

Plasticity2DProblem *plasticity2DProblem_new_with_data(int dim, int nc, NumericsMatrix *M,
                                                         double *q, double *eta, double *theta) {
  /* Create model from eta and theta for backward compatibility */
  Plasticity_DruckerPrager_model *model = NULL;
  if (eta || theta) {
    model = plasticity_DruckerPrager_model_new(eta, theta);
  }
  return plasticityProblem_new_with_data(dim, nc, M, q, model);
}

Plasticity2DProblem *plasticity2D_copy(Plasticity2DProblem *problem) {
  return plasticity_copy(problem);
}

void plasticity2D_rescaling(Plasticity2DProblem *problem, double alpha, double gamma) {
  plasticity_rescaling(problem, alpha, gamma);
}
