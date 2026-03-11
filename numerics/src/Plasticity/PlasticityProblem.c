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

/* ============================================================================
 * Helper macros for accessing model parameters with type checking
 * ============================================================================ */

#define CHECK_MODEL_TYPE(problem, expected_type) \
  do { \
    if ((problem)->model_type != (expected_type)) { \
      numerics_error("PlasticityProblem", "Expected model type %s but got %s", \
                     plasticity_model_type_to_string(expected_type), \
                     plasticity_model_type_to_string((problem)->model_type)); \
    } \
  } while(0)

#define GET_DRUCKER_PRAGER(problem) \
  ((problem)->model.drucker_prager)

#define GET_VON_MISES(problem) \
  ((problem)->model.von_mises)

#define GET_ETA(problem) \
  (GET_DRUCKER_PRAGER(problem)->eta)

#define GET_THETA(problem) \
  (GET_DRUCKER_PRAGER(problem)->theta)

#define GET_SIGMA_Y(problem) \
  (GET_VON_MISES(problem)->sigma_y)

/* ============================================================================
 * Model Type Functions
 * ============================================================================ */

const char *plasticity_model_type_to_string(PlasticityModelType model_type) {
  switch (model_type) {
    case PLASTICITY_MODEL_DRUCKER_PRAGER:
      return "Drucker-Prager";
    case PLASTICITY_MODEL_VON_MISES:
      return "Von Mises";
    case PLASTICITY_MODEL_UNKNOWN:
      return "Unknown";
    default:
      return "Invalid";
  }
}

/* ============================================================================
 * Drucker-Prager Model Functions
 * ============================================================================ */

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

/* ============================================================================
 * Von Mises Model Functions
 * ============================================================================ */

Plasticity_VonMises_model *plasticity_VonMises_model_new(double *sigma_y) {
  Plasticity_VonMises_model *model = (Plasticity_VonMises_model *)malloc(
      sizeof(Plasticity_VonMises_model));
  model->sigma_y = sigma_y;
  return model;
}

void plasticity_VonMises_model_free(Plasticity_VonMises_model *model) {
  if (model) {
    free(model);
  }
}

/* ============================================================================
 * Problem Functions
 * ============================================================================ */

void plasticity_display(PlasticityProblem *problem) {
  assert(problem);
  int n = problem->dimension * problem->numberOfCones;
  printf("PlasticityProblem Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfCones:%d \n", problem->numberOfCones);
  printf("model_type: %s\n", plasticity_model_type_to_string(problem->model_type));

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

  /* Display model-specific parameters based on model_type */
  switch (problem->model_type) {
    case PLASTICITY_MODEL_DRUCKER_PRAGER: {
      double *eta = GET_ETA(problem);
      double *theta = GET_THETA(problem);
      if (eta) {
        printf("eta vector (Drucker-Prager):\n");
        NM_vector_display(eta, problem->numberOfCones);
      } else
        printf("No eta vector\n");

      if (theta) {
        printf("theta vector (Drucker-Prager):\n");
        NM_vector_display(theta, problem->numberOfCones);
      } else
        printf("No theta vector\n");
      break;
    }
    case PLASTICITY_MODEL_VON_MISES: {
      double *sigma_y = GET_SIGMA_Y(problem);
      if (sigma_y) {
        printf("sigma_y vector (Von Mises):\n");
        NM_vector_display(sigma_y, problem->numberOfCones);
      } else
        printf("No sigma_y vector\n");
      break;
    }
    case PLASTICITY_MODEL_UNKNOWN:
    default:
      printf("No model parameters (model_type = %d)\n", problem->model_type);
      break;
  }
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
  
  /* TODO: Add file format version marker for multiple model types
   * For now, we don't write model_type to maintain backward compatibility */
  
  NM_write_in_file(problem->M, file);
  for (i = 0; i < problem->M->size1; i++) {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  
  /* Write model-specific parameters */
  switch (problem->model_type) {
    case PLASTICITY_MODEL_DRUCKER_PRAGER: {
      double *eta = GET_ETA(problem);
      double *theta = GET_THETA(problem);
      if (eta) {
        for (i = 0; i < nc; i++) {
          fprintf(file, "%32.24e ", eta[i]);
        }
        fprintf(file, "\n");
      }
      if (theta) {
        for (i = 0; i < nc; i++) {
          fprintf(file, "%32.24e ", theta[i]);
        }
        fprintf(file, "\n");
      }
      break;
    }
    case PLASTICITY_MODEL_VON_MISES: {
      double *sigma_y = GET_SIGMA_Y(problem);
      if (sigma_y) {
        for (i = 0; i < nc; i++) {
          fprintf(file, "%32.24e ", sigma_y[i]);
        }
        fprintf(file, "\n");
      }
      break;
    }
    default:
      break;
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
  
  /* For now, default to Drucker-Prager model for file reading
   * TODO: Add file format version marker to support multiple model types in files */
  problem->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
  
  problem->M = NM_new_from_file(file);

  problem->q = (double *)malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++) {
    check_io(fscanf(file, "%lf ", &(problem->q[i])));
  }
  
  /* Read model-specific parameters based on model_type */
  switch (problem->model_type) {
    case PLASTICITY_MODEL_DRUCKER_PRAGER: {
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
      problem->model.drucker_prager = model;
      break;
    }
    case PLASTICITY_MODEL_VON_MISES: {
      Plasticity_VonMises_model *model = (Plasticity_VonMises_model *)malloc(
          sizeof(Plasticity_VonMises_model));
      model->sigma_y = (double *)malloc(nc * sizeof(double));
      for (i = 0; i < nc; i++) {
        check_io(fscanf(file, "%lf ", &(model->sigma_y[i])));
      }
      problem->model.von_mises = model;
      break;
    }
    default:
      problem->model.generic = NULL;
      break;
  }
  
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

  /* Free model-specific data based on model_type */
  switch (problem->model_type) {
    case PLASTICITY_MODEL_DRUCKER_PRAGER: {
      Plasticity_DruckerPrager_model *model = GET_DRUCKER_PRAGER(problem);
      if (model) {
        if (model->theta) {
          free(model->theta);
          model->theta = NULL;
        }
        if (model->eta) {
          free(model->eta);
          model->eta = NULL;
        }
        free(model);
        problem->model.drucker_prager = NULL;
      }
      break;
    }
    case PLASTICITY_MODEL_VON_MISES: {
      Plasticity_VonMises_model *model = GET_VON_MISES(problem);
      if (model) {
        if (model->sigma_y) {
          free(model->sigma_y);
          model->sigma_y = NULL;
        }
        free(model);
        problem->model.von_mises = NULL;
      }
      break;
    }
    default:
      if (problem->model.generic) {
        free(problem->model.generic);
        problem->model.generic = NULL;
      }
      break;
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
  problem->model_type = PLASTICITY_MODEL_UNKNOWN;
  problem->model.generic = NULL;

  return problem;
}

PlasticityProblem *plasticityProblem_new_with_data(int dim, int nc, NumericsMatrix *M,
                                                    double *q, Plasticity_DruckerPrager_model *model) {
  PlasticityProblem *problem = (PlasticityProblem *)malloc(sizeof(PlasticityProblem));

  problem->dimension = dim;
  problem->numberOfCones = nc;
  problem->M = M;
  problem->q = q;
  if (model) {
    problem->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
    problem->model.drucker_prager = model;
  } else {
    problem->model_type = PLASTICITY_MODEL_UNKNOWN;
    problem->model.generic = NULL;
  }

  return problem;
}

PlasticityProblem *plasticity_copy(PlasticityProblem *problem) {
  if (!problem) return NULL;

  int nc = problem->numberOfCones;
  int n = problem->M->size0;
  PlasticityProblem *new_problem = (PlasticityProblem *)malloc(sizeof(PlasticityProblem));
  new_problem->dimension = problem->dimension;
  new_problem->numberOfCones = problem->numberOfCones;
  new_problem->model_type = problem->model_type;
  new_problem->M = NM_new();
  NM_copy(problem->M, new_problem->M);
  new_problem->q = (double *)malloc(n * sizeof(double));
  memcpy(new_problem->q, problem->q, n * sizeof(double));
  
  /* Copy model-specific data based on model_type */
  switch (problem->model_type) {
    case PLASTICITY_MODEL_DRUCKER_PRAGER: {
      Plasticity_DruckerPrager_model *model = GET_DRUCKER_PRAGER(problem);
      if (model) {
        new_problem->model.drucker_prager = (Plasticity_DruckerPrager_model *)malloc(
            sizeof(Plasticity_DruckerPrager_model));
        new_problem->model.drucker_prager->eta = (double *)malloc(nc * sizeof(double));
        memcpy(new_problem->model.drucker_prager->eta, model->eta, nc * sizeof(double));
        new_problem->model.drucker_prager->theta = (double *)malloc(nc * sizeof(double));
        memcpy(new_problem->model.drucker_prager->theta, model->theta, nc * sizeof(double));
      } else {
        new_problem->model.drucker_prager = NULL;
      }
      break;
    }
    case PLASTICITY_MODEL_VON_MISES: {
      Plasticity_VonMises_model *model = GET_VON_MISES(problem);
      if (model) {
        new_problem->model.von_mises = (Plasticity_VonMises_model *)malloc(
            sizeof(Plasticity_VonMises_model));
        new_problem->model.von_mises->sigma_y = (double *)malloc(nc * sizeof(double));
        memcpy(new_problem->model.von_mises->sigma_y, model->sigma_y, nc * sizeof(double));
      } else {
        new_problem->model.von_mises = NULL;
      }
      break;
    }
    default:
      new_problem->model.generic = NULL;
      break;
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
