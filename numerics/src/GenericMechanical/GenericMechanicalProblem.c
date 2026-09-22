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

#include "GenericMechanicalProblem.h"

#include <assert.h>  // for assert
#include <stdlib.h>  // for malloc, free

#include "FrictionContactProblem.h"        // IWYU pragma: keep
#include "LinearComplementarityProblem.h"  // IWYU pragma: keep
#include "NumericsMatrix.h"                // for NumericsMatrix, NM_new
#include "RelayProblem.h"                  // IWYU pragma: keep
#include "SparseBlockMatrix.h"             // IWYU pragma: keep
#include "numerics_verbose.h"  // for check_io, numerics_error, numerics_warning
#include "safe_casts.h"
#include "siconos_debug.h"  // for DEBUG_EXPR

GenericMechanicalProblem* genericMechanicalProblem_new(void) {
  GenericMechanicalProblem* paux =
      (GenericMechanicalProblem*)malloc(sizeof(GenericMechanicalProblem));
  *paux = (GenericMechanicalProblem){0};
  return paux;
}

void genericMechanicalProblem_free(GenericMechanicalProblem* problem, unsigned int level) {
  if (!problem) return;
  while (problem->lastLocal) {
    GMP_LocalProblem* local = problem->lastLocal;

    switch (local->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP: {
        /* q_local is owned by the LCP problem and will be freed there. */
        local->q_local = NULL;
        freeLinearComplementarityProblem((LinearComplementarityProblem*)(local->problem));
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_RELAY: {
        /* q_local is owned by the Relay problem and will be freed there. */
        local->q_local = NULL;
        freeRelay_problem((RelayProblem*)(local->problem));
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_FC2D:
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        /* q_local is owned by the FC problem and will be freed there. */
        local->q_local = NULL;
        frictionContactProblem_free((FrictionContactProblem*)(local->problem));
        break;
      }
      default:
        (void)numerics_warning("genericMechanicalProblem_free",
                               "unknown problem type, possible memory leak");
    }

    if (local->q_local) free(local->q_local);
    problem->lastLocal = local->previous;
    *local = (GMP_LocalProblem){0};
    free(local);
  }

  if (level == GMP_FREE_MATRIX) {
    problem->M = NM_free(problem->M);
    if (problem->q) free(problem->q);
    problem->q = NULL;
  }
  free(problem);
}

void* gmp_add(GenericMechanicalProblem* problem, SICONOS_NUMERICS_PROBLEM_TYPE problemType,
              size_t size) {
  if (!problem) {
    (void)numerics_error("gmp_add", "NULL GenericMechanicalProblem");
    return NULL;
  }

  GMP_LocalProblem* newProblem = (GMP_LocalProblem*)malloc(sizeof(GMP_LocalProblem));
  *newProblem = (GMP_LocalProblem){0};
  newProblem->type = problemType;
  newProblem->size = size;

  problem->globalSize += size;
  if (size > problem->maxLocalSize) problem->maxLocalSize = size;

  if (!problem->lastLocal) {
    problem->firstLocal = newProblem;
    problem->lastLocal = newProblem;
  } else {
    problem->lastLocal->next = newProblem;
    newProblem->previous = problem->lastLocal;
    problem->lastLocal = newProblem;
  }

  switch (problemType) {
    case SICONOS_NUMERICS_PROBLEM_LCP: {
      newProblem->problem = (void*)malloc(sizeof(LinearComplementarityProblem));
      LinearComplementarityProblem* pLCP = (LinearComplementarityProblem*)newProblem->problem;
      pLCP->M = NM_new();
      newProblem->q_local = (double*)malloc(size * sizeof(double));
      pLCP->q = newProblem->q_local;
      pLCP->M->storageType = NM_DENSE; /* local problem is dense */
      pLCP->M->size0 = to_int(size);
      pLCP->M->size1 = to_int(size);
      pLCP->size = to_int(size);
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_RELAY: {
      newProblem->problem = (void*)malloc(sizeof(RelayProblem));
      RelayProblem* pRelay = (RelayProblem*)newProblem->problem;
      pRelay->M = NM_new();
      newProblem->q_local = (double*)malloc(size * sizeof(double));
      pRelay->q = newProblem->q_local;
      pRelay->M->storageType = NM_DENSE; /* local problem is dense */
      pRelay->M->size0 = to_int(size);
      pRelay->M->size1 = to_int(size);
      pRelay->size = to_int(size);
      pRelay->lb = (double*)malloc(size * sizeof(double));
      pRelay->ub = (double*)malloc(size * sizeof(double));
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
      newProblem->problem = NULL;
      newProblem->q_local = (double*)malloc(size * sizeof(double));
      return newProblem->q_local;
    }
    case SICONOS_NUMERICS_PROBLEM_FC3D: {
      newProblem->problem = (void*)malloc(sizeof(FrictionContactProblem));
      FrictionContactProblem* pFC3D = (FrictionContactProblem*)newProblem->problem;
      pFC3D->mu = (double*)malloc(sizeof(double));
      pFC3D->M = NM_new();
      pFC3D->M->storageType = NM_DENSE; /* local problem is dense */
      pFC3D->M->size0 = to_int(size);
      pFC3D->M->size1 = to_int(size);
      pFC3D->numberOfContacts = 1;
      newProblem->q_local = (double*)malloc(size * sizeof(double));
      pFC3D->q = newProblem->q_local;
      pFC3D->dimension = 3;
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_FC2D: {
      newProblem->problem = (void*)malloc(sizeof(FrictionContactProblem));
      FrictionContactProblem* pFC2D = (FrictionContactProblem*)newProblem->problem;
      pFC2D->mu = (double*)malloc(sizeof(double));
      pFC2D->M = NM_new();
      pFC2D->M->storageType = NM_DENSE; /* local problem is dense */
      pFC2D->M->size0 = to_int(size);
      pFC2D->M->size1 = to_int(size);
      pFC2D->numberOfContacts = 1;
      newProblem->q_local = (double*)malloc(size * sizeof(double));
      pFC2D->q = newProblem->q_local;
      pFC2D->dimension = 2;
      break;
    }
    default:
      (void)numerics_error("gmp_add", "unknown problem type");
      free(newProblem);
      return NULL;
  }
  return newProblem->problem;
}

void genericMechanicalProblem_display(GenericMechanicalProblem* problem) {
  if (!problem) {
    (void)numerics_warning("genericMechanicalProblem_display", "NULL problem");
    return;
  }
  GMP_LocalProblem* local = problem->firstLocal;
  printf("\nBEGIN Display a GenericMechanicalProblem(Numerics):\n");

  while (local) {
    printf("-->A sub-problem %s.\n", ns_problem_id_to_name(local->type));
    local = local->next;
  }
  printf("The sparse block matrix is :\n");
  NM_display(problem->M);
  printf("The q vector is :\n");
  for (size_t ii = 0; ii < problem->globalSize; ii++) printf("%e ", problem->q[ii]);

  printf("\nEND Display a GenericMechanicalProblem:\n");
}

void genericMechanicalProblem_printInFile(GenericMechanicalProblem* problem, FILE* file) {
  if (!problem || !file) return;
  GMP_LocalProblem* local_problem = problem->firstLocal;
  /* Print M */
  NM_write_in_file(problem->M, file);
  fprintf(file, "\n");
  /* Print Q */
  for (size_t ii = 0; ii < problem->globalSize; ii++) fprintf(file, "%e\n", problem->q[ii]);
  fprintf(file, "\n");
  /* Print the type and options (mu) */
  while (local_problem) {
    fprintf(file, "%d\n", local_problem->type);
    if (local_problem->type == SICONOS_NUMERICS_PROBLEM_FC3D)
      fprintf(file, "%e\n", ((FrictionContactProblem*)local_problem->problem)->mu[0]);
    local_problem = local_problem->next;
  }
}

GenericMechanicalProblem* genericMechanical_newFromFile(FILE* file) {
  if (!file) {
    (void)numerics_error("genericMechanical_newFromFile", "NULL file descriptor");
    return NULL;
  }

  GenericMechanicalProblem* problem = genericMechanicalProblem_new();
  size_t nsubProb = 0;
  int prbType = 0;
  size_t global_offset, local_size;
  void* prb;

  problem->M = NM_new_from_file(file);
  if (!problem->M) {
    genericMechanicalProblem_free(problem, GMP_FREE_GMP);
    return NULL;
  }
  SparseBlockStructuredMatrix* m = problem->M->matrix1;

  problem->q = (double*)malloc(to_size_t(problem->M->size1) * sizeof(double));
  for (size_t i = 0; i < to_size_t(problem->M->size1); i++) {
    check_io(fscanf(file, "%lf ", problem->q + i));
  }
  nsubProb = m->filled1 - 1;
  global_offset = 0;
  for (size_t blockIndex = 0; blockIndex < nsubProb; blockIndex++) {
    if (blockIndex) global_offset = m->blocksize0[blockIndex - 1];
    local_size = m->blocksize0[blockIndex] - global_offset;
    check_io(fscanf(file, "%d\n", &prbType));
    prb = gmp_add(problem, prbType, local_size);
    if (prbType == SICONOS_NUMERICS_PROBLEM_FC3D) {
      check_io(fscanf(file, "%lf ", ((FrictionContactProblem*)prb)->mu));
    }
  }

  DEBUG_EXPR(genericMechanicalProblem_display(problem););
  return problem;
}

GenericMechanicalProblem* genericMechanical_new_from_filename(const char* filename) {
  if (!filename) {
    (void)numerics_error("genericMechanical_new_from_filename", "NULL filename");
    return NULL;
  }

  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    (void)numerics_error("genericMechanical_new_from_filename", "could not open file");
    return NULL;
  }

  GenericMechanicalProblem* problem = genericMechanical_newFromFile(file);
  fclose(file);
  return problem;
}

/** Return the non-smooth problem formulation name from its id number. */
const char* ns_problem_id_to_name(SICONOS_NUMERICS_PROBLEM_TYPE id) {
  switch (id) {
    case SICONOS_NUMERICS_PROBLEM_LCP:
      return "LCP";
    case SICONOS_NUMERICS_PROBLEM_MLCP:
      return "MLCP";
    case SICONOS_NUMERICS_PROBLEM_NCP:
      return "NCP";
    case SICONOS_NUMERICS_PROBLEM_MCP:
      return "MCP";
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
      return "EQUALITY";
    case SICONOS_NUMERICS_PROBLEM_FC2D:
      return "FC2D";
    case SICONOS_NUMERICS_PROBLEM_FC3D:
      return "FC3D";
    case SICONOS_NUMERICS_PROBLEM_VI:
      return "VI";
    case SICONOS_NUMERICS_PROBLEM_AVI:
      return "AVI";
    case SICONOS_NUMERICS_PROBLEM_RELAY:
      return "RELAY";
    default:
      (void)numerics_warning("ns_problem_id_to_name", "unknown problem id");
      return NULL;
  }
}
