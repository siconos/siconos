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
#include "CohesiveFrictionContactProblem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NumericsMatrix.h"
#include "numerics_verbose.h"

CohesiveFrictionContactProblem* cohesiveFrictionContactProblem_new(void) {
  CohesiveFrictionContactProblem* problem =
      (CohesiveFrictionContactProblem*)malloc(sizeof(CohesiveFrictionContactProblem));
  problem->dimension = 0;
  problem->numberOfContacts = 0;
  problem->numberOfCohesivePoints = 0;
  problem->M = NULL;
  problem->q = NULL;
  problem->mu = NULL;
  problem->c_n = NULL;
  problem->c_t = NULL;
  problem->q_v = NULL;
  problem->q_u = NULL;
  problem->W = NULL;
  problem->V = NULL;
  problem->X = NULL;
  problem->U = NULL;
  problem->internal_state = NULL;
  problem->internal_state_size = 0;
  return problem;
}

CohesiveFrictionContactProblem* cohesiveFrictionContactProblem_new_with_data(
    int dim, int nc, int n_coh, NumericsMatrix* M, double* q, double* mu, NumericsMatrix* V,
    NumericsMatrix* X, NumericsMatrix* U) {
  CohesiveFrictionContactProblem* problem = cohesiveFrictionContactProblem_new();
  problem->dimension = dim;
  problem->numberOfContacts = nc;
  problem->numberOfCohesivePoints = n_coh;
  problem->M = M;
  problem->q = q;
  problem->mu = mu;
  problem->V = V;
  problem->X = X;
  problem->U = U;
  problem->internal_state = NULL;
  problem->internal_state_size = 0;
  return problem;
}

void cohesiveFrictionContactProblem_free(CohesiveFrictionContactProblem* problem) {
  if (problem) {
    if (problem->M) {
      NM_free(problem->M);
      problem->M = NULL;
    }
    if (problem->W) {
      NM_free(problem->W);
      problem->W = NULL;
    }
    if (problem->V) {
      NM_free(problem->V);
      problem->V = NULL;
    }
    if (problem->X) {
      NM_free(problem->X);
      problem->X = NULL;
    }
    if (problem->U) {
      NM_free(problem->U);
      problem->U = NULL;
    }
    free(problem);
  }
}

void cohesiveFrictionContact_display(CohesiveFrictionContactProblem* problem) {
  if (!problem) {
    fprintf(stderr, "CohesiveFrictionContact_display :: NULL pointer");
    return;
  }
  printf("CohesiveFrictionContactProblem Display:\n");
  printf("=====================================\n");
  printf("Dimension: %d\n", problem->dimension);
  printf("Number of contacts: %d\n", problem->numberOfContacts);
  printf("Number of cohesive points: %d\n", problem->numberOfCohesivePoints);

  if (problem->M) {
    printf("Matrix M:\n");
    NM_display(problem->M);
  }

  if (problem->q) {
    printf("Vector q (size %d):\n",
           problem->dimension * (problem->numberOfContacts + problem->numberOfCohesivePoints));
    for (int i = 0; i < problem->dimension *
                            (problem->numberOfContacts + problem->numberOfCohesivePoints);
         i++) {
      printf("  q[%d] = %e\n", i, problem->q[i]);
    }
  }

  if (problem->mu) {
    printf("Friction coefficients mu:\n");
    for (int i = 0; i < problem->numberOfContacts; i++) {
      printf("  mu[%d] = %e\n", i, problem->mu[i]);
    }
  }

  if (problem->c_n) {
    printf("Normal cohesion intensity c_n:\n");
    for (int i = 0; i < problem->numberOfCohesivePoints; i++) {
      printf("  c_n[%d] = %e\n", i, problem->c_n[i]);
    }
  }

  if (problem->c_t) {
    printf("Tangent cohesion intensity c_t:\n");
    for (int i = 0; i < problem->numberOfCohesivePoints; i++) {
      printf("  c_t[%d] = %e\n", i, problem->c_t[i]);
    }
  }

  if (problem->q_v) {
    printf("Contact vector q_v:\n");
    for (int i = 0; i < problem->dimension * problem->numberOfContacts; i++) {
      printf("  q_v[%d] = %e\n", i, problem->q_v[i]);
    }
  }

  if (problem->q_u) {
    printf("Cohesive vector q_u:\n");
    for (int i = 0; i < problem->dimension * problem->numberOfCohesivePoints; i++) {
      printf("  q_u[%d] = %e\n", i, problem->q_u[i]);
    }
  }

  if (problem->W) {
    printf("Matrix W (contact-contact block):\n");
    NM_display(problem->W);
  }

  if (problem->V) {
    printf("Matrix V (cohesion mapping):\n");
    NM_display(problem->V);
  }

  if (problem->X) {
    printf("Matrix X (additional cohesive contributions):\n");
    NM_display(problem->X);
  }

  if (problem->U) {
    printf("Matrix U (coupling terms):\n");
    NM_display(problem->U);
  }

  printf("=====================================\n");
}

int cohesiveFrictionContact_printInFile(CohesiveFrictionContactProblem* problem, FILE* file) {
  if (!problem || !file) return -1;

  fprintf(file, "%d\n", problem->dimension);
  fprintf(file, "%d\n", problem->numberOfContacts);
  fprintf(file, "%d\n", problem->numberOfCohesivePoints);

  NM_write_in_file(problem->M, file);

  int m = problem->dimension * (problem->numberOfContacts + problem->numberOfCohesivePoints);
  for (int i = 0; i < m; i++) {
    fprintf(file, "%.32e\n", problem->q[i]);
  }

  for (int i = 0; i < problem->numberOfContacts; i++) {
    fprintf(file, "%.32e\n", problem->mu[i]);
  }

  // Write cohesion intensities
  if (problem->c_n) {
    for (int i = 0; i < problem->numberOfCohesivePoints; i++) {
      fprintf(file, "%.32e\n", problem->c_n[i]);
    }
  }

  if (problem->c_t) {
    for (int i = 0; i < problem->numberOfCohesivePoints; i++) {
      fprintf(file, "%.32e\n", problem->c_t[i]);
    }
  }

  // Write q_u
  if (problem->q_u) {
    for (int i = 0; i < problem->dimension * problem->numberOfCohesivePoints; i++) {
      fprintf(file, "%.32e\n", problem->q_u[i]);
    }
  }

  // Write q_v
  if (problem->q_v) {
    for (int i = 0; i < problem->dimension * problem->numberOfContacts; i++) {
      fprintf(file, "%.32e\n", problem->q_v[i]);
    }
  }

  // Write required matrices W, V, X, U
  NM_write_in_file(problem->W, file);
  NM_write_in_file(problem->V, file);
  NM_write_in_file(problem->X, file);
  NM_write_in_file(problem->U, file);

  return 0;
}

int cohesiveFrictionContact_printInFilename(CohesiveFrictionContactProblem* problem,
                                            char* filename) {
  FILE* file = fopen(filename, "w");
  if (!file) return -1;
  int info = cohesiveFrictionContact_printInFile(problem, file);
  fclose(file);
  return info;
}

CohesiveFrictionContactProblem* cohesiveFrictionContact_newFromFile(FILE* file) {
  if (!file) return NULL;

  CohesiveFrictionContactProblem* problem = cohesiveFrictionContactProblem_new();

  int dim, nc, n_coh;
  if (fscanf(file, "%d", &dim) != 1) goto fail;
  if (fscanf(file, "%d", &nc) != 1) goto fail;
  if (fscanf(file, "%d", &n_coh) != 1) goto fail;

  problem->dimension = dim;
  problem->numberOfContacts = nc;
  problem->numberOfCohesivePoints = n_coh;

  problem->M = NM_new_from_file(file);
  if (!problem->M) goto fail;

  int m = dim * (nc + n_coh);
  problem->q = (double*)malloc(m * sizeof(double));
  for (int i = 0; i < m; i++) {
    if (fscanf(file, "%lf", &problem->q[i]) != 1) goto fail;
  }

  problem->mu = (double*)malloc(nc * sizeof(double));
  for (int i = 0; i < nc; i++) {
    if (fscanf(file, "%lf", &problem->mu[i]) != 1) goto fail;
  }

  // Read cohesion intensities
  problem->c_n = (double*)malloc(n_coh * sizeof(double));
  for (int i = 0; i < n_coh; i++) {
    if (fscanf(file, "%lf", &problem->c_n[i]) != 1) goto fail;
  }

  problem->c_t = (double*)malloc(n_coh * sizeof(double));
  for (int i = 0; i < n_coh; i++) {
    if (fscanf(file, "%lf", &problem->c_t[i]) != 1) goto fail;
  }

  // Read q_u
  problem->q_u = (double*)malloc(dim * n_coh * sizeof(double));
  for (int i = 0; i < dim * n_coh; i++) {
    if (fscanf(file, "%lf", &problem->q_u[i]) != 1) goto fail;
  }

  // Read q_v
  problem->q_v = (double*)malloc(dim * nc * sizeof(double));
  for (int i = 0; i < dim * nc; i++) {
    if (fscanf(file, "%lf", &problem->q_v[i]) != 1) goto fail;
  }

  // Read required matrices W, V, X, U
  problem->W = NM_new_from_file(file);
  if (!problem->W) goto fail;

  problem->V = NM_new_from_file(file);
  if (!problem->V) goto fail;

  problem->X = NM_new_from_file(file);
  if (!problem->X) goto fail;

  problem->U = NM_new_from_file(file);
  if (!problem->U) goto fail;

  return problem;

fail:
  cohesiveFrictionContactProblem_free(problem);
  return NULL;
}

CohesiveFrictionContactProblem* cohesiveFrictionContact_newFromFilename(char* filename) {
  FILE* file = fopen(filename, "r");
  if (!file) return NULL;
  CohesiveFrictionContactProblem* problem = cohesiveFrictionContact_newFromFile(file);
  fclose(file);
  return problem;
}

void cohesiveFrictionContactProblem_set_cohesion(CohesiveFrictionContactProblem* problem,
                                                 double* c_n, double* c_t) {
  if (!problem) return;
  problem->c_n = c_n;
  problem->c_t = c_t;
}

void cohesiveFrictionContactProblem_compute_effective_q(
    CohesiveFrictionContactProblem* problem, double* q_eff) {
  if (!problem || !q_eff) return;

  int m = problem->dimension * problem->numberOfContacts;

  // Start with q
  memcpy(q_eff, problem->q, m * sizeof(double));

  // Add cohesive contribution using V matrix and cohesion intensities
  // Note: The actual computation depends on how c_n and c_t are mapped
  // through V, X, U matrices. This is a placeholder implementation.
  if (problem->c_n && problem->c_t && problem->V) {
    // Combine c_n and c_t into a single cohesion vector if needed
    // Then: q_eff = q + V * c
    // For now, just a placeholder - actual implementation depends on
    // the specific formulation
  }
}

int cohesiveFrictionContactProblem_build_M_q_from_blocks(
    CohesiveFrictionContactProblem* problem) {
  if (!problem) return -1;

  int dim = problem->dimension;
  int nc = problem->numberOfContacts;
  int n_coh = problem->numberOfCohesivePoints;
  int storageType = -1;

  int contact_friction_case = -1;

  /* // Check that all required blocks are present */
  if ((nc > 0) && (n_coh > 0)) {
    storageType = problem->W->storageType;
    contact_friction_case = 0;
    if (!problem->W || !problem->V || !problem->U || !problem->X) {
      fprintf(stderr,
              "cohesiveFrictionContactProblem_build_M_q_from_blocks: "
              "Missing block matrices (W, V, U, X must all be set)\n");
      return -1;
    }
  } else if ((nc == 0) && (n_coh > 0)) {
    contact_friction_case = 1;
    if (!problem->X) {
      storageType = problem->X->storageType;
      fprintf(stderr,
              "cohesiveFrictionContactProblem_build_M_q_from_blocks: "
              "Missing block matrix (X must  be set)\n");
      return -1;
    }
  } else if ((nc > 0) && (n_coh == 0)) {
    storageType = problem->W->storageType;
    contact_friction_case = 2;
    if (!problem->W) {
      fprintf(stderr,
              "cohesiveFrictionContactProblem_build_M_q_from_blocks: "
              "Missing block matrix (W must  be set)\n");
      return -1;
    }
  }

  if ((nc > 0) && (n_coh > 0)) {
    if (!problem->q_v || !problem->q_u) {
      fprintf(stderr,
              "cohesiveFrictionContactProblem_build_M_q_from_blocks: "
              "Missing block vectors (q_v, q_u must all be set)\n");
      return -1;
    }
  }

  // Compute sizes
  int m = dim * nc;     // size of velocity/contact block
  int n = dim * n_coh;  // size of cohesive/displacement block
  int M_size = m + n;   // total size of M

  // Create or reuse M matrix
  if (problem->M) {
    // Check if existing M has correct size
    if (problem->M->size0 != M_size || problem->M->size1 != M_size) {
      NM_free(problem->M);
      problem->M = NULL;
    }
  }

  if (!problem->M) {
    // Create new M matrix with appropriate size
    // Use sparse representation if the storage is sparse, otherwise dense
    if (storageType == NM_SPARSE) {
      problem->M = NM_create(NM_SPARSE, M_size, M_size);
      NM_triplet_alloc(problem->M, M_size);
    } else {
      problem->M = NM_create(NM_DENSE, M_size, M_size);
    }
  }

  // Build block matrix M = [[W, V], [U, X]]
  // Using NM_insert to place each block
  if (contact_friction_case == 0) {
    // Insert W at position (0, 0)
    NM_insert(problem->M, problem->W, 0, 0);

    // Insert V at position (0, m)
    NM_insert(problem->M, problem->V, 0, m);

    // Insert U at position (m, 0)
    NM_insert(problem->M, problem->U, m, 0);

    // Insert X at position (m, m)
    NM_insert(problem->M, problem->X, m, m);
  } else if (contact_friction_case == 1) {
    // Insert X at position (m, m)
    NM_insert(problem->M, problem->X, m, m);
  } else if (contact_friction_case == 2) {
    // Insert W at position (0, 0)
    NM_insert(problem->M, problem->W, 0, 0);
  }

  // Create or reuse q vector
  if (problem->q) {
    free(problem->q);
  }
  problem->q = (double*)malloc(M_size * sizeof(double));

  // Build q = [q_v, q_u] by concatenation
  memcpy(problem->q, problem->q_v, m * sizeof(double));
  memcpy(problem->q + m, problem->q_u, n * sizeof(double));

  if (verbose > 0) {
    printf(
        "cohesiveFrictionContactProblem_build_M_q_from_blocks: "
        "Built M (%dx%d) and q (%d) from blocks\n",
        M_size, M_size, M_size);
    if (contact_friction_case == 0) {
      printf("  W: %dx%d, V: %dx%d, U: %dx%d, X: %dx%d\n", problem->W->size0,
             problem->W->size1, problem->V->size0, problem->V->size1, problem->U->size0,
             problem->U->size1, problem->X->size0, problem->X->size1);
      printf("  q_v: %d, q_u: %d\n", m, n);
    } else if (contact_friction_case == 1) {
      printf("  X: %dx%d\n", problem->X->size0, problem->X->size1);
    } else if (contact_friction_case == 2) {
      printf("  W: %dx%d\n", problem->W->size0, problem->W->size1);
    }
  }

  return 0;
}
