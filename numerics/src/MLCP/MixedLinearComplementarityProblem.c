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
#ifndef MLCP_PROBLEM_C
#define MLCP_PROBLEM_C

#include "MixedLinearComplementarityProblem.h"

#include <assert.h>  // for assert
#include <stdio.h>   // for printf, fprintf, fscanf, NULL, FILE
#include <stdlib.h>  // for malloc, atof, free, exit, EXIT_FAILURE

#include "NumericsMatrix.h"    // for NM_dense_display, NumericsMatrix, NM_c...
#include "numerics_verbose.h"  // for CHECK_IO

void mixedLinearComplementarity_free(MixedLinearComplementarityProblem* problem) {
  if (problem->isStorageType1) {
    NM_clear(problem->M);
    free(problem->M);
    free(problem->q);
    free(problem->blocksRows);
    free(problem->blocksIsComp);
  }
  if (problem->isStorageType2) {
    free(problem->A);
    free(problem->B);
    free(problem->C);
    free(problem->D);
    free(problem->a);
    free(problem->b);
  }
  free(problem);
}

MixedLinearComplementarityProblem* mixedLinearComplementarity_new(void) {
  MixedLinearComplementarityProblem* mlcp =
      (MixedLinearComplementarityProblem*)malloc(sizeof(MixedLinearComplementarityProblem));

  mlcp->isStorageType1 = 0;
  mlcp->isStorageType2 = 0;
  mlcp->m = 0;
  mlcp->n = 0;

  mlcp->blocksRows = NULL;
  mlcp->blocksIsComp = NULL;

  mlcp->M = NULL;
  mlcp->q = NULL;
  mlcp->A = NULL;
  mlcp->B = NULL;
  mlcp->C = NULL;
  mlcp->D = NULL;
  mlcp->a = NULL;
  mlcp->b = NULL;

  return mlcp;
}

void mixedLinearComplementarity_display(MixedLinearComplementarityProblem* p) {
  int n = p->n;
  int m = p->m;
  numerics_printf("MLCP DISPLAY:\n-------------\n");
  numerics_printf_verbose(0, "n :%d m: %d\n", p->n, p->m);

  numerics_printf_verbose(0, p->isStorageType1 ? "using (M)\n" : "not using (M)\n");
  numerics_printf_verbose(0, p->isStorageType2 ? "using (ABCD)\n" : "not using (ABCD)\n");
  if (p->blocksRows) {
    numerics_printf_verbose(0, "blocks are:\n");
    int NumBlock = 0;
    while (p->blocksRows[NumBlock] < n + m) {
      if (p->blocksIsComp[NumBlock]) {
        numerics_printf_verbose(
            0, "->block of complementarity condition (type %d), from line %d, to line %d.\n",
            p->blocksIsComp[NumBlock], p->blocksRows[NumBlock],
            p->blocksRows[NumBlock + 1] - 1);
      } else {
        numerics_printf_verbose(
            0, "->block of equality type (type %d), from line %d, to line %d.\n",
            p->blocksIsComp[NumBlock], p->blocksRows[NumBlock],
            p->blocksRows[NumBlock + 1] - 1);
      }
      NumBlock++;
    }
  }

  if (p->M) {
    numerics_printf_verbose(0, "M matrix:");
    NM_display(p->M);
  } else
    numerics_printf_verbose(0, "No M matrix:");

  if (p->q) {
    numerics_printf_verbose(0, "q matrix:");
    NM_dense_display(p->q, n + m, 1, 0);
  } else
    numerics_printf_verbose(0, "No q matrix:");

  if (p->A) {
    numerics_printf_verbose(0, "A matrix:");
    NM_dense_display(p->A, n, n, 0);
  } else {
    numerics_printf_verbose(0, "No A matrix:");
    if (p->M && !p->M->storageType) {
      numerics_printf_verbose(0, "A matrix from M:");
      NM_dense_display(p->M->matrix0, n, n, n + m);
    }
  }
  if (p->B) {
    numerics_printf_verbose(0, "B matrix:");
    NM_dense_display(p->B, m, m, 0);
  } else {
    numerics_printf_verbose(0, "No B matrix:");
    if (p->M && !p->M->storageType) {
      numerics_printf_verbose(0, "B matrix from M:");
      NM_dense_display(p->M->matrix0 + n * (n + m) + n, m, m, n + m);
    }
  }

  if (p->C) {
    numerics_printf_verbose(0, "C matrix:");
    NM_dense_display(p->C, n, m, 0);
  } else {
    numerics_printf_verbose(0, "No C matrix:");
    if (p->M && !p->M->storageType) {
      numerics_printf_verbose(0, "C matrix from M:");
      NM_dense_display(p->M->matrix0 + n * (n + m), n, m, n + m);
    }
  }

  if (p->D) {
    numerics_printf_verbose(0, "D matrix:");
    NM_dense_display(p->D, m, n, 0);
  } else {
    numerics_printf_verbose(0, "No D matrix:");
    if (p->M && !p->M->storageType) {
      numerics_printf_verbose(0, "D matrix from M:");
      NM_dense_display(p->M->matrix0 + n, m, n, n + m);
    }
  }
  if (p->a) {
    numerics_printf_verbose(0, "a matrix:");
    NM_dense_display(p->a, n, 1, 0);
  } else
    numerics_printf_verbose(0, "No a matrix:");
  if (p->b) {
    numerics_printf_verbose(0, "b matrix:");
    NM_dense_display(p->b, m, 1, 0);
  } else
    numerics_printf_verbose(0, "No b matrix:");
}
int mixedLinearComplementarity_printInFile(MixedLinearComplementarityProblem* problem,
                                           FILE* file) {
  int info = 0;
  if (!problem) {
    fprintf(stderr,
            "Numerics, MixedLinearComplementarityProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i, j;
  fprintf(file, "%d\n", problem->isStorageType1);
  fprintf(file, "%d\n", problem->isStorageType2);
  int n = problem->n;
  fprintf(file, "%d\n", n);
  int m = problem->m;
  fprintf(file, "%d\n", m);

  if (problem->isStorageType1) {
    int nbBlocks = 0;
    fprintf(file, "%d ", problem->blocksRows[nbBlocks]);
    while (problem->blocksRows[nbBlocks] < (n + m)) {
      nbBlocks++;
      fprintf(file, "%d ", problem->blocksRows[nbBlocks]);
    }
    fprintf(file, "\n");
    for (i = 0; i < nbBlocks; i++) {
      fprintf(file, "%d ", problem->blocksIsComp[i]);
    }
    fprintf(file, "\n");
    NM_write_in_file(problem->M, file);

    for (i = 0; i < problem->M->size1; i++) {
      fprintf(file, "%32.24e ", problem->q[i]);
    }
    fprintf(file, "\n");
    /* return 1; */
    /* if (problem->isStorageType2)  */
    /* { */
    /*   printf("Numerics, MixedLinearComplementarityProblem printInFile only Storage1 has been
     * printed.\n"); */
    /* } */
  }
  if (problem->isStorageType2) {
    for (i = 0; i < problem->n; i++) {
      for (j = 0; j < problem->n; j++) {
        fprintf(file, "%32.24e ", problem->A[i + j * n]);
      }
      fprintf(file, "\n");
    }

    for (i = 0; i < problem->m; i++) {
      for (j = 0; j < problem->m; j++) {
        fprintf(file, "%32.24e ", problem->B[i + j * m]);
      }
      fprintf(file, "\n");
    }
    for (i = 0; i < problem->n; i++) {
      for (j = 0; j < problem->m; j++) {
        fprintf(file, "%32.24e ", problem->C[i + j * n]);
      }
      fprintf(file, "\n");
    }
    for (i = 0; i < problem->m; i++) {
      for (j = 0; j < problem->n; j++) {
        fprintf(file, "%32.24e ", problem->D[i + j * m]);
      }
      fprintf(file, "\n");
    }

    for (i = 0; i < problem->n; i++) {
      fprintf(file, "%32.24e ", problem->a[i]);
    }
    fprintf(file, "\n");
    for (i = 0; i < problem->m; i++) {
      fprintf(file, "%32.24e ", problem->b[i]);
    }
  }
  return info;
}

int mixedLinearComplementarity_newFromFile(MixedLinearComplementarityProblem* problem,
                                           FILE* file) {
  int info = 0;
  assert(file);
  if (!problem) {
    fprintf(stderr,
            "Numerics, MixedLinearComplementarityProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i, j;
  int st1, st2;
  CHECK_IO(fscanf(file, "%d\n", &st1));
  problem->isStorageType1 = st1;
  CHECK_IO(fscanf(file, "%d\n", &st2));
  problem->isStorageType2 = st2;

  int n;
  CHECK_IO(fscanf(file, "%d\n", &n));
  problem->n = n;
  int m;
  CHECK_IO(fscanf(file, "%d\n", &m));
  problem->m = m;

  if (problem->isStorageType1) {
    int* blocksRows = (int*)malloc((n + m + 1) * sizeof(int));
    int nbBlocks = 0;
    CHECK_IO(fscanf(file, "%d ", &(blocksRows[nbBlocks])));
    while (blocksRows[nbBlocks] < (n + m)) {
      nbBlocks++;
      CHECK_IO(fscanf(file, "%d ", &(blocksRows[nbBlocks])));
    }
    problem->blocksRows = (int*)malloc((nbBlocks + 1) * sizeof(int));
    // CHECK_IO(fscanf(file,"\n"));
    for (i = 0; i <= nbBlocks; i++) {
      problem->blocksRows[i] = blocksRows[i];
    }
    free(blocksRows);
    problem->blocksIsComp = (int*)malloc((nbBlocks) * sizeof(int));

    // fprintf(file,"\n");
    for (i = 0; i < nbBlocks; i++) {
      CHECK_IO(fscanf(file, "%d ", &(problem->blocksIsComp[i])));
    }

    // fprintf(file,"\n");
    problem->M = NM_new_from_file(file);
    problem->q = (double*)malloc(problem->M->size1 * sizeof(double));

    for (i = 0; i < problem->M->size1; i++) {
      CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
    }
    // fprintf(file,"\n");
    /* return 1; */
    /* if (problem->isStorageType2)  */
    /* { */
    /*   printf("Numerics, MixedLinearComplementarityProblem printInFile only Storage1 has been
     * printed.\n"); */
    /* } */
  }
  if (problem->isStorageType2) {
    problem->A = (double*)malloc(n * n * sizeof(double));
    problem->B = (double*)malloc(m * m * sizeof(double));
    problem->C = (double*)malloc(n * m * sizeof(double));
    problem->D = (double*)malloc(m * n * sizeof(double));
    problem->a = (double*)malloc(n * sizeof(double));
    problem->b = (double*)malloc(m * sizeof(double));

    for (i = 0; i < problem->n; i++) {
      for (j = 0; j < problem->n; j++) {
        CHECK_IO(fscanf(file, "%lf ", &(problem->A[i + j * n])));
      }
      /* CHECK_IO(fscanf(file,"\n")); */
    }

    for (i = 0; i < problem->m; i++) {
      for (j = 0; j < problem->m; j++) {
        CHECK_IO(fscanf(file, "%lf ", &(problem->B[i + j * m])));
      }
      /* fprintf(file,"\n"); */
    }
    for (i = 0; i < problem->n; i++) {
      for (j = 0; j < problem->m; j++) {
        CHECK_IO(fscanf(file, "%lf ", &(problem->C[i + j * n])));
      }
      /* fprintf(file,"\n"); */
    }
    for (i = 0; i < problem->m; i++) {
      for (j = 0; j < problem->n; j++) {
        CHECK_IO(fscanf(file, "%lf ", &(problem->D[i + j * m])));
      }
      /* fprintf(file,"\n"); */
    }

    for (i = 0; i < problem->n; i++) {
      CHECK_IO(fscanf(file, "%lf ", &(problem->a[i])));
    }
    /* fprintf(file,"\n"); */
    for (i = 0; i < problem->m; i++) {
      CHECK_IO(fscanf(file, "%lf ", &(problem->b[i])));
    }
  }
  return info;
}

int mixedLinearComplementarity_newFromFileOld(MixedLinearComplementarityProblem* problem,
                                              FILE* file) {
  int n = 0, m = 0, NbLines = 0;
  int i, j, m2;
  char val[128];

  double *vecA, *vecB, *vecC, *vecD, *vecM, *vecQ;
  double *a, *b;
  CHECK_IO(fscanf(file, "%d", &n));
  CHECK_IO(fscanf(file, "%d", &m));
  CHECK_IO(fscanf(file, "%d", &NbLines));

  m2 = m * m;

  vecM = (double*)malloc((n + m) * (NbLines) * sizeof(double));

  vecQ = (double*)malloc((NbLines) * sizeof(double));
  vecA = (double*)malloc(n * (NbLines - m) * sizeof(double));
  vecB = (double*)malloc(m2 * sizeof(double));
  vecC = (double*)malloc((NbLines - m) * m * sizeof(double));
  vecD = (double*)malloc(m * n * sizeof(double));
  a = (double*)malloc((NbLines - m) * sizeof(double));
  b = (double*)malloc(m * sizeof(double));

  problem->blocksRows = (int*)malloc(3 * sizeof(int));
  problem->blocksIsComp = (int*)malloc(2 * sizeof(int));
  problem->blocksRows[0] = 0;
  problem->blocksRows[1] = n;
  problem->blocksRows[2] = n + m;
  problem->blocksIsComp[0] = 0;
  problem->blocksIsComp[1] = 1;

  problem->M = NM_create_from_data(NM_DENSE, NbLines, n + m, vecM);

  problem->isStorageType1 = 1;  // Both problems seems to be stored
  problem->isStorageType2 = 1;  // Both problems seems to be stored

  problem->q = vecQ;
  problem->A = vecA;
  problem->B = vecB;
  problem->C = vecC;
  problem->D = vecD;
  problem->a = a;
  problem->b = b;
  problem->blocksRows[1] = n;
  problem->blocksRows[2] = n + m;
  problem->n = n;
  problem->m = m;

  for (i = 0; i < NbLines - m; ++i) {
    for (j = 0; j < n; ++j) {
      CHECK_IO(fscanf(file, "%s", val));
      vecA[(NbLines - m) * j + i] = atof(val);
      vecM[(NbLines)*j + i] = atof(val);
    }
  }
  for (i = 0; i < m; ++i) {
    for (j = 0; j < m; ++j) {
      CHECK_IO(fscanf(file, "%s", val));
      vecB[m * j + i] = atof(val);
      /*  vecM[ n*(m+n)+(n+m)*j+n+i ] = atof(val);*/
      vecM[n * (NbLines) + (NbLines)*j + (NbLines - m) + i] = atof(val);
    }
  }
  for (i = 0; i < NbLines - m; ++i) {
    for (j = 0; j < m; ++j) {
      CHECK_IO(fscanf(file, "%s", val));
      vecC[(NbLines - m) * j + i] = atof(val);
      vecM[(NbLines) * (n + j) + i] = atof(val);
    }
  }
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      CHECK_IO(fscanf(file, "%s", val));
      vecD[m * j + i] = atof(val);
      vecM[(NbLines)*j + i + (NbLines - m)] = atof(val);
    }
  }

  for (i = 0; i < NbLines - m; ++i) {
    CHECK_IO(fscanf(file, "%s", val));
    a[i] = atof(val);
    vecQ[i] = atof(val);
  }
  for (i = 0; i < m; ++i) {
    CHECK_IO(fscanf(file, "%s", val));
    b[i] = atof(val);
    vecQ[i + NbLines - m] = atof(val);
  }

  return 0;
}

int mixedLinearComplementarity_newFromFilename(MixedLinearComplementarityProblem* problem,
                                               const char* filename) {
  int info = 0;
  FILE* file = fopen(filename, "r");

  info = mixedLinearComplementarity_newFromFile(problem, file);

  fclose(file);
  return info;
}

MixedLinearComplementarityProblem* mixedLinearComplementarity_fromMtoABCD(
    MixedLinearComplementarityProblem* problem) {
  MixedLinearComplementarityProblem* mlcp = mixedLinearComplementarity_new();

  if (!problem->isStorageType1) {
    numerics_printf_verbose(
        0, "mixedLinearComplementarity_new: Wrong Storage (!isStorageType1)\n");
    exit(EXIT_FAILURE);
  }

  mlcp->isStorageType1 = 0;
  mlcp->isStorageType2 = 1;

  mlcp->n = problem->n; /* Equalities */
  mlcp->m = problem->m; /* Inequalities */

  int m = problem->m;
  int n = problem->n;
  mlcp->A = (double*)calloc(n * n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      mlcp->A[i + j * n] = problem->M->matrix0[i + j * (m + n)];
    }
  }

  mlcp->B = (double*)calloc(m * m, sizeof(double));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      mlcp->B[i + j * m] = problem->M->matrix0[i + n + (j + n) * (m + n)];
    }
  }

  mlcp->C = (double*)calloc(n * m, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      mlcp->C[i + j * m] = problem->M->matrix0[i + (j + n) * (m + n)];
    }
  }

  mlcp->D = (double*)calloc(m * n, sizeof(double));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      mlcp->D[i + j * n] = problem->M->matrix0[i + n + j * (m + n)];
    }
  }

  mlcp->a = (double*)calloc(n, sizeof(double));
  for (int i = 0; i < n; i++) {
    mlcp->a[i] = problem->q[i];
  }

  mlcp->b = (double*)calloc(m, sizeof(double));
  for (int i = 0; i < m; i++) {
    mlcp->b[i] = problem->q[i + n];
  }

  return mlcp;
}

#endif
