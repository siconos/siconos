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

#include "SiconosConfig.h"  // for HAVE_PATHFERRIS // IWYU pragma: keep
#include "SimpleLCP.h"      // for ABCDtoM, FortranToPathSparse, MCP_Termination

#ifdef HAVE_PATHFERRIS
#include <limits.h>
#include <stdio.h>

#define PATHFERRIS_LOG_IN_FILE

#if defined(__cplusplus)
extern "C" {
#endif

#include "PATH_SDK/include/Error.h"
#include "PATH_SDK/include/MCP_Interface.h"
#include "PATH_SDK/include/Macros.h"
#include "PATH_SDK/include/Memory.h"
#include "PATH_SDK/include/Options.h"
#include "PATH_SDK/include/Output.h"
#include "PATH_SDK/include/Output_Interface.h"
#include "PATH_SDK/include/Path.h"
#include "PATH_SDK/include/PathOptions.h"

#if defined(__cplusplus)
}
#endif

typedef struct {
  int variables;

  int n;
  int nnz;

  double *z;
  double *lb;
  double *ub;

  int *m_start;
  int *m_len;
  int *m_row;
  double *m_data;

  double *q;
} Problem;

static Problem problem;
static int filled;

static CB_FUNC(void) start(void *v) {
  filled = 0;
  return;
}

static CB_FUNC(void) problem_size(void *v, int *n, int *nnz) {
  *n = problem.n;
  *nnz = problem.nnz + 1;
  return;
}

static CB_FUNC(void) bounds(void *v, int n, double *z, double *lb, double *ub) {
  int i;

  for (i = 0; i < n; i++) {
    z[i] = problem.z[i];
    lb[i] = problem.lb[i];
    ub[i] = problem.ub[i];
  }
  return;
}

static CB_FUNC(int) function_evaluation(void *v, int n, double *z, double *f) {
  int col, colStart, colEnd, row;
  double value;

  for (col = 0; col < n; col++) {
    f[col] = problem.q[col];
  }

  for (col = 0; col < n; col++) {
    value = z[col];

    if (value != 0) {
      colStart = problem.m_start[col] - 1;
      colEnd = colStart + problem.m_len[col];

      while (colStart < colEnd) {
        row = problem.m_row[colStart] - 1;
        f[row] += problem.m_data[colStart] * value;
        colStart++;
      }
    }
  }

  return 0;
}

static CB_FUNC(int)
    jacobian_evaluation(void *v, int n, double *z, int wantf, double *f, int *nnz,
                        int *col_start, int *col_len, int *row, double *data) {
  int element;

  if (wantf) {
    function_evaluation(v, n, z, f);
  }

  if (!filled) {
    for (element = 0; element < problem.n; element++) {
      col_start[element] = problem.m_start[element];
      col_len[element] = problem.m_len[element];
    }

    for (element = 0; element < problem.nnz; element++) {
      row[element] = problem.m_row[element];
      data[element] = problem.m_data[element];
    }

    filled = 1;
  }

  *nnz = problem.nnz;
  return 0;
}

static CB_FUNC(void) mcp_typ(void *d, int nnz, int *typ) {
  int i;

  for (i = 0; i < nnz; i++) {
    typ[i] = PRESOLVE_LINEAR;
  }
  return;
}

static MCP_Interface mcp_interface = {
    NULL, problem_size, bounds, function_evaluation, jacobian_evaluation, NULL, start, NULL,
    NULL, NULL,         NULL};

static Presolve_Interface mcp_presolve = {NULL, NULL, NULL, NULL, NULL, mcp_typ, NULL};

static void install_interface(MCP *m) {
  MCP_SetInterface(m, &mcp_interface);
  MCP_SetPresolveInterface(m, &mcp_presolve);
  return;
}

static void sort(int rows, int cols, int elements, int *row, int *col, double *data) {
  double *m_data;
  int *m_start;
  int *m_len;
  int *m_row;

  int i, cs, ce;

  m_start = (int *)Memory_Allocate(sizeof(int) * (cols + 1));
  m_len = (int *)Memory_Allocate(sizeof(int) * (cols + 1));
  m_row = (int *)Memory_Allocate(sizeof(int) * (elements + 1));
  m_data = (double *)Memory_Allocate(sizeof(double) * (elements + 1));

  for (i = 0; i < cols; i++) {
    m_len[i] = 0;
  }

  for (i = 0; i < elements; i++) {
    if ((col[i] < 1) || (col[i] > cols)) {
      Error("column incorrect.\n");
    }

    if ((row[i] < 1) || (row[i] > rows)) {
      Error("column incorrect.\n");
    }

    m_len[col[i] - 1]++;
  }

  m_start[0] = 0;
  for (i = 1; i < cols; i++) {
    m_start[i] = m_start[i - 1] + m_len[i - 1];
    m_len[i - 1] = 0;
  }
  m_len[i - 1] = 0;

  for (i = 0; i < elements; i++) {
    cs = col[i] - 1;
    ce = m_start[cs] + m_len[cs];
    m_row[ce] = row[i];
    m_data[ce] = data[i];
    m_len[cs]++;
  }

  elements = 0;
  for (i = 0; i < cols; i++) {
    cs = m_start[i];
    ce = cs + m_len[i];

    while (cs < ce) {
      row[elements] = m_row[cs];
      col[elements] = i + 1;
      data[elements] = m_data[cs];
      elements++;
      cs++;
    }
  }

  Memory_Free(m_data);
  Memory_Free(m_row);
  Memory_Free(m_len);
  Memory_Free(m_start);
  return;
}

static void create(int variables, int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
                   double *z, double *lb, double *ub) {
  int m_index;
  int m_count;
  int i;

  problem.n = variables;
  problem.nnz = m_nnz;

  problem.z = (double *)Memory_Allocate(sizeof(double) * problem.n);
  problem.lb = (double *)Memory_Allocate(sizeof(double) * problem.n);
  problem.ub = (double *)Memory_Allocate(sizeof(double) * problem.n);

  problem.m_start = (int *)Memory_Allocate(sizeof(int) * problem.n);
  problem.m_len = (int *)Memory_Allocate(sizeof(int) * problem.n);
  problem.m_row = (int *)Memory_Allocate(sizeof(int) * problem.nnz + 1);
  problem.m_data = (double *)Memory_Allocate(sizeof(double) * problem.nnz + 1);

  problem.q = (double *)Memory_Allocate(sizeof(double) * problem.n);

  sort(variables, variables, m_nnz, m_i, m_j, m_ij);

  for (i = 0; i < variables; i++) {
    problem.z[i] = z[i];

    problem.q[i] = q[i];
    problem.lb[i] = lb[i];
    problem.ub[i] = ub[i];
  }

  m_index = 0;
  m_count = 0;
  for (i = 0; i < variables; i++) {
    problem.m_start[i] = m_count + 1;
    problem.m_len[i] = 0;

    while ((m_index < m_nnz) && (m_j[m_index] <= i + 1)) {
      if (m_ij[m_index] != 0) {
        problem.m_len[i]++;
        problem.m_row[m_count] = m_i[m_index];
        problem.m_data[m_count] = m_ij[m_index];
        m_count++;
      }
      m_index++;
    }
  }
  problem.nnz = m_count;
  return;
}

static void destroy(void) {
  Memory_Free(problem.z);
  Memory_Free(problem.lb);
  Memory_Free(problem.ub);
  Memory_Free(problem.m_start);
  Memory_Free(problem.m_len);
  Memory_Free(problem.m_row);
  Memory_Free(problem.m_data);
  Memory_Free(problem.q);
  return;
}
void printLCP(int variables, int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
              double *lb, double *ub) {
  int i = 0;
  printf("********PRINT PATH INPUT***************\n");
  printf("***************************************\n");
  printf("dim: %d val non nul: %d\n", variables, m_nnz);
  for (i = 0; i < m_nnz; i++) {
    printf("%d %d %10.7f\n", m_i[i], m_j[i], m_ij[i]);
  }
  printf("bounds:\n");
  for (i = 0; i < variables; i++) {
    printf("%d | %10.7f | %10.7f\n", i, lb[i], ub[i]);
  }
  printf("***************************************\n");
}

int nbNonNulElems(int n, double *M, double tol) {
  int i, j;
  int cmp = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if (fabs(M[i + j * n]) > tol) {
        cmp++;
      }
    }
  return cmp;
}
void FortranToPathSparse(int n, double *M, double tol, int *m_i, int *m_j, double *m_ij) {
  int i, j;
  int cmp = 0;
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++) {
      if (fabs(M[i + j * n]) > tol) {
        m_i[cmp] = i + 1;
        m_j[cmp] = j + 1;
        m_ij[cmp] = M[i + j * n];
        cmp++;
      }
    }
}
void ABCDtoM(int n, int m, double *A, double *B, double *C, double *D, double *a, double *b,
             double *M, double *q) {
  int i = 0, j = 0;
  int NM = n + m;
  /*A in M*/
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) M[i + j * NM] = A[i + j * n];
  /*C in M*/
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) M[i + (j + n) * NM] = C[i + j * n];
  /*D in M*/
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) M[i + n + (j)*NM] = D[i + j * m];
  /*B in M*/
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) M[i + n + (j + n) * NM] = B[i + j * m];
  /*a in q*/
  for (i = 0; i < n; i++) q[i] = a[i];
  /*b in q*/
  for (i = 0; i < m; i++) q[n + i] = b[i];
}

void SimpleLCP(int variables, int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
               double *lb, double *ub, MCP_Termination *status, double *z) {
  Options_Interface *o;
  MCP *m;
  Information info;

  double *x;
  double dnnz;
  int i;
#ifdef PATHFERRIS_LOG_IN_FILE
  FILE *f;
  f = fopen("path.log", "w");
  Output_SetLog(f);
#endif
  o = Options_Create();
  Path_AddOptions(o);
  Options_Default(o);

  Output_Printf(Output_Log | Output_Status | Output_Listing, "%s: LCP Link\n", Path_Version());

  create(variables, m_nnz, m_i, m_j, m_ij, q, z, lb, ub);

  if (problem.n == 0) {
    Output_Printf(Output_Log | Output_Status,
                  "\n ** EXIT - solution found (degenerate model).\n");
    (*status) = MCP_Solved;
    Options_Destroy(o);
    return;
  }

  dnnz = MIN(1.0 * problem.nnz, 1.0 * problem.n * problem.n);
  if (dnnz > INT_MAX) {
    Output_Printf(Output_Log | Output_Status, "\n ** EXIT - model too large.\n");
    (*status) = MCP_Error;
    Options_Destroy(o);
    return;
  }
  problem.nnz = (int)dnnz;

  Output_Printf(Output_Log | Output_Status | Output_Listing,
                "%d row/cols, %d non-zeros, %3.2f%% dense.\n\n", problem.n, problem.nnz,
                100.0 * problem.nnz / (1.0 * problem.n * problem.n));

  m = MCP_Create(problem.n, problem.nnz + 1);
  MCP_Jacobian_Structure_Constant(m, 1);
  install_interface(m);

  Options_Read(o, "path.opt");
  Options_SetDouble(o, "con_tol", 1e-12);
  Options_Display(o);

  info.generate_output = Output_Log | Output_Status | Output_Listing;
  info.use_start = True;
  info.use_basics = True;

  (*status) = Path_Solve(m, &info);

  x = MCP_GetX(m);

  for (i = 0; i < problem.n; i++) {
    z[i] = x[i];
  }

  MCP_Destroy(m);
  destroy();

  Options_Destroy(o);
#ifdef PATHFERRIS_LOG_IN_FILE
  fclose(f);
#endif
  return;
}
#else

void SimpleLCP(int variables, int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
               double *lb, double *ub, MCP_Termination *status, double *z) {
  ;
}
void printLCP(int variables, int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
              double *lb, double *ub) {
  ;
}

int nbNonNulElems(int n, double *M, double tol) { return 0; }
void FortranToPathSparse(int n, double *M, double tol, int *m_i, int *m_j, double *m_ij) { ; }
void ABCDtoM(int n, int m, double *A, double *B, double *C, double *D, double *a, double *b,
             double *M, double *q) {
  ;
}

#endif
