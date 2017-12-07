/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "NumericsMatrix_internal.h"
#include "NumericsSparseMatrix.h"
#include "debug.h"
#include "numerics_verbose.h"

#ifdef WITH_SUPERLU

#include <slu_ddefs.h>

/** \struct NM_SuperLU_WS NumericsMatrix_internal.h
 * Structure for holding the data SuperLU needs
 */
struct NM_SuperLU_WS {
  SuperMatrix* L;
  SuperMatrix* U;
  int_t* perm_r;
  int_t* perm_c;
  superlu_options_t* options;
#ifdef SUPERLU_MAJOR_VERSION
  GlobalLU_t* Glu;
#endif
};

NM_SuperLU_WS* NM_SuperLU_factorize(NumericsMatrix* A)
{
  SuperMatrix SA, SAC;
  SuperLUStat_t stat;
  int_t *etree;

  int status;

  NumericsSparseLinearSolverParams* params = NM_linearSolverParams(A);

  if (params->solver_data)
  {
    return (NM_SuperLU_WS*) params->solver_data;
  }

  params->solver_data = calloc(1, sizeof(NM_SuperLU_WS));
  NM_SuperLU_WS* superlu_ws = (NM_SuperLU_WS*) params->solver_data;

  if (!superlu_ws->options) superlu_ws->options = (superlu_options_t*)malloc(sizeof(superlu_options_t));

#ifdef SUPERLU_MAJOR_VERSION
  if (!superlu_ws->Glu) superlu_ws->Glu = (GlobalLU_t*)malloc(sizeof(GlobalLU_t));
#endif

  set_default_options(superlu_ws->options);

  if (verbose > 1)
    superlu_ws->options->PrintStat = YES;
/* TODO SuperLU_PIVOT_TOLERANCE, SuperLU_ORDERING, SuperLU_SCALE
 * SuperLU_DROPTOL, SuperLU_STRATEGY, SuperLU_IRSTEP*/

  StatInit(&stat);

  CSparseMatrix* C = NM_csc(A);

  superlu_ws->L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
  superlu_ws->U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
  if ( !(superlu_ws->perm_r = intMalloc(C->m)) ) ABORT("Malloc fails for perm_r[].");
  if ( !(superlu_ws->perm_c = intMalloc(C->n)) ) ABORT("Malloc fails for perm_c[].");
  if ( !(etree = intMalloc(C->n)) ) ABORT("Malloc fails for etree[].");

  /* Symbolic part  */
  int_t* indices;
  int_t* pointers;
  size_t nnz = NM_sparse_nnz(C);

  if (sizeof(*C->i) != sizeof(*indices))
  {
    int_t* iWork = (int_t*)NM_iWork(A, (size_t)(nnz + C->n) + 1, sizeof(int_t));
    indices = iWork;
    pointers = &iWork[nnz];

    for (size_t i = 0; i < nnz  ; ++i) indices[i] = (int_t)C->i[i];
    for (size_t i = 0; i <= C->n; ++i) pointers[i] = (int_t)C->p[i];
  }
  else
  {
    indices = (int_t*)C->i;
    pointers = (int_t*)C->p;
  }

  dCreate_CompCol_Matrix(&SA, C->m, C->n, nnz, C->x, indices, pointers, SLU_NC, SLU_D, SLU_GE);

  int permc_spec = 3;
  get_perm_c(permc_spec, &SA, superlu_ws->perm_c);

  sp_preorder(superlu_ws->options, &SA, superlu_ws->perm_c, etree, &SAC);

  /* Numerical part */
  int panel_size = sp_ienv(1);
  int relax = sp_ienv(2);
  double drop_tol = 0.0;

  /* Versions 4.x and earlier do not include a #define'd version numbers  */
#ifndef SUPERLU_MAJOR_VERSION
  dgstrf(superlu_ws->options, &SAC, drop_tol, relax, &panel_size, etree, NULL, 0, superlu_ws->perm_c, superlu_ws->perm_r, superlu_ws->L, superlu_ws->U, &stat, &status);
#else
  dgstrf(superlu_ws->options, &SAC, relax, panel_size, etree, NULL, 0, superlu_ws->perm_c, superlu_ws->perm_r, superlu_ws->L, superlu_ws->U, superlu_ws->Glu, &stat, &status);
#endif

  if (status)
  {
    fprintf(stderr, "dgstrf() error returns INFO= %d\n", status);
    return NULL;
  }

  SUPERLU_FREE(etree);
  Destroy_SuperMatrix_Store(&SA);
  Destroy_CompCol_Permuted(&SAC);
  StatFree(&stat);

  return superlu_ws;
}

int NM_SuperLU_solve(NumericsMatrix* A, double* b, NM_SuperLU_WS* superlu_ws)
{
  SuperMatrix B;
  SuperLUStat_t stat;

  int status;

  CSparseMatrix* C = NM_csc(A);

  dCreate_Dense_Matrix(&B, C->n, 1, b, C->n, SLU_DN, SLU_D, SLU_GE);

  StatInit(&stat);

  dgstrs(NOTRANS, superlu_ws->L, superlu_ws->U, superlu_ws->perm_c, superlu_ws->perm_r, &B, &stat, &status);

  Destroy_SuperMatrix_Store(&B);
  StatFree(&stat);

  return status;
}

void NM_SuperLU_free(void* p)
{
  assert(p);
  NumericsSparseLinearSolverParams* params = (NumericsSparseLinearSolverParams*) p;
  assert(params);
  NM_SuperLU_WS* superlu_ws = (NM_SuperLU_WS*) params->solver_data;
  assert(superlu_ws);

  SUPERLU_FREE (superlu_ws->perm_r);
  SUPERLU_FREE (superlu_ws->perm_c);
  Destroy_SuperNode_Matrix(superlu_ws->L);
  Destroy_CompCol_Matrix(superlu_ws->U);
  SUPERLU_FREE (superlu_ws->L);
  SUPERLU_FREE (superlu_ws->U);

  superlu_ws->perm_r = NULL;
  superlu_ws->perm_c = NULL;
  superlu_ws->L = NULL;
  superlu_ws->U = NULL;

  if (superlu_ws->options)
  {
    free(superlu_ws->options);
    superlu_ws->options = NULL;
  }

#ifdef SUPERLU_MAJOR_VERSION
  if (superlu_ws->Glu)
  {
    free(superlu_ws->Glu);
    superlu_ws->Glu = NULL;
  }
#endif
  /* Here we free superlu_ws ...  */
  free(superlu_ws);
  params->solver_data = NULL;

}

void NM_SuperLU_extra_display(NM_SuperLU_WS* superlu_ws)
{
#if 0
  if (verbose > 2)
  {
    SuperLU_FN(report_info) (superlu_ws->control, superlu_ws->info);

    if (verbose > 3)
    {
      SuperLU_FN(report_control) (superlu_ws->control);
    }
  }
  else if (verbose > 1)
  {
    if (superlu_ws->control[SuperLU_IRSTEP] > 0)
    {
      printf("SuperLU : backward error estimate omega1 %g\n", superlu_ws->info[SuperLU_OMEGA1]);
      printf("SuperLU : backward error estimate omega2 %g\n", superlu_ws->info[SuperLU_OMEGA2]);
    }
    printf("SuperLU : solve FLOPS %g\n", superlu_ws->info[SuperLU_SOLVE_FLOPS]);
    printf("SuperLU : solve time %g\n", superlu_ws->info[SuperLU_SOLVE_TIME]);
    printf("SuperLU : wall time %g\n", superlu_ws->info[SuperLU_SOLVE_WALLTIME]);

  }
#endif
}

#endif /* WITH_SuperLU */
