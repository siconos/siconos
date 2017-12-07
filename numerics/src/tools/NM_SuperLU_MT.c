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

#ifdef WITH_SUPERLU_MT

#include <slu_mt_ddefs.h>

/** \struct NM_SuperLU_MT_WS NumericsMatrix_internal.h
 * Structure for holding the data SuperLU needs
 */
struct NM_SuperLU_MT_WS {
  SuperMatrix* L;
  SuperMatrix* U;
  int_t* perm_r;
  int_t* perm_c;
  superlumt_options_t* options;
  Gstat_t* stat;
  SuperMatrix* A;
};

NM_SuperLU_MT_WS* NM_SuperLU_MT_factorize(NumericsMatrix* A)
{
  SuperMatrix SAC;

  int status;

  NumericsSparseLinearSolverParams* params = NM_linearSolverParams(A);

  if (params->solver_data)
  {
    return (NM_SuperLU_MT_WS*) params->solver_data;
  }

  params->solver_data = calloc(1, sizeof(NM_SuperLU_MT_WS));
  NM_SuperLU_MT_WS* superlu_mt_ws = (NM_SuperLU_MT_WS*) params->solver_data;

  if (!superlu_mt_ws->options) superlu_mt_ws->options = (superlumt_options_t*)SUPERLU_MALLOC(sizeof(superlumt_options_t));
  if (!superlu_mt_ws->stat) superlu_mt_ws->stat = (Gstat_t*)SUPERLU_MALLOC(sizeof(Gstat_t));

  CSparseMatrix* C = NM_csc(A);

  superlu_mt_ws->L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
  superlu_mt_ws->U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
  if ( !(superlu_mt_ws->perm_r = intMalloc(C->m)) ) SUPERLU_ABORT("Malloc fails for perm_r[].");
  if ( !(superlu_mt_ws->perm_c = intMalloc(C->n)) ) SUPERLU_ABORT("Malloc fails for perm_c[].");

  /* Symbolic part  */
  int_t* indices;
  int_t* pointers;
  size_t nnz = NM_sparse_nnz(C);
  assert(C->n > 0);

  if (sizeof(*C->i) != sizeof(*indices))
  {
    int_t* iWork = (int_t*)NM_iWork(A, (size_t)(nnz + C->n) + 1, sizeof(int_t));
    indices = iWork;
    pointers = &iWork[nnz];

    for (size_t i = 0; i < nnz  ; ++i) indices[i] = (int_t)C->i[i];
    for (size_t i = 0; i <= (size_t)C->n; ++i) pointers[i] = (int_t)C->p[i];
  }
  else
  {
    indices = (int_t*)C->i;
    pointers = (int_t*)C->p;
  }

  superlu_mt_ws->A = (SuperMatrix*)malloc(sizeof(SuperMatrix));
  dCreate_CompCol_Matrix(superlu_mt_ws->A, C->m, C->n, nnz, C->x, indices, pointers, SLU_NC, SLU_D, SLU_GE);

/* TODO SuperLU_MT_PIVOT_TOLERANCE, SuperLU_MT_ORDERING, SuperLU_MT_SCALE
 * SuperLU_MT_DROPTOL, SuperLU_MT_STRATEGY, SuperLU_MT_IRSTEP*/

  /* Numerical part */
  int panel_size = sp_ienv(1);
  int relax = sp_ienv(2);
  double diag_pivot_thresh = 1.0;
  yes_no_t usepr = NO;
  double drop_tol = 0.0;
  double* work = NULL;
  int lwork = 0;

  // XXX fix that --xhub
  int_t nprocs = 0;
  char* nprocs_str = getenv("SICONOS_SPARSE_SOLVER_PROCS");
  if (nprocs_str)
  {
    nprocs = atoi(nprocs_str);
  }
  else
  {
    nprocs = 1;
  }

  /*  XXX allow the user to choose that and resuse perm_c */
  int_t permc_spec = 3;
  get_perm_c(permc_spec, superlu_mt_ws->A, superlu_mt_ws->perm_c);

  /* ------------------------------------------------------------
     Allocate storage and initialize statistics variables. 
     ------------------------------------------------------------*/
  StatAlloc(C->n, nprocs, panel_size, relax, superlu_mt_ws->stat);
  StatInit(C->n, nprocs, superlu_mt_ws->stat);

  pdgstrf_init(nprocs, EQUILIBRATE, NOTRANS, NO, panel_size, relax,
    diag_pivot_thresh, usepr, drop_tol, superlu_mt_ws->perm_c, superlu_mt_ws->perm_r,
    work, lwork, superlu_mt_ws->A, &SAC, superlu_mt_ws->options, superlu_mt_ws->stat);

  if (verbose > 1)
    superlu_mt_ws->options->PrintStat = YES;

  pdgstrf(superlu_mt_ws->options, &SAC, superlu_mt_ws->perm_r, superlu_mt_ws->L, superlu_mt_ws->U, superlu_mt_ws->stat, &status);

  if (status)
  {
    fprintf(stderr, "dgstrf() error returns INFO= %d\n", status);
    return NULL;
  }

  /* ------------------------------------------------------------
     Deallocate storage after factorization.
     ------------------------------------------------------------*/
  pxgstrf_finalize(superlu_mt_ws->options, &SAC);
  SUPERLU_FREE(superlu_mt_ws->options);
  superlu_mt_ws->options = NULL;

  return superlu_mt_ws;
}

int NM_SuperLU_MT_solve(NumericsMatrix* A, double* b, NM_SuperLU_MT_WS* superlu_mt_ws)
{
  SuperMatrix B, X;

  int status;
  equed_t equed = NOEQUIL;

  CSparseMatrix* C = NM_csc(A);

  bool iter_refinement = false;

  double* bcopy = NULL;
  if (iter_refinement)
  {
    bcopy = (double*)malloc(C->n * sizeof(double));
    memcpy(bcopy, b, C->n * sizeof(double));
    dCreate_Dense_Matrix(&B, C->n, 1, bcopy, C->n, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, C->n, 1, b, C->n, SLU_DN, SLU_D, SLU_GE);

    dgstrs(NOTRANS, superlu_mt_ws->L, superlu_mt_ws->U, superlu_mt_ws->perm_r, superlu_mt_ws->perm_c, &X, superlu_mt_ws->stat, &status);
  }
  else
  {
    dCreate_Dense_Matrix(&B, C->n, 1, b, C->n, SLU_DN, SLU_D, SLU_GE);

    dgstrs(NOTRANS, superlu_mt_ws->L, superlu_mt_ws->U, superlu_mt_ws->perm_r, superlu_mt_ws->perm_c, &B, superlu_mt_ws->stat, &status);
  }


#if 0
  if (superlu_mt_ws->options->fact == DOFACT || superlu_mt_ws->options->fact == EQUILIBRATE)
  {
    equed = NOEQUIL;
  }
#endif

  if (iter_refinement)
  {
    double ferr;
    double berr;

    if (!status)
      dgsrfs(NOTRANS, superlu_mt_ws->A, superlu_mt_ws->L, superlu_mt_ws->U, superlu_mt_ws->perm_r, superlu_mt_ws->perm_c, equed, NULL, NULL, &B, &X, &ferr, &berr, superlu_mt_ws->stat, &status);

    Destroy_SuperMatrix_Store(&X);
    free(bcopy);
  }

  Destroy_SuperMatrix_Store(&B);

  return status;
}

void NM_SuperLU_MT_free(void* p)
{
  assert(p);
  NumericsSparseLinearSolverParams* params = (NumericsSparseLinearSolverParams*) p;
  assert(params);
  NM_SuperLU_MT_WS* superlu_mt_ws = (NM_SuperLU_MT_WS*) params->solver_data;
  assert(superlu_mt_ws);

  SUPERLU_FREE (superlu_mt_ws->perm_r);
  SUPERLU_FREE (superlu_mt_ws->perm_c);
  Destroy_SuperNode_Matrix(superlu_mt_ws->L);
  Destroy_CompCol_Matrix(superlu_mt_ws->U);
  SUPERLU_FREE (superlu_mt_ws->L);
  SUPERLU_FREE (superlu_mt_ws->U);

  superlu_mt_ws->perm_r = NULL;
  superlu_mt_ws->perm_c = NULL;
  superlu_mt_ws->L = NULL;
  superlu_mt_ws->U = NULL;

  if (superlu_mt_ws->A)
  {
    Destroy_SuperMatrix_Store(superlu_mt_ws->A);
    SUPERLU_FREE(superlu_mt_ws->A);
    superlu_mt_ws->A = NULL;
  }

  if (superlu_mt_ws->stat)
  {
    StatFree(superlu_mt_ws->stat);
    SUPERLU_FREE(superlu_mt_ws->stat);
    superlu_mt_ws->stat = NULL;
  }

  if (superlu_mt_ws->options)
  {
    SUPERLU_FREE(superlu_mt_ws->options);
    superlu_mt_ws->options = NULL;
  }

  /* Here we free superlu_mt_ws ...  */
  free(superlu_mt_ws);
  params->solver_data = NULL;

}

void NM_SuperLU_MT_extra_display(NM_SuperLU_MT_WS* superlu_mt_ws)
{
  if (superlu_mt_ws->stat)
  {
    PrintStat(superlu_mt_ws->stat);
  }
#if 0
  if (verbose > 2)
  {
    SuperLU_MT_FN(report_info) (superlu_mt_ws->control, superlu_mt_ws->info);

    if (verbose > 3)
    {
      SuperLU_MT_FN(report_control) (superlu_mt_ws->control);
    }
  }
  else if (verbose > 1)
  {
    if (superlu_mt_ws->control[SuperLU_MT_IRSTEP] > 0)
    {
      printf("SuperLU_MT : backward error estimate omega1 %g\n", superlu_mt_ws->info[SuperLU_MT_OMEGA1]);
      printf("SuperLU_MT : backward error estimate omega2 %g\n", superlu_mt_ws->info[SuperLU_MT_OMEGA2]);
    }
    printf("SuperLU_MT : solve FLOPS %g\n", superlu_mt_ws->info[SuperLU_MT_SOLVE_FLOPS]);
    printf("SuperLU_MT : solve time %g\n", superlu_mt_ws->info[SuperLU_MT_SOLVE_TIME]);
    printf("SuperLU_MT : wall time %g\n", superlu_mt_ws->info[SuperLU_MT_SOLVE_WALLTIME]);

  }
#endif
}

#endif /* WITH_SuperLU_MT */
