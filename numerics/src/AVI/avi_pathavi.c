/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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


#include <string.h>

#include "SparseMatrix_internal.h"
#include "AffineVariationalInequalities.h"
#include "AVI_Solvers.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"
#include "SiconosSets.h"

#ifdef HAVE_PATHVI

#include "PATHVI_helpers.h"
#include "PATHVI_SDK/include/vi_consts.h"
#include "PATHVI_SDK/include/csc_matrix.h"
#include "PATHVI_SDK/include/vi_scheduler.h"
#include "PATHVI_SDK/include/vi_solver_pathvi.h"
#include "PATHVI_SDK/include/printv.h"

static void pathvi_csc_transfert(struct csc_matrix *primjac, CSparseMatrix* M)
{
  size_t n = primjac->n;
  size_t nnz = primjac->nnz;

  /*  Easy copy */
  memcpy(primjac->x, M->x, nnz * sizeof(double));

  /*  we may have to change the behavior based on the type used for integers */
  if (sizeof(primjac->i) != sizeof(CS_INT))
  {
    for (size_t i = 0; i < nnz; ++i)
    {
      primjac->i[i] = (PATHVI_INDX_TYPE)M->i[i];
    }

    for (size_t i = 0; i <= n; ++i)
    {
      primjac->j[i] = (PATHVI_INDX_TYPE)M->p[i];
    }
  }
  else
  {
    memcpy(primjac->i, M->i, nnz * sizeof(CS_INT));
    memcpy(primjac->j, M->p, n+1 * sizeof(CS_INT));
  }

}


static int pathvi_evaluate_function(struct vi_desc *desc, double *primvar, double *primfunc)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  AffineVariationalInequalities* AVI = (AffineVariationalInequalities*) env->problem;

  memcpy(primfunc, AVI->q, AVI->size * sizeof(double));
  NM_gemv(1.0, AVI->M, primvar, 1., primfunc);

  return 0;
}

static int pathvi_evaluate_jacobian(struct vi_desc *desc, double *primvar, double *primfunc, struct csc_matrix *primjac)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  AffineVariationalInequalities* AVI = (AffineVariationalInequalities*) env->problem;

  size_t n = AVI->size;

  switch (AVI->M->storageType)
  {
  case NM_DENSE:
  {
    size_t nnz = AVI->M->size0 * AVI->M->size1;

    break;
  }
  case NM_SPARSE:
  {
    CSparseMatrix* M = NM_csc(AVI->M);
    CS_INT nnz = M->p[n];

    /* check dimenstions */
    if (M->n > primjac->max_n) { primjac->max_n = M->n; primjac->j = (PATHVI_INDX_TYPE*)realloc(primjac->j, (M->n+1) * sizeof(primjac->j)); }

    if (nnz > primjac->max_nnz)
    {
      primjac->max_nnz = nnz;
      primjac->i = (PATHVI_INDX_TYPE*)realloc(primjac->i, nnz * sizeof(primjac->i));
      primjac->x = (double*)realloc(primjac->x, nnz * sizeof(double));
    }

    /*  Update the size */
    csc_size(primjac, M->m, M->n, nnz);

    break;
  }
  default:
    numerics_error_nonfatal("evaluate_jacobian", "Unsupported matrix storage type %d\n", AVI->M->storageType);
    return -1;
  }

  return 0;
}


static int pathavi_get_jacobian_nnz(struct vi_desc *desc, int *nnz)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  *nnz = (int)NM_nnz(((AffineVariationalInequalities*) env->problem)->M);
  return 0;
}

static int pathavi_get_jacobian_structure(struct vi_desc *desc, struct csc_matrix *primjac)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  NumericsMatrix *M = ((AffineVariationalInequalities*) env->problem)->M;

  CSparseMatrix* Mcsc = NM_csc(M);

  if ((primjac->struct_mask & CSC_STRUCT_FIXED) &&
      !(primjac->struct_mask & CSC_STRUCT_FILLED)) {
    for (size_t i = 0; i <= env->n; ++i) {
      primjac->j[i] = Mcsc->p[i];
    }

    for (size_t i = 0; i < Mcsc->p[env->n]; ++i) {
      primjac->i[i] = Mcsc->i[i];
    }

    csc_size(primjac, env->n, env->n, Mcsc->p[env->n]);
    primjac->struct_mask |= CSC_STRUCT_FILLED;
  }

  return 0;
}

int avi_pathavi(AffineVariationalInequalities* problem, double *z, double *w, SolverOptions* options)
{
  int info = 0;
  bool use_scheduler = false;

  int nb_cstr = 0;
  int nnz_H = 0;
  double* lambda = NULL;

  if (problem->poly.set->id == SICONOS_SET_POLYHEDRON_UNIFIED)
  {
    nb_cstr = problem->poly.unif->A->size0;
    nnz_H = (int)NM_nnz(problem->poly.unif->A);
    lambda = (double*) calloc(nb_cstr, sizeof(double));
  }
  else
  {
    numerics_error_nonfatal("avi_pathavi", "unsupported set type %d", problem->poly.set->id);
    return -1;
  }

  SN_generic_pathvi_env env = {
    .problem = problem,
    .n = problem->size,
    .m = nb_cstr,
    .z = z,
    .F = w,
    .lambda = lambda
  };

  struct vi_desc_operations vi_ops = {
    .get_row_name           = &pathvi_get_row_name,
    .get_col_name           = &pathvi_get_col_name,
    .get_primvar            = &pathvi_get_z,
    .set_primvar            = &pathvi_set_z,
    .get_prim_marginal      = &pathvi_get_F,
    .set_prim_marginal      = &pathvi_set_F,
    .get_dualvar            = &pathvi_get_lambda,
    .set_dualvar            = &pathvi_set_lambda,
    .evaluate_function      = &pathvi_evaluate_function,
    .evaluate_jacobian      = &pathvi_evaluate_jacobian,
    .get_jacobian_nnz       = &pathavi_get_jacobian_nnz,
    .get_jacobian_structure = &pathavi_get_jacobian_structure
  };

  struct vi_desc * pathvi_obj = vi_desc_create(nb_cstr, problem->size, NM_nnz(problem->M), nnz_H, &env, &vi_ops);
  pathvi_obj->nlflag = 0;

  if (problem->lb)
  {
    memcpy(pathvi_obj->primlower, problem->lb, env.n*sizeof(double));
  }
  else
  {
    for (size_t i = 0; i < env.n; ++i) pathvi_obj->primlower[i] = -INF;
  }

  if (problem->ub)
  {
    memcpy(pathvi_obj->primupper, problem->ub, env.n*sizeof(double));
  }
  else
  {
    for (size_t i = 0; i < env.n; ++i) pathvi_obj->primupper[i] = INF;
  }

  if (problem->poly.set)
  {
    if (problem->poly.set->id == SICONOS_SET_POLYHEDRON_UNIFIED)
    {
      polyhedron_unified* p = problem->poly.unif;
      pathvi_csc_transfert(pathvi_obj->A, NM_csc(p->A));

      memcpy(pathvi_obj->b, p->b, env.m*sizeof(double));

      for (size_t i = 0; i < env.m; ++i)
      {
        switch (p->type[i])
        {
        case SICONOS_LE:
        {
          pathvi_obj->sense[i] = 'L';
          pathvi_obj->duallower[i] = -INF;
          pathvi_obj->dualupper[i] = 0.;
          break;
        }
        case SICONOS_EQ:
        {
          pathvi_obj->sense[i] = 'E';
          pathvi_obj->duallower[i] = -INF;
          pathvi_obj->dualupper[i] = INF;
          break;
        }
        case SICONOS_GE:
        {
          pathvi_obj->sense[i] = 'G';
          pathvi_obj->duallower[i] = 0.;
          pathvi_obj->dualupper[i] = INF;
          break;
        }
        default:
        {
          numerics_error_nonfatal("avi_pathavi", "unsupported constraint type %c", p->type[i]);
          info = EINVAL;
          goto exit;
        }
        }
      }

      // We should not do this
      csc_transpose(pathvi_obj->negAt, pathvi_obj->A);
      csc_negate(pathvi_obj->negAt);
    }
    else
    {
      numerics_error_nonfatal("avi_pathavi", "unsupported set type %d", problem->poly.set->id);
      return -1;
    }

  }

  // bug youngdae
  pathvi_register_options(pathvi_obj->opt);

  option_set_d(pathvi_obj->opt, "convergence_tolerance", options->dparam[0]);

  struct printv_operations printv_ops = {
    .print       = &pathvi_print
  };

  set_printv_operations(&printv_ops);


  // Solve the problem
  if (use_scheduler)
  {
    struct vi_scheduler *sched = vi_scheduler_create(VI_SOLVER_PATHVI, pathvi_obj);
    int sinfo = vi_scheduler_run(sched);
    info = vi_scheduler_get_status(sched);
    vi_scheduler_free(&sched);
  }
  else
  {
    struct vi_solver *avi = vi_solver_create(VI_SOLVER_PATHVI, pathvi_obj);
    int sinfo = avi->ops->solve(avi);
    info = avi->status;
    vi_solver_free(&avi);
  }

exit:

  return info;
}

#else

int avi_pathavi(AffineVariationalInequalities* problem, double *z, double *w, SolverOptions* options)
{
  numerics_error_nonfatal("avi_pathavi", "PATHVI not enabled");
  return -1;
}


#endif /*  HAVE_PATHVI */
