/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include <assert.h>                        // for assert
#include "Friction_cst.h"                  // for SICONOS_GLOBAL_FRICTION_3D...
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "NumericsSparseMatrix.h"                // for NSM_TRIPLET ...
#include "numerics_verbose.h"              // for numerics_printf_verbose
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_dscal
#include "CSparseMatrix_internal.h"                // for NSM_TRIPLET ...
#include "gfc3d_balancing.h"


/* #define DEWBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"                         // for DEBUG_EXPR
#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
#endif

void gfc3d_rescaling(
  GlobalFrictionContactProblem* problem,
  double alpha,
  double beta,
  double gamma)
{
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;

  /* scaling of M */
  NM_scal(alpha*gamma*gamma, problem->M);
  /* scaling of H */
  NM_scal(beta*gamma, problem->H);
  /* scaling of q */
  cblas_dscal(n,alpha*gamma,problem->q,1);
  /* scaling of b */
  cblas_dscal(m,beta,problem->b,1);

}

void gfc3d_balancing_MHHT(
  GlobalFrictionContactProblem* problem,
  BalancingMatrices * B_for_M,
  BalancingMatrices * B_for_H)
{
  assert(B_for_M);
  assert(B_for_H);

  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;

  /* scaling of M */
  NumericsMatrix* M_tmp = NM_create(NM_SPARSE, n, n);
  NM_triplet_alloc(M_tmp, n);
  NM_gemm(1.0,problem->M,B_for_M->D2,0.0, M_tmp) ; // M * D2M
  NM_gemm(1.0,B_for_M->D2,M_tmp, 0.0, problem->M); // D1M *M  * D2M

  /* scaling of H */
  /* Warning the matrix H must be scaled such
   * that the cone constraint is respected */
  for (int contact =0, i=0; contact < nc; contact++, i++, i++, i++)
  {
    /* A choice among other */
    NM_triplet(B_for_H->D2)->x[i+1] = NM_triplet(B_for_H->D2)->x[i];
    NM_triplet(B_for_H->D2)->x[i+2] = NM_triplet(B_for_H->D2)->x[i];
  }
  NumericsMatrix* H_tmp = NM_create(NM_SPARSE, n, m);
  NM_triplet_alloc(H_tmp, n);
  NM_gemm(1.0,problem->H,B_for_H->D2,0.0, H_tmp); // H * D2H
  NM_gemm(1.0,B_for_M->D1,H_tmp, 0.0, problem->H); // D1M * H* D2H
  //NM_gemm(1.0,B_for_M->D1,problem->M, 0.0, M_tmp);
  //NM_copy(M_tmp, problem->M);

  /* scaling of q */
  double * q_tmp = (double *) malloc(n*sizeof(double));
  NM_gemv(1.0, B_for_M->D2, problem->q, 0.0, q_tmp);
  cblas_dcopy(n, q_tmp, 1, problem->q, 1);

  /* scaling of b */
  double * b_tmp = (double *) malloc(m*sizeof(double));
  NM_gemv(1.0, B_for_H->D2, problem->b, 0.0, b_tmp);
  cblas_dcopy(m, b_tmp, 1, problem->b, 1);


  NM_clear(M_tmp);
  free(M_tmp);
  NM_clear(H_tmp);
  free(H_tmp);

  free(q_tmp);
  free(b_tmp);
}


void gfc3d_balancing_M(
  GlobalFrictionContactProblem* problem,
  BalancingMatrices * B_for_M)
{
  assert(B_for_M);

  NM_compute_balancing_matrices(problem->M, 1e-03, 100, B_for_M);

  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;
  /* scaling of M */
  /* NM_scal(alpha*gamma*gamma, problem->M); */
  NumericsMatrix* M_tmp = NM_create(NM_SPARSE, n, n);
  NM_triplet_alloc(M_tmp, n);
  NM_gemm(1.0,problem->M,B_for_M->D2,0.0, M_tmp);
  NM_gemm(1.0,B_for_M->D1,M_tmp, 0.0, problem->M);

  /* scaling of H */
  /* NM_scal(beta*gamma, problem->H);*/
  /* Warning the matrix H must be scaled such
   * that the cone constraint is respected */
  NumericsMatrix* H_tmp = NM_create(NM_SPARSE, n, m);
  NM_triplet_alloc(H_tmp, n);
  NM_gemm(1.0, B_for_M->D2, problem->H, 0.0, H_tmp);
  NM_copy(H_tmp, problem->H);

  /* scaling of q */
  /* cblas_dscal(n,alpha*gamma,problem->q,1); */
  double * q_tmp = (double *) malloc(n*sizeof(double));
  NM_gemv(1.0, B_for_M->D2, problem->q, 0.0, q_tmp);
  cblas_dcopy(n, q_tmp, 1, problem->q, 1);

  /* scaling of b */
  /* cblas_dscal(m,beta,problem->b,1); */

  NM_clear(M_tmp);
  free(M_tmp);
  NM_clear(H_tmp);
  free(H_tmp);

  free(q_tmp);
  //free(b_tmp);
}


GlobalFrictionContactProblem*  gfc3d_balancing_problem(GlobalFrictionContactProblem* problem,
                                                               SolverOptions* options)
{
  GlobalFrictionContactProblem * rescaled_problem = NULL;
  GlobalFrictionContactProblem_balancing_data  *data = NULL;

  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
  {
    rescaled_problem =  globalFrictionContact_copy(problem);
    data = gfc3d_balancing_data_new();
    rescaled_problem->env = (void*) data;
    data->original_problem = problem;
  }
  else
  {
    return problem;
  }

  size_t nc = problem->numberOfContacts;
  size_t n = problem->M->size0;
  size_t m = 3 * nc;
  // double* q = problem->q;
  // double* b = problem->b;
  // double* mu = problem->mu;

  NumericsMatrix *M = problem->M;
  NumericsMatrix *H = problem->H;


  data->original_problem = problem;

  double alpha_r=0.0, beta_r=0.0;
  // BalancingMatrices * B_for_M = NULL;
  // BalancingMatrices * B_for_H = NULL;

  /* NumericsMatrix *Htrans =  NM_transpose(H); */
  /* /\* Compute M + rho H H^T (storage in W)*\/ */
  /* NumericsMatrix *W = NM_create(NM_SPARSE,n,n); */
  /* NM_triplet_alloc(W, n); */
  /* W->matrix2->origin = NSM_TRIPLET; */


  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_SCALAR)
  {
    alpha_r = NM_norm_inf(M);
    beta_r = NM_norm_inf(H);
    numerics_printf_verbose(1,"---- GFC3D - BALANCING - Scalar rescaling of the problem");
    numerics_printf_verbose(1,"---- GFC3D - BALANCING - alpha_r = %e\t beta_r= %e\n", alpha_r, beta_r);

    gfc3d_rescaling(rescaled_problem, 1./alpha_r, 1.0/beta_r, 1.0);

    data->alpha= 1.0/alpha_r;
    data->beta=  1.0/beta_r;
    data->gamma= 1.0 ;
  }
  else if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_M)
  {
    numerics_printf_verbose(1,"---- GFC3D - BALANCING - Rescaling of the problem by balancing M");
    data->B_for_M  = NM_BalancingMatrices_new(problem->M);
    NM_compute_balancing_matrices(problem->M, 1e-2, 5, data->B_for_M);
    gfc3d_balancing_M(rescaled_problem, data->B_for_M);
  }
  /* else if (options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_H) */
  /* { */
  /*   numerics_printf_verbose(1,"---- GFC3D - BALANCING - Rescaling of the problem by balancing H"); */
  /*   data->B_for_H  = NM_BalancingMatrices_new(problem->H); */
  /*   globalFrictionContact_balancing_H(rescaled_problem, data->B_for_H); */
  /* } */
  else if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT)
  {

    numerics_printf_verbose(1,"---- GFC3D - BALANCING - Rescaling of the problem by balancing MHHT");

    /****** compute balancing matrices **************************************/
    NumericsMatrix * MHHT = NM_create(NM_SPARSE, n+m, n+m);
    NM_triplet_alloc(MHHT, n+m);
    MHHT->matrix2->origin = NSM_TRIPLET;
    NM_insert(MHHT, problem->M, 0, 0);
    NM_insert(MHHT, problem->H, 0, n);
    NumericsMatrix *HT =  NM_transpose(H);
    NM_insert(MHHT, HT, n, 0);

    //NM_display(MHHT);
    BalancingMatrices * B_for_MHHT = NM_BalancingMatrices_new(MHHT);
    NM_compute_balancing_matrices(MHHT, 1e-2, 5, B_for_MHHT);
    DEBUG_EXPR(NM_display(B_for_MHHT->D1););
    DEBUG_EXPR(NM_display(B_for_MHHT->D2););

    /* to simplify the use of the balancing matrices, we split it */
    data->B_for_M = NM_BalancingMatrices_new(problem->M); // lazy mode
    data->B_for_H = NM_BalancingMatrices_new(problem->H);

    for(size_t i =0; i < n ; i++)
    {
      NM_triplet(data->B_for_M->D1)->x[i] = NM_triplet(B_for_MHHT->D1)->x[i]; // D1M
      NM_triplet(data->B_for_M->D2)->x[i] = NM_triplet(B_for_MHHT->D2)->x[i]; // D2M
    }

    for(size_t i =0; i < n ; i++)
    {
      NM_triplet(data->B_for_H->D1)->x[i] = NM_triplet(B_for_MHHT->D1)->x[i]; // D1M
    }

    for(size_t i =0; i < m ; i++)
    {
      NM_triplet(data->B_for_H->D2)->x[i] = NM_triplet(B_for_MHHT->D2)->x[i+n]; //D2H
    }



    /****** balanced problem          ***************************************/
    gfc3d_balancing_MHHT(rescaled_problem, data->B_for_M, data->B_for_H);
  }
  else
  {
    numerics_printf_verbose(1,"---- GFC3D - BALANCING - No rescaling of the problem");
  }

  return rescaled_problem;
}


void gfc3d_balancing_go_to_balanced_variables(GlobalFrictionContactProblem* balanced_problem,
                                SolverOptions* options,
                                double *r, double *u, double* v)
{
  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
  {
    size_t nc = balanced_problem->numberOfContacts;
    size_t n = balanced_problem->M->size0;
    size_t m = 3 * nc;

    GlobalFrictionContactProblem_balancing_data  *data = (GlobalFrictionContactProblem_balancing_data * ) balanced_problem->env;
    assert(data);

    /* rescale */
    if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_SCALAR)
    {
      cblas_dscal(m, data->alpha/data->beta, r, 1);
      cblas_dscal(m, data->beta, u, 1);
      cblas_dscal(n, 1.0/data->gamma, v, 1);
    }
    else if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_M)
    {
      for(size_t i =0; i < n ; i++)
      {
        v[i] = v[i]/NM_triplet(data->B_for_M->D2)->x[i];
      }
      //nothing for u and r
    }
    else if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT)
    {
      for(size_t i =0; i < n ; i++)
      {
        v[i] = v[i]/NM_triplet(data->B_for_M->D2)->x[i];
      }
      for(size_t i =0; i < m ; i++)
      {
        r[i] = r[i]/NM_triplet(data->B_for_H->D2)->x[i];
        u[i] = u[i]*NM_triplet(data->B_for_H->D2)->x[i];
      }
    }
    else
      numerics_printf_verbose(1,"---- GFC3D - BALANCING - rescaling type is not implemented");

  }
  //else continue;

}
void gfc3d_balancing_back_to_original_variables(GlobalFrictionContactProblem* balanced_problem,
                                    SolverOptions* options,
                                    double *r, double *u, double *v)
{

  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
  {
    size_t nc = balanced_problem->numberOfContacts;
    size_t n = balanced_problem->M->size0;
    size_t m = 3 * nc;

    GlobalFrictionContactProblem_balancing_data  *data = (GlobalFrictionContactProblem_balancing_data * ) balanced_problem->env;
    assert(data);

    /* rescale */
    if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_SCALAR)
    {
      cblas_dscal(m, data->beta/data->alpha, r, 1);
      cblas_dscal(m, 1.0/data->beta, u, 1);
      cblas_dscal(n, data->gamma, v, 1);
    }
    else if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_M)
    {
      for(size_t i =0; i < n ; i++)
      {
        v[i] = v[i]*NM_triplet(data->B_for_M->D2)->x[i];
      }
      //nothing for u and r
    }
    else if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT)
    {
      for(size_t i =0; i < n ; i++)
      {
        v[i] = v[i]*NM_triplet(data->B_for_M->D2)->x[i];
      }
      for(size_t i =0; i < m ; i++)
      {
        r[i] = r[i]*NM_triplet(data->B_for_H->D2)->x[i];
        u[i] = u[i]/NM_triplet(data->B_for_H->D2)->x[i];
      }
    }
    else
      numerics_printf_verbose(1,"---- GFC3D - BALANCING - rescaling type is not implemented");

  }

  //else continue;
}


GlobalFrictionContactProblem* gfc3d_balancing_free(GlobalFrictionContactProblem* balanced_problem,
                                                   SolverOptions* options)
{
  assert(balanced_problem);
  GlobalFrictionContactProblem_balancing_data  *balancing_data = (GlobalFrictionContactProblem_balancing_data * ) balanced_problem->env;
  if (balancing_data)
  {
    balanced_problem->env = gfc3d_balancing_data_free(balancing_data);
    globalFrictionContact_free(balanced_problem);
    return NULL;
  }
  else
    return balanced_problem;
}

GlobalFrictionContactProblem_balancing_data   * gfc3d_balancing_data_new()
{
  GlobalFrictionContactProblem_balancing_data  * data = malloc(sizeof(GlobalFrictionContactProblem_balancing_data));
  data->B_for_M =NULL;
  data->B_for_H =NULL;
  return data;
}

GlobalFrictionContactProblem_balancing_data  * gfc3d_balancing_data_free
(GlobalFrictionContactProblem_balancing_data * data)
{
  if (data->B_for_M)
    data->B_for_M = NM_BalancingMatrices_free(data->B_for_M);
  if (data->B_for_H)
    data->B_for_H = NM_BalancingMatrices_free(data->B_for_H);
  free(data);
  return NULL;
}
