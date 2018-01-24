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


#include "ConvexQP_Solvers.h"
#include "ConvexQP_computeError.h"
#include "NumericsMatrix.h"
#include "SparseBlockMatrix.h"
#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include "NumericsSparseMatrix.h"
#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ConvexQP_cst.h"
#include "numerics_verbose.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

const char* const   SICONOS_CONVEXQP_ADMM_STR = "CONVEXQP ADMM";
typedef struct {
  double * xi_hat;
  double * u_hat;
  double * u_k;
  double * xi_k;
}
  ConvexQP_ADDM_data;

void convexQP_ADMM_init(ConvexQP* problem, SolverOptions* options)
{
  int n = problem->size;
  int m = problem->m;
  if (!options->dWork || options->dWorkSize != 2*m+n)
  {
    options->dWork = (double*)calloc(2*m+n,sizeof(double));
    options->dWorkSize = 2*m+n;
  }
  if  (options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] == SICONOS_CONVEXQP_ADMM_ACCELERATION ||
       options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] == SICONOS_CONVEXQP_ADMM_ACCELERATION_AND_RESTART )
  {
    options->solverData=(ConvexQP_ADDM_data *)malloc(sizeof(ConvexQP_ADDM_data));
    ConvexQP_ADDM_data * data = (ConvexQP_ADDM_data *)options->solverData;
    data->xi_hat =  (double*)calloc(m,sizeof(double));
    data->u_hat =  (double*)calloc(m,sizeof(double));
    data->xi_k =  (double*)calloc(m,sizeof(double));
    data->u_k =  (double*)calloc(m,sizeof(double));
  }
}
void convexQP_ADMM_free(ConvexQP* problem, SolverOptions* options)
{
  if (options->dWork)
  {
    free(options->dWork);
    options->dWork=NULL;
    options->dWorkSize = 0;
  }
  if (options->solverData)
  {
    ConvexQP_ADDM_data * data = (ConvexQP_ADDM_data *)options->solverData;
    free(data->xi_hat);
    free(data->u_hat);
    free(data->xi_k);
    free(data->u_k);
    free(data);
  }

}

void convexQP_ADMM(ConvexQP* problem,
                   double *z, double * w,
                   double * xi, double *u,
                   int* info, SolverOptions* options)
{

  //DEBUG_EXPR(convexQP_display(problem););
/* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  if (verbose > 0)
  {
    solver_options_print(options);
  }


  double* q = problem->q;
  NumericsMatrix* M = problem->M;

  NumericsMatrix* A = problem->A;
  double * b = problem->b;

  /* Dimension of the problem */
  int n =  problem->size;
  int AisIdentity = 0;
  if (!problem->istheNormConvexQPset)
  {
    problem->normConvexQP= cblas_dnrm2(n , problem->q , 1);
    DEBUG_PRINTF("problem->norm ConvexQP= %12.8e\n", problem->normConvexQP);
    problem->istheNormConvexQPset=1;
  }
  double norm_q =problem->normConvexQP;
  DEBUG_PRINTF("norm_q = %12.8e\n", norm_q);

  if (!A) /* A is considered to be the identity and b =0 */
  {
    DEBUG_PRINT("A is not given and then considered as Identity and b =0\n");
    problem->m = n;
    problem->A= NM_create(NM_SPARSE,n,n);
    NM_triplet_alloc(problem->A,0);
    problem->A->matrix2->origin= NSM_TRIPLET;

    for (int k =0; k< n; k++)
    {
      NM_zentry(problem->A, k, k, 1);
    }
    DEBUG_EXPR(NM_display(problem->A));

    if (!b)
    {
      problem->b = (double * ) calloc(n,sizeof(double));
    }
    A = problem->A;
    b = problem->b;
    AisIdentity = 1;
  }


  int m =  problem->m;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  /*****  ADMM iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double rho = 0.0;
  rho = dparam[SICONOS_CONVEXQP_ADMM_RHO];
  if (rho == 0.0)
    numerics_error("ConvexQP_ADMM", "dparam[SICONOS_CONVEXQP_PGOC_RHO] must be nonzero");

  /* double tau=1; */
  int internal_allocation=0;
  if (!options->dWork || options->dWorkSize != 2*m+n)
  {
    convexQP_ADMM_init(problem, options);
    internal_allocation = 1;
  }

  double * tmp =  options->dWork;

  int accelerated = iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION];

  /* z_k = (double *)malloc(n * sizeof(double)); */
  /* w_k = (double *)malloc(n * sizeof(double)); */
  /* u_k = (double *)calloc(m,sizeof(double)); */
  /* xi_k = (double *)calloc(m,sizeof(double)); */
  /* xi_tmp = (double *)calloc(m,sizeof(double)); */


  /* Compute M + rho A^T A (storage in M)*/
  NumericsMatrix *Atrans;
  if (!A)
  {
    if (M->storageType != A->storageType)
    {
      numerics_error("ConvexQP_ADMM", "M and A must have the same storage");
    }
  }


  NumericsMatrix *W = NM_new();

  NM_copy(M, W);

  if (AisIdentity)
  {
    NM_add_to_diag3(W, rho);
  }
  else
  {
    Atrans = NM_transpose(A);
    NM_gemm(rho, Atrans, A, 1.0, W);
  }


  if (accelerated == SICONOS_CONVEXQP_ADMM_NO_ACCELERATION)
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      /********************/
      /*  1 - Compute z */
      /********************/

      /* compute the rhs */
      /* q --> z */
      cblas_dcopy(n , q , 1 , z, 1);
      cblas_dscal(n, -1, z,1);

      DEBUG_EXPR(NV_display(u,m));
      DEBUG_EXPR(NV_display(xi,m));

      /*  u -b + xi_k --> u */
      cblas_daxpy(m, -1.0, b, 1, u , 1);
      DEBUG_EXPR(NV_display(u,m));
      cblas_daxpy(m, 1.0, xi, 1, u , 1);
      DEBUG_EXPR(NV_display(u,m));

      if (AisIdentity)
      {
        cblas_daxpy(m, rho, u, 1, z, 1);
      }
      else
      {
        NM_gemv(rho, Atrans, u, 1.0, z);
      }
      DEBUG_PRINT("rhs:")
        DEBUG_EXPR(NV_display(z,n));

      /* Linear system solver */
      /* cblas_dcopy(n , w_k , 1 , z, 1); */
      NM_gesv_expert(W,z,NM_KEEP_FACTORS);
      DEBUG_PRINT("z:")
        DEBUG_EXPR(NV_display(z,n));

      /********************/
      /*  2 - Compute u */
      /********************/

      /* A z_k - xi_k + b */
      cblas_dcopy(m , b , 1 , tmp, 1);
      cblas_daxpy(m, -1, xi, 1, tmp , 1);
      if (AisIdentity)
      {
        cblas_daxpy(m, 1.0, z, 1, tmp , 1);
      }
      else
      {
        NM_gemv(1.0, A, z, 1.0, tmp);
      }
      DEBUG_PRINT("before projection")
        DEBUG_EXPR(NV_display(tmp,m));

      problem->ProjectionOnC(problem,tmp,u);

      DEBUG_EXPR(NV_display(u,m));

      /**********************/
      /*  3 - Compute xi */
      /**********************/

      /* xi_k + A z_k -u_k +b -> xi_k */
      cblas_daxpy(m, -1, b, 1, xi , 1);

      if (AisIdentity)
      {
        cblas_daxpy(m, -1.0, z, 1, xi , 1);
      }
      else
      {
        NM_gemv(-1.0, A, z, 1.0, xi);
      }

      cblas_daxpy(m, 1, u, 1, xi , 1);
      DEBUG_EXPR(NV_display(xi,m));

      DEBUG_EXPR(NV_display(u,m));

      /* **** Criterium convergence **** */
      convexQP_compute_error(problem, z , xi, w, u, tolerance, rho, options, norm_q, &error);

      numerics_printf_verbose(1,"---- ConvexQP - ADMM  - Iteration %i rho = %14.7e \tError = %14.7e", iter, rho, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }

    cblas_dscal(m, rho, xi, 1);
  }
  else if (accelerated == SICONOS_CONVEXQP_ADMM_ACCELERATION ||
           accelerated == SICONOS_CONVEXQP_ADMM_ACCELERATION_AND_RESTART)
  {
    double eta = dparam[SICONOS_CONVEXQP_ADMM_RESTART_ETA];

    ConvexQP_ADDM_data * data = (ConvexQP_ADDM_data *)options->solverData;

    double * u_k = data->u_k;
    double * xi_k =  data->xi_k;
    double * u_hat =  data->u_hat;
    double * xi_hat = data->xi_hat;

    cblas_dcopy(m , xi , 1 , xi_k, 1);
    cblas_dcopy(m , u , 1 , u_k, 1);
    cblas_dcopy(m , xi , 1 , xi_hat, 1);
    cblas_dcopy(m , u , 1 , u_hat, 1);

    double e_k = INFINITY, e, d, alpha;
    double tau , tau_k = 1.0;
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      /********************/
      /*  1 - Compute z */
      /********************/

      /* compute the rhs */
      /* q --> z */
      cblas_dcopy(n , q , 1 , z, 1);
      cblas_dscal(n, -1, z,1);
      /*  u -b + xi_k --> u */
      cblas_dcopy(m , u_hat , 1 , tmp, 1);
      cblas_daxpy(m, -1.0, b, 1, tmp , 1);
      cblas_daxpy(m, 1.0, xi_hat, 1, tmp , 1);
      if (AisIdentity)
      {
        cblas_daxpy(m, rho, tmp, 1, z, 1);
      }
      else
      {
        NM_gemv(rho, Atrans, tmp, 1.0, z);
      }
      DEBUG_PRINT("rhs:");
      DEBUG_EXPR(NV_display(z,n));

      /* Linear system solver */
      /* cblas_dcopy(n , w_k , 1 , z, 1); */
      NM_gesv_expert(W,z,NM_KEEP_FACTORS);
      DEBUG_PRINT("z:");
      DEBUG_EXPR(NV_display(z,n));

      /********************/
      /*  2 - Compute u */
      /********************/

      /* A z_k - xi_k + b */
      cblas_dcopy(m , b , 1 , tmp, 1);
      cblas_daxpy(m, -1, xi_hat, 1, tmp , 1);
      if (AisIdentity)
      {
        cblas_daxpy(m, 1.0, z, 1, tmp , 1);
      }
      else
      {
        NM_gemv(1.0, A, z, 1.0, tmp);
      }
      DEBUG_PRINT("before projection");
      DEBUG_EXPR(NV_display(tmp,m));

      problem->ProjectionOnC(problem,tmp,u);

      DEBUG_EXPR(NV_display(u,m));

      /**********************/
      /*  3 - Compute xi */
      /**********************/

      /* xi_k + A z_k -u_k +b -> xi_k */
      cblas_dcopy(m , xi_hat , 1 , xi, 1);
      cblas_daxpy(m, -1, b, 1, xi , 1);

      if (AisIdentity)
      {
        cblas_daxpy(m, -1.0, z, 1, xi , 1);
      }
      else
      {
        NM_gemv(-1.0, A, z, 1.0, xi);
      }

      cblas_daxpy(m, 1, u, 1, xi , 1);
      DEBUG_EXPR(NV_display(xi,m));



      /**********************/
      /*  3 - Acceleration  */
      /**********************/
      DEBUG_EXPR(NV_display(u_hat,m));
      cblas_dcopy(m , u_hat , 1 , tmp, 1);
      cblas_daxpy(m, -1, u, 1, tmp , 1);
      d = cblas_dnrm2(m , tmp , 1);
      e = d * d;
      DEBUG_EXPR(NV_display(xi_hat,m));
      cblas_dcopy(m , xi_hat , 1 , tmp, 1);
      cblas_daxpy(m, -1, xi, 1, tmp , 1);
      d = cblas_dnrm2(m , tmp , 1);
      e += d * d;

      DEBUG_PRINTF("residual e = %e \n", e);
      DEBUG_PRINTF("residual e_k = %e \n", e_k);
      DEBUG_PRINTF("eta  = %e \n", eta);

      if ((e <  eta * e_k) || (accelerated == SICONOS_CONVEXQP_ADMM_ACCELERATION))
      {
        tau  = 0.5 *(1 +sqrt(1.0+4.0*tau_k*tau_k));
        alpha = (tau_k-1.0)/tau;

        cblas_dcopy(m , u , 1 , u_hat, 1);
        cblas_dscal(m, 1+alpha, u_hat,1);
        cblas_daxpy(m, -alpha, u_k, 1, u_hat , 1);
        DEBUG_EXPR(NV_display(u_hat,m));

        cblas_dcopy(m , xi , 1 , xi_hat, 1);
        cblas_dscal(m, 1+alpha, xi_hat,1);
        cblas_daxpy(m, -alpha, xi_k, 1, xi_hat , 1);
        DEBUG_EXPR(NV_display(xi_hat,m));
        DEBUG_PRINTF("tau  = %e, \t tau_k  = %e \t alpha  = %e   \n", tau, tau_k, alpha);
        tau_k=tau;
        e_k=e;
      }
      else
      {
        tau_k=1.0;
        e_k = e_k /eta;
        DEBUG_PRINTF("tau_k  = %e \t alpha  = %e   \n", tau_k);
        cblas_dcopy(m , xi_k , 1 , xi_hat, 1);
        cblas_dcopy(m , u_k , 1 , u_hat, 1);
      }

      cblas_dcopy(m , xi , 1 , xi_k, 1);
      cblas_dcopy(m , u , 1 , u_k, 1);


      error = sqrt(e);
      if (fabs(norm_q) > DBL_EPSILON)
        error /= norm_q;

      numerics_printf_verbose(1,"---- ConvexQP - ADMM  - Iteration %i rho = %14.7e \t residual = %14.7e", iter, rho, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }

    /* recompute a last value of z with u and xi */
    /* compute the rhs */
    /* q --> z */
    cblas_dcopy(n , q , 1 , z, 1);
    cblas_dscal(n, -1, z,1);
    /*  u -b + xi_k --> u */
    cblas_dcopy(m , u , 1 , tmp, 1);
    cblas_daxpy(m, -1.0, b, 1, tmp , 1);
    cblas_daxpy(m, 1.0, xi, 1, tmp , 1);
    if (AisIdentity)
    {
      cblas_daxpy(m, rho, tmp, 1, z, 1);
    }
    else
    {
      NM_gemv(rho, Atrans, tmp, 1.0, z);
    }
    DEBUG_PRINT("rhs:");
    DEBUG_EXPR(NV_display(z,n));

    /* Linear system solver */
    /* cblas_dcopy(n , w_k , 1 , z, 1); */
    NM_gesv_expert(W,z,NM_KEEP_FACTORS);
    DEBUG_PRINT("z:");
    DEBUG_EXPR(NV_display(z,n));
    
    /* check the full criterion */
    /* **** Criterium convergence **** */
    convexQP_compute_error(problem, z , xi, w, u, tolerance, rho, options, norm_q, &error);
    numerics_printf_verbose(1,"---- ConvexQP - ADMM  - Iteration %i rho = %14.7e \t error = %14.7e", iter, rho, error);

    if (error < tolerance) hasNotConverged = 0;
    else hasNotConverged = 1;
    *info = hasNotConverged;



    cblas_dscal(m, rho, xi, 1);
  }





  //verbose=1;

  numerics_printf_verbose(1,"---- ConvexQP - ADMM - #Iteration %i Final Residual = %14.7e", iter, error);
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /* free(z_k); */
  /* free(w_k); */
  if (internal_allocation)
  {
    convexQP_ADMM_free(problem,options);
  }
  /* free(xi_k); */

}


int convexQP_ADMM_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ADMM Solver\n");
  }

  options->solverId = SICONOS_CONVEXQP_ADMM;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;

  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
  options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] = SICONOS_CONVEXQP_ADMM_NO_ACCELERATION; /* 0 Acceleration */


  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;

  options->dparam[SICONOS_CONVEXQP_ADMM_RHO] = 1.0;
  options->dparam[SICONOS_CONVEXQP_ADMM_RESTART_ETA] = 0.999;

  options->internalSolvers = NULL;


  return 0;
}
