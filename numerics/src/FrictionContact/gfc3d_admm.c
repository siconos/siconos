/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "fc3d_projection.h"
//#include "gfc3d_projection.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
#include "projectionOnCone.h"
#include "SiconosLapack.h"
#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "sanitizer.h"
#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "float.h"
/* #define DEBUG_STDOUT */
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
#include "debug.h"
const char* const   SICONOS_GLOBAL_FRICTION_3D_ADMM_STR = "GFC3D ADMM";

typedef struct {
  double * reaction_hat;
  double * u_hat;
  double * u_k;
  double * reaction_k;
  double * b;
}
  Gfc3d_ADDM_data;




void gfc3d_ADMM_init(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;
  if (!options->dWork || options->dWorkSize != 2*m+n)
  {
    options->dWork = (double*)calloc(2*m+n,sizeof(double));
    options->dWorkSize = 2*m+n;
  }
  if  (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION ||
       options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART )
  {
    options->solverData=(Gfc3d_ADDM_data *)malloc(sizeof(Gfc3d_ADDM_data));
    Gfc3d_ADDM_data * data = (Gfc3d_ADDM_data *)options->solverData;
    data->reaction_hat =  (double*)calloc(m,sizeof(double));
    data->u_hat =  (double*)calloc(m,sizeof(double));
    data->reaction_k =  (double*)calloc(m,sizeof(double));
    data->u_k =  (double*)calloc(m,sizeof(double));
    data->b =  (double*)calloc(m,sizeof(double));
  }
}
void gfc3d_ADMM_free(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  if (options->dWork)
  {
    free(options->dWork);
    options->dWork=NULL;
    options->dWorkSize = 0;
  }
  if (options->solverData)
  {
    Gfc3d_ADDM_data * data = (Gfc3d_ADDM_data *)options->solverData;
    free(data->reaction_hat);
    free(data->u_hat);
    free(data->reaction_k);
    free(data->u_k);
    free(data->b);
    free(data);
  }

}

void gfc3d_ADMM(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
                double* restrict velocity, double* restrict globalVelocity,
                int* restrict info, SolverOptions* restrict options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;
  NumericsMatrix* M = problem->M;
  NumericsMatrix* H = problem->H;
  double* q = problem->q;
  double* mu = problem->mu;

  assert((int)H->size1 == problem->numberOfContacts * problem->dimension);
  assert((int)M->size0 == M->size1);
  assert((int)M->size0 == H->size0); /* size(velocity) ==
                                      * Htrans*globalVelocity */



  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  *info = gfc3d_checkTrivialCaseGlobal(n, q, velocity, reaction, globalVelocity, options);

  if (*info == 0)
    return;

  int contact; /* Number of the current row of blocks in M */

  double norm_q = cblas_dnrm2(n , problem->q , 1);

  double norm_b = cblas_dnrm2(m , problem->b , 1);


  
  numerics_printf_verbose(1,"---- GFC3D - ADMM - Problem information");
  numerics_printf_verbose(1,"---- GFC3D - ADMM - 1-norm of M = %g norm of q = %g ", NM_norm_1(problem->M), norm_q);
  numerics_printf_verbose(1,"---- GFC3D - ADMM - inf-norm of M = %g ", NM_norm_inf(problem->M));

  numerics_printf_verbose(1,"---- GFC3D - ADMM - 1-norm of H = %g norm of b = %g ", NM_norm_1(problem->H), norm_b);
  numerics_printf_verbose(1,"---- GFC3D - ADMM - inf-norm of H = %g ", NM_norm_inf(problem->H));
  numerics_printf_verbose(1,"---- GFC3D - ADMM -  M is symmetric = %i ", NM_is_symmetric(problem->M));


  int internal_allocation=0;
  if (!options->dWork || options->dWorkSize != 2*m+n)
  {
    gfc3d_ADMM_init(problem, options);
    internal_allocation = 1;
  }
  /*****  ADMM iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double rho = 0.0;

  if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
      SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT)
  {
    rho = dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }
  else if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
           SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_NORM_INF)
  {
    double norm_1_M =   NM_norm_1(problem->M);
    double norm_1_H =   NM_norm_1(problem->H);
    if ((fabs(norm_1_H) > DBL_EPSILON) &&  (fabs(norm_1_M) > DBL_EPSILON))
      rho = norm_1_M/norm_1_H;
    else
      rho = dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }
  else if  (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
            SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_ADAPTIVE)
  {
    numerics_error("gfc3d_ADMM", "Adaptive rho stratgey not yet implemented.");
  }
  
  
  if (rho <= DBL_EPSILON)
    numerics_error("gfc3d_ADMM", "dparam[SICONOS_FRICTION_3D_ADMM_RHO] must be nonzero");

  double * tmp =  options->dWork;

  /* Compute M + rho H H^T (storage in M)*/
  NumericsMatrix *Htrans;

  NumericsMatrix *W = NM_new();

  NM_copy(M, W);
  Htrans = NM_transpose(H);
  NM_gemm(rho, H, Htrans, 1.0, W);

  double eta = dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA];

  Gfc3d_ADDM_data * data = (Gfc3d_ADDM_data *)options->solverData;

  double * u = velocity;
  double * z = globalVelocity;

  double * u_k = data->u_k;
  double * reaction_k =  data->reaction_k;
  double * u_hat =  data->u_hat;
  double * reaction_hat = data->reaction_hat;

  double * b = data->b;

  cblas_dscal(m, 1.0/rho, reaction, 1);

  cblas_dcopy(m , reaction , 1 , reaction_k, 1);
  cblas_dcopy(m , u , 1 , u_k, 1);

  cblas_dcopy(m , reaction , 1 , reaction_hat, 1);
  cblas_dcopy(m , u , 1 , u_hat, 1);

  double e_k = INFINITY, e, d, alpha, r, s, residual;
  double tau , tau_k = 1.0;
  int pos;
  double normUT;
  while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      /********************/
      /*  0 - Compute b   */
      /********************/

      cblas_dcopy(m,problem->b,1,b,1);
      for (contact = 0 ; contact < nc ; ++contact)
      {
        pos = contact * 3;
        normUT = sqrt(u[pos + 1] * u[pos + 1] + u[pos + 2] * u[pos + 2]);
        b[pos] +=  problem->mu[contact]*normUT;

      }


      /********************/
      /*  1 - Compute z */
      /********************/

      /* compute the rhs */
      /* q --> z */
      cblas_dcopy(n , q , 1 , z, 1);
      //cblas_dscal(n, -1, z,1);

      /*  u -b + reaction_k --> u */
      cblas_dcopy(m , u_hat , 1 , tmp, 1);
      cblas_daxpy(m, -1.0, b, 1, tmp , 1);
      cblas_daxpy(m, 1.0, reaction_hat, 1, tmp , 1);
      NM_gemv(rho, H, tmp, 1.0, z);


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

      /* H^ z_k - reaction_k + b */
      cblas_dcopy(m , b , 1 , u, 1);
      cblas_daxpy(m, -1, reaction_hat, 1, u , 1);
      NM_gemv(1.0, Htrans, z, 1.0, u);

      DEBUG_PRINT("before projection");
      DEBUG_EXPR(NV_display(u,m));

      /* Loop through the contact points */
      for (contact = 0 ; contact < nc ; ++contact)
      {
        pos = contact * 3;
        projectionOnCone(&u[pos], 1.0/mu[contact]);
      }


      DEBUG_EXPR(NV_display(u,m));

      /**********************/
      /*  3 - Compute reaction */
      /**********************/


      /* - H^T z_k + u_k -b ->  reaction (residual) */
      cblas_dcopy(m , u, 1 , reaction, 1);
      cblas_daxpy(m, -1, b, 1, reaction , 1);
      NM_gemv(-1.0, Htrans, z, 1.0, reaction);
      d = cblas_dnrm2(m , reaction , 1);
      r = d * d;

      /* reaction_hat -  A z_k + u_k -b ->  xi */
      cblas_daxpy(m, 1, reaction_hat, 1, reaction , 1);

      /**********************/
      /*  3 - Acceleration  */
      /**********************/

      DEBUG_EXPR(NV_display(u_hat,m));
      cblas_dcopy(m , u_hat , 1 , tmp, 1);
      cblas_daxpy(m, -1, u, 1, tmp , 1);
      d = cblas_dnrm2(m , tmp , 1);
      s = d * d;

      e =r+s;

      DEBUG_PRINTF("residual e = %e \n", e);
      DEBUG_PRINTF("residual r = %e \n", r);
      DEBUG_PRINTF("residual s = %e \n", s);
      DEBUG_PRINTF("residual e_k = %e \n", e_k);
      DEBUG_PRINTF("eta  = %e \n", eta);

      if ((e <  eta * e_k))
      {
        tau  = 0.5 *(1 +sqrt(1.0+4.0*tau_k*tau_k));
        alpha = (tau_k-1.0)/tau;

        cblas_dcopy(m , u , 1 , u_hat, 1);
        cblas_dscal(m, 1+alpha, u_hat,1);
        cblas_daxpy(m, -alpha, u_k, 1, u_hat , 1);
        DEBUG_EXPR(NV_display(u_hat,m));

        cblas_dcopy(m , reaction , 1 , reaction_hat, 1);
        cblas_dscal(m, 1+alpha, reaction_hat,1);
        cblas_daxpy(m, -alpha, reaction_k, 1, reaction_hat , 1);
        DEBUG_EXPR(NV_display(reaction_hat,m));
        DEBUG_PRINTF("tau  = %e, \t tau_k  = %e \t alpha  = %e   \n", tau, tau_k, alpha);
        tau_k=tau;
        e_k=e;
      }
      else
      {
        tau_k=1.0;
        e_k = e_k /eta;
        DEBUG_PRINTF("tau_k  = %e \t alpha  = %e   \n", tau_k);
        cblas_dcopy(m , reaction_k , 1 , reaction_hat, 1);
        cblas_dcopy(m , u_k , 1 , u_hat, 1);
      }

      cblas_dcopy(m , reaction , 1 , reaction_k, 1);
      cblas_dcopy(m , u , 1 , u_k, 1);

      residual = sqrt(e);
      if (fabs(norm_q) > DBL_EPSILON)
        residual /= norm_q;

      numerics_printf_verbose(1,"---- GFC3D - ADMM  - Iteration %i rho = %14.7e, residual = %14.7e, tol = %14.7e", iter, rho, residual, tolerance);

      if (residual < tolerance)
      {
        /* check the full criterion */
        cblas_dscal(m, rho, reaction, 1);
        gfc3d_compute_error(problem,  reaction, u, z,  tolerance, options, norm_q, &error);
        if (error < dparam[SICONOS_DPARAM_TOL])
        {
          hasNotConverged = 0;
          numerics_printf_verbose(1,"---- GFC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
        }
        else
        {
          numerics_printf_verbose(1,"---- GFC3D - ADMM  - The tolerance on the  residual is not sufficient to reach accuracy (error =  %14.7e)", error);
          tolerance = tolerance * residual/error;
          numerics_printf_verbose(1,"---- GFC3D - ADMM  - We reduce the tolerance on the residual to %14.7e", tolerance);
          cblas_dscal(m, 1.0/rho, reaction, 1);
        }
      }
      *info = hasNotConverged;
    }

  if (iter==itermax)
  {
    cblas_dscal(m, rho, reaction, 1);
    gfc3d_compute_error(problem,  reaction, u, z,  tolerance, options, norm_q, &error);
    numerics_printf_verbose(1,"---- GFC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
  }
  
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;



  /***** Free memory *****/
  if (internal_allocation)
  {
    gfc3d_ADMM_free(problem,options);
  }
}



int gfc3d_ADMM_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ADMM Solver\n");
  }

  options->solverId = SICONOS_GLOBAL_FRICTION_3D_ADMM;

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
  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] = SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;


  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
    SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_NORM_INF;
  
  options->dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 1.0;
  options->dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA] = 0.999;

  options->internalSolvers = NULL;


  return 0;
}
