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
//#include "fc3d_projection.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
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
const char* const   SICONOS_FRICTION_3D_ADMM_STR = "FC3D ADMM";

typedef struct
{
  double * xi;
  double * xi_k;
  double * xi_hat;
  double * z;
  double * z_k;
  double * z_hat;

  double * q;
}
Fc3d_ADDM_data;




void fc3d_admm_init(FrictionContactProblem* problem, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  /* int n = problem->M->size0; */
  int m = 3 * nc;
  /* if(!options->dWork || options->dWorkSize != m) */
  /* { */
  /*   options->dWork = (double*)calloc(m,sizeof(double)); */
  /*   options->dWorkSize = m; */
  /* } */
  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION ||
      options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART)
  {
    options->solverData=(Fc3d_ADDM_data *)malloc(sizeof(Fc3d_ADDM_data));
    Fc3d_ADDM_data * data = (Fc3d_ADDM_data *)options->solverData;
    data->xi_hat = (double*)calloc(m,sizeof(double));
    data->xi = (double*)calloc(m,sizeof(double));
    data->xi_k = (double*)calloc(m,sizeof(double));

    data->z_hat = (double*)calloc(m,sizeof(double));
    data->z = (double*)calloc(m,sizeof(double));
    data->z_k = (double*)calloc(m,sizeof(double));
    data->q = (double*)calloc(m,sizeof(double));
  }
}
void fc3d_admm_free(FrictionContactProblem* problem, SolverOptions* options)
{
  if(options->dWork)
  {
    free(options->dWork);
    options->dWork=NULL;
    options->dWorkSize = 0;
  }
  if(options->solverData)
  {
    Fc3d_ADDM_data * data = (Fc3d_ADDM_data *)options->solverData;
    free(data->xi);
    free(data->xi_hat);
    free(data->z_hat);
    free(data->xi_k);
    free(data->z_k);
    free(data->z);
    free(data->q);
    free(data);
  }

}

void fc3d_admm(FrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity,
               int* restrict info, SolverOptions* restrict options)
{

  /* verbose=1; */
  /* frictionContact_display(problem); */
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* int n = problem->M->size0; */
  int m = 3 * nc;
  NumericsMatrix* M = problem->M;
  double* q = problem->q;
  double* mu = problem->mu;

  assert((int)M->size0 == M->size1);

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  /* *info = fc3d_checkTrivialCase(problem, velocity, reaction, options); */

  if(*info == 0)
    return;

  int contact; /* Number of the current row of blocks in M */

  double norm_q = cblas_dnrm2(m , problem->q , 1);


  if(!(NM_is_symmetric(M)))

  {
    numerics_warning("fc3d_admm","---- FC3D - ADMM - M is not symmetric.\n");
  }



  numerics_printf_verbose(1,"---- FC3D - ADMM - Problem information");
  numerics_printf_verbose(1,"---- FC3D - ADMM - 1-norm of M = %g norm of q = %g ", NM_norm_1(problem->M), norm_q);
  numerics_printf_verbose(1,"---- FC3D - ADMM - inf-norm of M = %g ", NM_norm_inf(problem->M));

  int internal_allocation=0;
  if(!(Fc3d_ADDM_data *)options->solverData)
  {
    fc3d_admm_init(problem, options);
    internal_allocation = 1;
  }
  /*****  ADMM iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double rho = 0.0, rho_k=0.0, rho_ratio=0.0;
  int is_rho_variable=0, has_rho_changed = 0 ;

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
     SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT)
  {
    rho = dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }
  else if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
           SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_NORM_INF)
  {
    double norm_1_M =   NM_norm_1(problem->M);
    double norm_1_H =   1.0;
    if ((fabs(norm_1_H) > DBL_EPSILON) &&  (fabs(norm_1_M) > DBL_EPSILON))
      rho = norm_1_M/norm_1_H;
    else
      rho = dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }
  else if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
          SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING)
  {
    rho = dparam[SICONOS_FRICTION_3D_ADMM_RHO];
    is_rho_variable = 1 ;
  }
  rho_k=rho;
  has_rho_changed = 1;

  if(rho <= DBL_EPSILON)
    numerics_error("fc3d_admm", "dparam[SICONOS_FRICTION_3D_ADMM_RHO] (rho) must be nonzero");

  /* Compute M + rho I (storage in W)*/
  NumericsMatrix *W = NM_new();

  /* NM_copy(M, W); */
  /* NM_add_to_diag3(W, rho); */

  double eta = dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA];
  double br_tau = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU];
  double br_phi = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI];

  assert(br_tau > 1);
  assert(br_phi > 1);

  Fc3d_ADDM_data * data = (Fc3d_ADDM_data *)options->solverData;

  /* we use velocity as a tmp */
  double * tmp = velocity;

  double * z = data->z;;
  double * z_k = data->z_k;
  double * z_hat =  data->z_hat;


  double * xi =  data->xi;
  double * xi_k =  data->xi_k;
  double * xi_hat = data->xi_hat;

  double * q_s = data->q;



  cblas_dcopy(m , reaction , 1 , z_k, 1);
  cblas_dcopy(m , reaction , 1 , z_hat, 1);

  cblas_dscal(m, 1.0/rho, velocity, 1);
  cblas_dcopy(m , velocity , 1 , xi, 1);
  cblas_dcopy(m , velocity , 1 , xi_k, 1);
  cblas_dcopy(m , velocity , 1 , xi_hat, 1);

  double e_k = INFINITY, e,  alpha, r, s, residual;
  double tau , tau_k = 1.0;
  int pos;
  double normUT;

  while((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    DEBUG_PRINTF("\n\n\n############### iteration:%i\n", iter);

    if (has_rho_changed)
    {
      NM_copy(M,W);
      NM_add_to_diag3(W, rho);
    }

    /********************/
    /*  0 - Compute q(s)   */
    /********************/

    cblas_dcopy(m,q,1,q_s,1);
    for(contact = 0 ; contact < nc ; ++contact)
    {
      pos = contact * 3;
      normUT = rho*sqrt(xi[pos + 1] * xi[pos + 1] + xi[pos + 2] * xi[pos + 2]);
      q_s[pos] +=  problem->mu[contact]*normUT;
    }

    /********************/
    /*  1 - Compute r */
    /********************/

    /* compute the rhs */
    /* -q_s  --> reaction */
    cblas_dcopy(m , q_s, 1 , reaction, 1);
    cblas_dscal(m, -1.0, reaction,1);

    /* -q_s - rho * (xi_hat - z_hat )--> reaction */
    cblas_dcopy(m , xi_hat , 1 , tmp, 1);
    cblas_daxpy(m, -1.0, z_hat, 1, tmp , 1);
    cblas_daxpy(m, -1.0*rho, tmp, 1, reaction, 1);

    DEBUG_PRINT("rhs:");
    DEBUG_EXPR(NV_display(reaction,m));

    /* Linear system solver */
    NM_gesv_expert(W,reaction, NM_KEEP_FACTORS);
    DEBUG_PRINT("reaction:");
    DEBUG_EXPR(NV_display(reaction,m));

    /********************/
    /*  2 - Compute z */
    /********************/

    /* reaction  + xi_hat  --> z */
    cblas_dcopy(m , xi_hat  , 1 , z, 1);
    cblas_daxpy(m, 1, reaction, 1, z , 1);


    DEBUG_PRINT("Before projection :");
    DEBUG_EXPR(NV_display(z,m));

    /* Loop through the contact points */
    for(contact = 0 ; contact < nc ; ++contact)
    {
      projectionOnCone(&z[contact * 3], mu[contact]);
    }
    DEBUG_PRINT("After projection :");
    DEBUG_EXPR(NV_display(z,m));

    /**********************/
    /*  3 - Compute xi */
    /**********************/

    /* r - z --> residual  */
    cblas_dcopy(m , reaction, 1 , xi, 1);
    cblas_daxpy(m, -1.0, z, 1, xi , 1);
    r = cblas_dnrm2(m , xi , 1);


    cblas_daxpy(m, 1.0, xi_hat, 1, xi , 1);
    DEBUG_PRINT("xi : ")
    DEBUG_EXPR(NV_display(xi,m));

    /**********************/
    /*  3 - Residual  */
    /**********************/

    cblas_dcopy(m , z_hat , 1 , tmp, 1);
    cblas_daxpy(m, -1, z, 1, tmp , 1);

    s = rho*cblas_dnrm2(m , tmp , 1);
    /* s = cblas_dnrm2(m , tmp , 1); */

    e =r*r+s*s;

    DEBUG_PRINTF("residual e = %e \n", e);
    DEBUG_PRINTF("residual r = %e \n", r);
    DEBUG_PRINTF("residual s = %e \n", s);
    DEBUG_PRINTF("residual e_k = %e \n", e_k);
    DEBUG_PRINTF("eta  = %e \n", eta);

    /* printf("residual e = %e \n", e); */
    /* printf("residual r = %e \n", r); */
    /* printf("residual s = %e \n", s); */
    /* printf("residual e_k = %e \n", e_k); */
    
    /*********************************/
    /*  3 - Acceleration and restart */
    /*********************************/

    if((e <  eta * e_k))
    {
      tau  = 0.5 *(1 +sqrt(1.0+4.0*tau_k*tau_k));
      alpha = (tau_k-1.0)/tau;

      cblas_dcopy(m , z , 1 , z_hat, 1);
      cblas_dscal(m, 1+alpha, z_hat,1);
      cblas_daxpy(m, -alpha, z_k, 1, z_hat , 1);
      DEBUG_EXPR(NV_display(z_hat,m));

      cblas_dcopy(m , xi , 1 , xi_hat, 1);
      cblas_dscal(m, 1+alpha, xi_hat,1);
      cblas_daxpy(m, -alpha, xi_k, 1, xi_hat , 1);
      DEBUG_EXPR(NV_display(xi_hat,m));
      DEBUG_PRINTF("Accelerate :tau  = %e, \t tau_k  = %e, \t alpha  = %e   \n", tau, tau_k, alpha);
      /* printf("Accelerate :tau  = %e, \t tau_k  = %e, \t alpha  = %e   \n", tau, tau_k, alpha); */
      tau_k=tau;
      e_k=e;
    }
    else
    {
      tau_k=1.0;
      e_k = e_k /eta;
      DEBUG_PRINTF(" Restart tau_k  = %e  \n", tau_k);
      /* printf(" Restart tau_k  = %e  \n", tau_k); */
      cblas_dcopy(m , xi_k , 1 , xi_hat, 1);
      cblas_dcopy(m , z_k , 1 , z_hat, 1);
    }


    rho_k = rho ;

    if (is_rho_variable && iter >= 1)
    {
      if (r > br_phi * s)
      {
        rho = br_tau* rho_k;
        has_rho_changed = 1;
      }
      else if (s > br_phi * r)
      {
        rho = rho_k/br_tau;
        has_rho_changed = 1;
      }
      else
      {
        /* keep the value of rho */
        has_rho_changed = 0;
      }
    }
    rho_ratio = rho_k/rho;

    /* rho =1.0; */
    /* rho_k=1.0; */
    /* rho_ratio = rho_k/rho; */
    /* has_rho_changed = 0; */
    DEBUG_PRINTF("rho =%e\t,rho_k =%e \n", rho, rho_k);
    /* printf("rho =%e\t,rho_k =%e \n", rho, rho_k); */

    cblas_dscal(m, rho_ratio, xi,1);
    cblas_dscal(m, rho_ratio, xi_hat,1);

    /* Next step */
    cblas_dcopy(m , z , 1 , z_k, 1);
    cblas_dcopy(m , xi , 1 , xi_k, 1);

    residual = sqrt(e);
    if(fabs(norm_q) > DBL_EPSILON)
      residual /= norm_q;

    numerics_printf_verbose(1,"---- FC3D - ADMM  - Iteration %i rho = %14.7e, residual = %14.7e, tol = %14.7e", iter, rho, residual, tolerance);

    if(residual < tolerance)
    {
      fc3d_compute_error(problem,  reaction, velocity, tolerance, options, norm_q, &error);
      DEBUG_EXPR(NV_display(velocity,m));
      if(error < dparam[SICONOS_DPARAM_TOL])
      {
        hasNotConverged = 0;
        numerics_printf_verbose(1,"---- FC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
      }
      else
      {
        numerics_printf_verbose(1,"---- FC3D - ADMM  - The tolerance on the  residual is not sufficient to reach accuracy (error =  %14.7e)", error);
        tolerance = tolerance * residual/error;
        numerics_printf_verbose(1,"---- FC3D - ADMM  - We reduce the tolerance on the residual to %14.7e", tolerance);
      }
    }
    *info = hasNotConverged;
  }

  if(iter==itermax)
  {
    fc3d_compute_error(problem,  reaction, velocity, tolerance, options, norm_q, &error);
    numerics_printf_verbose(1,"---- FC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
  }

  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;



  /***** Free memory *****/
  if(internal_allocation)
  {
    fc3d_admm_free(problem,options);
  }
}



int fc3d_admm_setDefaultSolverOptions(SolverOptions* options)
{
  if(verbose > 0)
  {
    printf("Set the Default SolverOptions for the ADMM Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_ADMM;

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
  options->dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU]=1.5;
  options->dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI]=10.0;

  options->internalSolvers = NULL;
  options->solverData = NULL;

  return 0;
}
