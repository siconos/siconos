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
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "projectionOnCone.h"

#include "SiconosLapack.h"

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "sanitizer.h"


#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "NumericsSparseMatrix.h"


/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
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
  double * b;
}
Fc3d_ADMM_data;




void fc3d_admm_init(FrictionContactProblem* problem, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  /* int n = problem->M->size0; */
  int m = 3 * nc;


  int nb_constraints =m;
  options->solverData=(Fc3d_ADMM_data *)malloc(sizeof(Fc3d_ADMM_data));
  Fc3d_ADMM_data * data = (Fc3d_ADMM_data *)options->solverData;
  if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY||
      options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY )
  {
    nb_constraints =2*m;
  }

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION ||
     options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART||
     options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION)
  {

    data->xi_hat = (double*)calloc(nb_constraints,sizeof(double));
    data->xi = (double*)calloc(nb_constraints,sizeof(double));
    data->xi_k = (double*)calloc(nb_constraints,sizeof(double));

    data->z_hat = (double*)calloc(nb_constraints,sizeof(double));
    data->z = (double*)calloc(nb_constraints,sizeof(double));
    data->z_k = (double*)calloc(nb_constraints,sizeof(double));
    data->q = (double*)calloc(m,sizeof(double));
  }

  if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY||
      options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY )
  {
    data->b = (double*)calloc(nb_constraints,sizeof(double));
    if(!options->dWork || options->dWorkSize != 2*m)
    {
      options->dWork = (double*)calloc(2*m,sizeof(double));
      options->dWorkSize = 2*m;
  }
  }
  else
  {
    data->b =NULL;
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
    Fc3d_ADMM_data * data = (Fc3d_ADMM_data *)options->solverData;
    free(data->xi);
    free(data->xi_hat);
    free(data->z_hat);
    free(data->xi_k);
    free(data->z_k);
    free(data->z);
    free(data->q);
    free(data->b);
    free(data);
  }

}


static void fc3d_admm_symmetric(FrictionContactProblem* restrict problem,
                                double* restrict reaction,
                                double* restrict velocity,
                                int* restrict info, SolverOptions* restrict options,
                                double rho,  int is_rho_variable,
                                double norm_q)
{

  /* verbose=2;  */
  /* frictionContact_display(problem); */
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* int n = problem->M->size0; */
  int m = 3 * nc;
  NumericsMatrix* M=  NULL;

  /* if SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE = SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE,
     we force the copy into a NM_SPARSE storageType */
  if(iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE
     && problem->M->storageType == NM_SPARSE_BLOCK)
  {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    M = NM_create(NM_SPARSE,  problem->M->size0,  problem->M->size1);
    NM_copy_to_sparse(problem->M, M);
    NSM_diag_indices(M);
  }
  else
  {
    M = problem->M;
  }
  double* q = problem->q;
  double* mu = problem->mu;
  double alpha_r=0.0;
  FrictionContactProblem *  rescaled_problem =  problem;
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_YES)
  {
    alpha_r = NM_norm_inf(M);
    //alpha_r=1./10.0;

    rescaled_problem =  frictionContact_copy(problem);
    frictionContact_rescaling(rescaled_problem, 1./alpha_r, 1.0);

    M = rescaled_problem->M;
    q = rescaled_problem->q;
    norm_q = cblas_dnrm2(m , problem->q , 1);
    printf("norm_q = %e\n", norm_q);
    norm_q = cblas_dnrm2(m , rescaled_problem->q , 1);
    printf("norm_q (rescaled) = %e\n", norm_q);
  }

  /* Compute M + rho I (storage in W)*/
  NumericsMatrix *W = NM_new();


  /*****  ADMM iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[0];

  double eta = dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA];
  double br_tau = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU];
  double br_phi = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI];

  assert(br_tau > 1);
  assert(br_phi > 1);

  Fc3d_ADMM_data * data = (Fc3d_ADMM_data *)options->solverData;


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

  int contact; /* Number of the current row of blocks in M */

  double rho_k=0.0, rho_ratio=0.0;
  double e_k = INFINITY, e,  alpha, r, s, residual, r_scaled, s_scaled;
  double norm_r=0.0, norm_z=0.0;
  double tau , tau_k = 1.0;
  int pos;
  double normUT;
  int admm_has_converged=0;

  rho_k=rho;
  int has_rho_changed = 1;

  while((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    DEBUG_PRINTF("\n\n\n############### iteration:%i\n", iter);

    if (has_rho_changed)
    {
      /* NM_clear(W); */
      /* W= NM_new(); */
      NM_copy(M,W);
      NM_add_to_diag3(W, rho);
    }

    /********************/
    /*  0 - Compute q(s)   */
    /********************/

    cblas_dcopy(m,q,1,q_s,1);
    DEBUG_EXPR(fc3d_compute_error(problem,  reaction, velocity, tolerance, options, norm_q, &error););
    DEBUG_EXPR(NV_display(velocity,m););
    DEBUG_EXPR(NV_display(xi,m););

    cblas_dcopy(m , q, 1 , velocity, 1);
    NM_gemv(1.0, M, reaction, 1.0, velocity);

    for(contact = 0 ; contact < nc ; ++contact)
    {
      pos = contact * 3;
      /* xi is equal to velocity on the tangential at convergence */
      /* normUT = rho*sqrt(xi[pos + 1] * xi[pos + 1] + xi[pos + 2] * xi[pos + 2]); */
      normUT = sqrt(velocity[pos + 1] * velocity[pos + 1] + velocity[pos + 2] *velocity[pos + 2]);
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
    /*  3 - Compute xi    */
    /**********************/

    /* r - z --> residual  */
    cblas_dcopy(m , reaction, 1 , xi, 1);
    cblas_daxpy(m, -1.0, z, 1, xi , 1);
    r = cblas_dnrm2(m , xi , 1);

    norm_r = cblas_dnrm2(m , reaction , 1);
    norm_z = cblas_dnrm2(m , z , 1);


    cblas_daxpy(m, 1.0, xi_hat, 1, xi , 1);
    DEBUG_PRINT("xi : ")
    DEBUG_EXPR(NV_display(xi,m));

    /**********************/
    /*  3 - Residual      */
    /**********************/
    /* s = rho * (z_hat-z) */


    cblas_dcopy(m , z_hat , 1 , tmp, 1);
    cblas_daxpy(m, -1, z, 1, tmp , 1);

    s = rho*cblas_dnrm2(m , tmp , 1);
    double norm_rhoxi = rho*cblas_dnrm2(m , xi , 1);

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
    if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION ||
        options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART)
    {
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
        numerics_printf_verbose(2, "Accelerate :tau  = %e, \t tau_k  = %e, \t alpha  = %e ", tau, tau_k, alpha);
        tau_k=tau;
        e_k=e;
      }
      else
      {
        tau_k=1.0;
        e_k = e_k /eta;
        DEBUG_PRINTF(" Restart tau_k  = %e  \n", tau_k);
        numerics_printf_verbose(2," Restart tau_k  = %e", tau_k);
        cblas_dcopy(m , xi_k , 1 , xi_hat, 1);
        cblas_dcopy(m , z_k , 1 , z_hat, 1);
      }
    }
    else  if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION)
    {
      tau_k=1.0;
      e_k = e_k /eta;
      numerics_printf_verbose(2,"No acceleration and restart tau_k  = %e  \n", tau_k);
      cblas_dcopy(m , xi_k , 1 , xi_hat, 1);
      cblas_dcopy(m , z_k , 1 , z_hat, 1);
    }
    else
    {
      numerics_error("fc3d_admm", " options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] value is not recognize");
    }
    /*********************************/
    /*  4 - Updating rho             */
    /*********************************/

    rho_k = rho ;

    if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
    {
      r_scaled = r / (fmax(norm_z,norm_r));
      s_scaled = s / (norm_rhoxi);
      numerics_printf_verbose(2, "fc3d_admm. scaling : norm_r  = %e, \t norm_z  = %e, \t norm_rhoxi = %e, \t", norm_r,  norm_z, norm_rhoxi);
      numerics_printf_verbose(2, "fc3d_admm. residuals : r  = %e, \t  s = %e", r, s);
      numerics_printf_verbose(2, "fc3d_admm. scaled residuals : r_scaled  = %e, \t  s_scaled = %e", r_scaled, s_scaled);
    }
    else
    {
      r_scaled = r;
      s_scaled = s;
    }

    if (is_rho_variable)
    {
      if (r_scaled > br_phi * s_scaled)
      {
        rho = br_tau* rho_k;
        has_rho_changed = 1;
      }
      else if (s_scaled > br_phi * r_scaled)
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
    else
    {
      has_rho_changed = 0;
    }
    numerics_printf_verbose(2, "fc3d_admm. rho = %5.2e\t, rho_k = %5.2e\t ", rho, rho_k);
    rho_ratio = rho_k/rho;

    DEBUG_PRINTF("rho =%e\t,rho_k =%e \n", rho, rho_k);

    cblas_dscal(m, rho_ratio, xi,1);
    cblas_dscal(m, rho_ratio, xi_hat,1);

    /* Next step */
    cblas_dcopy(m , z , 1 , z_k, 1);
    cblas_dcopy(m , xi , 1 , xi_k, 1);


    /*********************************/
    /*  4 - Stopping criterium       */
    /*********************************/
    int stopping_criterion =0;
    residual = sqrt(e);
    /* if (fabs(norm_q) > DBL_EPSILON) */
    /*   residual /= norm_q; */
    /* if (residual < tolerance) */
    /*   stopping_criterion =1; */
    double scaling_error_primal  = fmax(norm_z,norm_r) +  sqrt(m);
    double epsilon_primal = tolerance *  scaling_error_primal ;
    double scaling_error_dual = norm_rhoxi + sqrt(m);
    double epsilon_dual =  tolerance * scaling_error_dual;
    if (r < epsilon_primal && s < epsilon_dual)
        stopping_criterion =1;


    numerics_printf_verbose(1,"---- FC3D - ADMM  - Iteration %i rho = %14.7e, residual = %14.7e, tol = %14.7e", iter, rho, residual, tolerance);
    numerics_printf_verbose(1,"---- FC3D - ADMM  -                            primal residual = %14.7e, epsilon_primal = %14.7e", r,  epsilon_primal);
    numerics_printf_verbose(1,"---- FC3D - ADMM  -                            dual residual = %14.7e, epsilon_dual = %14.7e", s,  epsilon_dual);

    if (verbose >1)
      frictionContactProblem_compute_statistics(problem,
                                                reaction,
                                                velocity,
                                                tolerance,
                                                0) ;
    if(stopping_criterion)
    {
      /* check the full criterion */
      if (options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_YES)
      {
        norm_q = cblas_dnrm2(m , problem->q , 1);
      }
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
        tolerance = tolerance * fmax(epsilon_dual/scaling_error_dual ,epsilon_primal/scaling_error_primal )/error;
        numerics_printf_verbose(1,"---- FC3D - ADMM  - We reduce the tolerance on the residual to %14.7e", tolerance);
        if (options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_YES)
        {
          norm_q = cblas_dnrm2(m , rescaled_problem->q , 1);
        }
      }
    }
    *info = hasNotConverged;
  }
  if(iter==itermax)
  {
    norm_q = cblas_dnrm2(m , problem->q , 1);
    fc3d_compute_error(problem,  reaction, velocity, tolerance, options, norm_q, &error);
    if (error < dparam[SICONOS_DPARAM_TOL])
    {
      *info = 0;
    }
    numerics_printf_verbose(1,"---- FC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
  }
  NM_clear(W);
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
}

static void fc3d_admm_asymmetric(FrictionContactProblem* restrict problem,
                                double* restrict reaction,
                                double* restrict velocity,
                                int* restrict info, SolverOptions* restrict options,
                                double rho,  int is_rho_variable,
                                double norm_q)
{

  /* verbose=2; */
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

  /*****  ADMM iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[0];

  double eta = dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA];
  double br_tau = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU];
  double br_phi = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI];

  assert(br_tau > 1);
  assert(br_phi > 1);

  Fc3d_ADMM_data * data = (Fc3d_ADMM_data *)options->solverData;

  double * tmp2 = options->dWork;

  /* we use velocity as a tmp */
  double * tmp = velocity;

  double * z = data->z;;
  double * z_k = data->z_k;
  double * z_hat =  data->z_hat;


  double * xi =  data->xi;
  double * xi_k =  data->xi_k;
  double * xi_hat = data->xi_hat;

  double * q_s = data->q;
  double * b_s = data->b;

  /*  constraints matrix */
  NumericsMatrix *A = NM_new();
  NM_copy(M,A);
  CSparseMatrix* A_triplet = NM_triplet(A);

  A->storageType=NM_SPARSE;
  NM_clearDense(A);
  NM_clearSparseBlock(A);
  for (int i=0; i<m; i++)
  {
    CHECK_RETURN(CSparseMatrix_zentry(A_triplet, i+m , i, 1.0));
  }
  A->size0 = 2*A->size0;

  DEBUG_EXPR(NM_display(A););
  /* NumericsMatrix *_A = NM_new(); */
  /* if (M->storageType == NM_DENSE) */
  /* { */
  /*   NM_copy(A,_A); */
  /*   NM_clear(A); */
  /*   NM_to_dense(_A,A); */
  /* } */
  /* DEBUG_EXPR(NM_display(A);); */

  NumericsMatrix *Atrans = NM_transpose(A);

  /*  cost Matrix */
  NumericsMatrix *M_T = NM_transpose(M);
  NumericsMatrix *M_s = NM_add(1.0, M, 1.0, M_T);
  NM_clear(M_T);
  DEBUG_EXPR(NM_display(M););
  DEBUG_EXPR(NM_display(M_s););


  /*iteration matrix */
  NumericsMatrix *W = NM_new();


  /* /\* initialization *\/ */
  /* cblas_dcopy(m , reaction , 1 , z_k, 1); */
  /* cblas_dcopy(m , reaction , 1 , z_hat, 1); */

  /* cblas_dscal(m, 1.0/rho, velocity, 1); */
  /* cblas_dcopy(m , velocity , 1 , xi, 1); */
  /* cblas_dcopy(m , velocity , 1 , xi_k, 1); */
  /* cblas_dcopy(m , velocity , 1 , xi_hat, 1); */

  int contact; /* Number of the current row of blocks in M */

  double rho_k=0.0, rho_ratio=0.0;
  double e_k = INFINITY, e,  alpha, r, s, residual, r_scaled, s_scaled;
  double norm_Ar=0.0, norm_z=0.0, norm_ATxi=0.0, norm_b_s=0.0;
  double tau , tau_k = 1.0;
  int pos;
  double normUT;
  int admm_has_converged=0;



  rho_k=rho;
  int has_rho_changed = 1;



  while((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    DEBUG_PRINTF("\n\n\n############### iteration:%i\n", iter);

    if (has_rho_changed)
    {
      /* NM_clear(W); */
      /* W= NM_new(); */
      NM_copy(M_s,W);
      NM_gemm(rho, Atrans, A, 1.0, W);
      DEBUG_EXPR(NM_display(W));
    }

    /*******************************/
    /*  0 - Compute q(s) and b(s)  */
    /*******************************/

    cblas_dcopy(m,q,1,q_s,1);
    DEBUG_EXPR(fc3d_compute_error(problem,  reaction, velocity, tolerance, options, norm_q, &error););
    DEBUG_EXPR(NV_display(velocity,m););
    DEBUG_EXPR(NV_display(&xi[m],m););

    cblas_dcopy(m , q, 1 , velocity, 1);
    NM_gemv(1.0, M, reaction, 1.0, velocity);

    for(contact = 0 ; contact < nc ; ++contact)
    {
      pos = contact * 3;
      normUT = sqrt(velocity[pos + 1] * velocity[pos+1] + velocity[pos + 2] * velocity[pos + 2]);
      q_s[pos] +=  problem->mu[contact]*normUT;
    }
    for (int i =0; i< m ; i++)
    {
      b_s[i] = q_s[i];
    }
    DEBUG_EXPR(NV_display(q_s,m));
    DEBUG_EXPR(NV_display(b_s,2*m));

    /********************/
    /*  1 - Compute r */
    /********************/
    /* compute the rhs */
    /* -q_s  --> reaction */
    cblas_dcopy(m , q_s, 1 , reaction, 1);
    cblas_dscal(m, -1.0, reaction,1);

    /* -q_s - rho * A^T(xi_hat +b_s - z_hat )--> reaction */

    cblas_dcopy(2*m , xi_hat , 1 , tmp2, 1);
    cblas_daxpy(2*m, 1.0, b_s, 1, tmp2 , 1);
    cblas_daxpy(2*m, -1.0, z_hat, 1, tmp2 , 1);

    NM_gemv(-1.0*rho, Atrans, tmp2, 1.0, reaction);

    DEBUG_PRINT("rhs:");
    DEBUG_EXPR(NV_display(reaction,m));

    /* Linear system solver */
    NM_gesv_expert(W,reaction, NM_KEEP_FACTORS);
    DEBUG_PRINT("reaction:");
    DEBUG_EXPR(NV_display(reaction,m));

    /********************/
    /*  2 - Compute z */
    /********************/

    /* A * reaction  + b_s + xi_hat  --> z */
    cblas_dcopy(2*m , xi_hat  , 1 , z, 1);
    cblas_daxpy(2*m, 1.0, b_s, 1, z , 1);
    NM_gemv(1.0, A, reaction, 1.0, z);

    DEBUG_PRINT("Before projection :");
    DEBUG_EXPR(NV_display(z,2*m));

    /* Loop through the contact points */
    for(contact = 0 ; contact < nc ; ++contact)
    {
      projectionOnDualCone(&z[contact * 3], mu[contact]);
    }
    for(contact = 0 ; contact < nc ; ++contact)
    {
      projectionOnCone(&z[contact * 3+m], mu[contact]);
    }

    DEBUG_PRINT("After projection :");
    DEBUG_EXPR(NV_display(z,2*m));

    /**********************/
    /*  3 - Compute xi    */
    /**********************/



    /* A * r - z + b_s--> residual  */
    if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
    {
      cblas_dscal(2*m , 0.0, xi, 1);
      NM_gemv(1.0, A, reaction, 1.0, xi);
      norm_Ar = cblas_dnrm2(m , xi , 1);

      cblas_daxpy(2*m , 1.0, b_s, 1 , xi, 1);
      cblas_daxpy(2*m, -1.0, z, 1, xi , 1);

      norm_z = cblas_dnrm2(2*m , z , 1);
      norm_b_s = cblas_dnrm2(2*m , b_s, 1);
    }
    else
    {
      cblas_dcopy(2*m , b_s, 1 , xi, 1);
      cblas_daxpy(2*m, -1.0, z, 1, xi , 1);
      NM_gemv(1.0, A, reaction, 1.0, xi);
    }

    r = cblas_dnrm2(2*m , xi , 1);

    cblas_daxpy(2*m, 1.0, xi_hat, 1, xi , 1);
    DEBUG_EXPR(NV_display(xi,2*m));

    /**********************/
    /*  3 - Residual      */
    /**********************/
    /* s = rho  A^T (z_hat-z) */

    cblas_dcopy(2*m , z_hat , 1 , tmp2, 1);
    cblas_daxpy(2*m, -1.0, z, 1, tmp2 , 1);

    cblas_dscal(m, 0.0, tmp, 1);

    NM_gemv(1.0*rho, Atrans, tmp2, 1.0, tmp);

    s = cblas_dnrm2(m , tmp , 1);
    if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
    {
      cblas_dscal(m, 0.0, tmp, 1);
      NM_gemv(1.0*rho, Atrans, xi, 1.0, tmp);
      norm_ATxi = cblas_dnrm2(m , tmp , 1);

    }
    e =r*r+s*s;

    DEBUG_PRINTF("residual e = %e \n", e);
    DEBUG_PRINTF("residual r = %e \n", r);
    DEBUG_PRINTF("residual s = %e \n", s);
    DEBUG_PRINTF("residual e_k = %e \n", e_k);
    DEBUG_PRINTF("eta  = %e \n", eta);

    /*********************************/
    /*  3 - Acceleration and restart */
    /*********************************/
    if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION ||
        options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART)
    {
      if((e <  eta * e_k))
      {
        tau  = 0.5 *(1 +sqrt(1.0+4.0*tau_k*tau_k));
        alpha = (tau_k-1.0)/tau;

        cblas_dcopy(2*m , z , 1 , z_hat, 1);
        cblas_dscal(2*m, 1+alpha, z_hat,1);
        cblas_daxpy(2*m, -alpha, z_k, 1, z_hat , 1);
        DEBUG_EXPR(NV_display(z_hat,2*m));

        cblas_dcopy(2*m , xi , 1 , xi_hat, 1);
        cblas_dscal(2*m, 1+alpha, xi_hat,1);
        cblas_daxpy(2*m, -alpha, xi_k, 1, xi_hat , 1);
        DEBUG_EXPR(NV_display(xi_hat,2*m));
        DEBUG_PRINTF("Accelerate :tau  = %e, \t tau_k  = %e, \t alpha  = %e   \n", tau, tau_k, alpha);
        numerics_printf_verbose(2, "Accelerate :tau  = %e, \t tau_k  = %e, \t alpha  = %e ", tau, tau_k, alpha);
        tau_k=tau;
        e_k=e;
      }
      else
      {
        tau_k=1.0;
        e_k = e_k /eta;
        DEBUG_PRINTF("Restart tau_k  = %e  \n", tau_k);
        numerics_printf_verbose(2,"Restart tau_k  = %e", tau_k);
        cblas_dcopy(2*m , xi_k , 1 , xi_hat, 1);
        cblas_dcopy(2*m , z_k , 1 , z_hat, 1);
      }
    }
    else  if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION)
    {
      tau_k=1.0;
      e_k = e_k /eta;
      numerics_printf_verbose(2,"Restart tau_k  = %e  \n", tau_k);
      cblas_dcopy(2*m , xi_k , 1 , xi_hat, 1);
      cblas_dcopy(2*m , z_k , 1 , z_hat, 1);
    }
    else
    {
      numerics_error("fc3d_admm", " options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] value is not recognize");
    }
    /*********************************/
    /*  4 - Updating rho             */
    /*********************************/

    rho_k = rho ;
    numerics_printf_verbose(2, "residuals : r  = %e, \t  s = %e", r, s);
    if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
    {
      r_scaled = r / (fmax(fmax(norm_b_s,norm_z),norm_Ar));
      s_scaled = s / (rho*norm_ATxi);
      /* s_scaled = s / (norm_ATxi); */
      numerics_printf_verbose(2, "scaling : norm_Ar  = %e, \t  norm_b_s = %e, \t norm_z  = %e, \t norm_ATxi = %e, \t", norm_Ar, norm_b_s, norm_z, norm_ATxi);

      numerics_printf_verbose(2, "scaled residuals : r_scaled  = %e, \t  s_scaled = %e", r_scaled, s_scaled);
    }
    else
    {
      r_scaled = r;
      s_scaled = s;
    }

    if (is_rho_variable)
    {
      if (r_scaled > br_phi * s_scaled)
      {
        rho = br_tau* rho_k;
        has_rho_changed = 1;
      }
      else if (s_scaled > br_phi * r_scaled)
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
    else
    {
      has_rho_changed = 0;
    }
    numerics_printf_verbose(2, "fc3d_admm. rho = %5.2e\t, rho_k = %5.2e\t, r = %5.2e\t,  s = %5.2e\t", rho, rho_k, r, s);

    rho_ratio = rho_k/rho;
    DEBUG_PRINTF("rho =%e\t,rho_k =%e \n", rho, rho_k);

    cblas_dscal(2*m, rho_ratio, xi,1);
    cblas_dscal(2*m, rho_ratio, xi_hat,1);

    /* Next step */
    cblas_dcopy(2*m , z , 1 , z_k, 1);
    cblas_dcopy(2*m , xi , 1 , xi_k, 1);


    /*********************************/
    /*  4 - Stopping criterium       */
    /*********************************/

    residual = sqrt(e);
    if(fabs(norm_q) > DBL_EPSILON)
      residual /= norm_q;

    numerics_printf_verbose(1,"---- FC3D - ADMM  - Iteration %i rho = %14.7e, residual = %14.7e, tol = %14.7e", iter, rho, residual, tolerance);

    /* 4-1 basic stopping criterion */
    admm_has_converged = 1;

    /* 4-2 Relative stopping criterion */
    /* admm_has_converged = (r < tolerance * fmax(r1,r2)) && (s < tolerance * s1 ) ; */
    /* numerics_printf_verbose(2,"---- FC3D - ADMM  - r1 = %14.7e, r2 = %14.7e, s1 = %14.7e , tolerance = %14.7e", r1, r2, s1, tolerance); */
    /* numerics_printf_verbose(2,"---- FC3D - ADMM  - r = %14.7e, tolerance * fmax(r1,r2) = %14.7e, s = %14.7e,  tolerance * s1  = %14.7e ", r, tolerance * fmax(r1,r2), s,  tolerance * s1); */

    if(admm_has_converged)
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

  NM_clear(A);
  NM_clear(M_s);
  NM_clear(W);

}

static double fc3d_admm_select_rho(NumericsMatrix* M, int * is_rho_variable, SolverOptions* restrict options)
{
  double rho=0.0;
  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] ==
     SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN)
  {
    rho = options->dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }
  else if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] ==
           SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF)
  {
    double norm_1_M =   NM_norm_1(M);
    double norm_1_H =   1.0;
    if ((fabs(norm_1_H) > DBL_EPSILON) &&  (fabs(norm_1_M) > DBL_EPSILON))
      rho = norm_1_M/norm_1_H;
    else
      rho =  options->dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }


  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
          SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING ||
          options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
          SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
  {
    *is_rho_variable = 1 ;
  }
  else
    *is_rho_variable = 0 ;
  return rho;
}

void fc3d_admm(FrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity,
               int* restrict info, SolverOptions* restrict options)
{

  verbose=1;
  /* frictionContact_display(problem); */

  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* int n = problem->M->size0; */
  int m = 3 * nc;

  NumericsMatrix* M = problem->M;
  assert((int)M->size0 == M->size1);


  /* Check for trivial case */
  *info = fc3d_checkTrivialCase(problem, velocity, reaction, options);

  if(*info == 0)
    return;


  double norm_q = cblas_dnrm2(m , problem->q , 1);
  if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO] ==
      SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES)
  {
    numerics_printf_verbose(1,"---- FC3D - ADMM - Problem information");
    numerics_printf_verbose(1,"---- FC3D - ADMM - 1-norm of M = %g norm of q = %g ", NM_norm_1(problem->M), norm_q);
    numerics_printf_verbose(1,"---- FC3D - ADMM - inf-norm of M = %g ", NM_norm_inf(problem->M));
    /* getchar(); */
  }
  int internal_allocation=0;
  if(!(Fc3d_ADMM_data *)options->solverData)
  {
    fc3d_admm_init(problem, options);
    internal_allocation = 1;
  }

  int is_rho_variable=0;
  double rho = fc3d_admm_select_rho(M, &is_rho_variable, options);

  if(rho <= DBL_EPSILON)
    numerics_error("fc3d_admm", "dparam[SICONOS_FRICTION_3D_ADMM_RHO] (rho) must be nonzero");



  if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY)
  {
    if(verbose >= 1)
    {
      if(!(NM_is_symmetric(M)))
      {
        double d= NM_symmetry_discrepancy(M);
        numerics_printf_verbose(1,"fc3d_admm ---- FC3D - ADMM - M is not symmetric (%e) but fc3d_admm_symmetric  \nis called",d);
      }
    }
    fc3d_admm_symmetric(problem, reaction, velocity, info, options, rho,  is_rho_variable, norm_q );
  }
  else if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY)
  {
    fc3d_admm_asymmetric(problem, reaction, velocity, info, options, rho,  is_rho_variable, norm_q );
  }
  else if (options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY)
  {
    if(!(NM_is_symmetric(M)))
    {
      /* double d= NM_symmetry_discrepancy(M); */
      /* numerics_warning("fc3d_admm","---- FC3D - ADMM - M is not symmetric (%e) but we assume it\n",d); */
       fc3d_admm_asymmetric(problem, reaction, velocity, info, options, rho,  is_rho_variable, norm_q );
    }
    else
    {
      fc3d_admm_symmetric(problem, reaction, velocity, info, options, rho,  is_rho_variable, norm_q );
    }

  }
  else
    numerics_error("fc3d_admm", "iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = %i is not implemented", options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY]);

  numerics_printf_verbose(1,"---- FC3D - ADMM - Solution information");
  numerics_printf_verbose(1,"---- FC3D - ADMM - norm of velocity = %e, norm of q = %e ",
                          cblas_dnrm2(m , velocity , 1) , norm_q);
  numerics_printf_verbose(1,"---- FC3D - ADMM - norm of reaction = %e",
                          cblas_dnrm2(m , reaction , 1) , norm_q);

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
  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY;
  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE] =  SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] =
    SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
    SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO] =
    SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options->dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 1.0;
  options->dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA] = 0.999;
  options->dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU]=2.0;
  options->dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI]=2.0;

  options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]=SICONOS_FRICTION_3D_RESCALING_NO;


  options->internalSolvers = NULL;
  options->solverData = NULL;

  return 0;
}
