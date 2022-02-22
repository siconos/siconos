/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include <math.h>                          // for fabs, sqrt, fmax, INFINITY
#include <stdio.h>                         // for NULL, printf
#include <stdlib.h>                        // for calloc, free, malloc
#include "CSparseMatrix_internal.h"
#include "Friction_cst.h"                  // for SICONOS_FRICTION_3D_ADMM_I...
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NumericsFwd.h"                   // for SolverOptions, GlobalFrict...
#include "NumericsMatrix.h"                // for NM_gemv, NumericsMatrix
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "float.h"                         // for DBL_EPSILON
#include "gfc3d_Solvers.h"                 // for gfc3d_checkTrivialCaseGlobal
#include "gfc3d_compute_error.h"           // for gfc3d_compute_error
#include "numerics_verbose.h"              // for numerics_printf_verbose
#include "projectionOnCone.h"              // for projectionOnDualCone
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_dscal
#include "NumericsSparseMatrix.h"                // for NSM_TRIPLET ...
#include "gfc3d_balancing.h"

//#define DEBUG_NOCOLOR
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif


const char* const   SICONOS_GLOBAL_FRICTION_3D_ADMM_STR = "GFC3D ADMM";

typedef struct
{
  double * reaction_hat;
  double * reaction_k;
  double * u_hat;
  double * u_k;
  double * u;
  double * b_full;
  double * u_old;
  double * sliding_direction;
  double * sliding_direction_old;
}
Gfc3d_ADDM_data;

/** pointer to function used to call local solver */
typedef int (*LinearSolverPtr)(NumericsMatrix *M, double *b, unsigned int nrhs);




void gfc3d_ADMM_init(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  size_t nc = problem->numberOfContacts;
  size_t n = problem->M->size0;
  size_t m = 3 * nc;
  if(!options->dWork || options->dWorkSize != m+n)
  {
    options->dWork = (double*)calloc(m+n,sizeof(double));
    options->dWorkSize = m+n;
  }
  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION ||
      options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART ||
      options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION) /* Could be optimized */
  {
    options->solverData=(Gfc3d_ADDM_data *)malloc(sizeof(Gfc3d_ADDM_data));
    Gfc3d_ADDM_data * data = (Gfc3d_ADDM_data *)options->solverData;
    data->reaction_hat = (double*)calloc(m,sizeof(double));
    data->reaction_k = (double*)calloc(m,sizeof(double));
    data->u_hat = (double*)calloc(m,sizeof(double));
    data->u_k = (double*)calloc(m,sizeof(double));
    data->u = (double*)calloc(m,sizeof(double));
    data->b_full = (double*)calloc(m,sizeof(double));
    if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H] ==
        SICONOS_FRICTION_3D_ADMM_FULL_H_YES)
    {
      data->u_old = (double*)calloc(m,sizeof(double));
      data->sliding_direction = (double*)calloc(2*nc,sizeof(double));
      data->sliding_direction_old = (double*)calloc(2*nc,sizeof(double));
    }
  }
}
void gfc3d_ADMM_free(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  if(options->dWork)
  {
    free(options->dWork);
    options->dWorkSize = 0;
  }
  options->dWork=NULL;

  if(options->solverData)
  {
    Gfc3d_ADDM_data * data = (Gfc3d_ADDM_data *)options->solverData;
    free(data->reaction_hat);
    free(data->u_hat);
    free(data->reaction_k);
    free(data->u_k);
    free(data->b_full);
    if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H] ==
        SICONOS_FRICTION_3D_ADMM_FULL_H_YES)
    {
      free(data->u_old);
      free(data->sliding_direction);
      free(data->sliding_direction_old);
    }
    free(data);
  }
  options->solverData = NULL;
}
static double gfc3d_admm_select_rho(NumericsMatrix* M, NumericsMatrix* H, int * is_rho_variable, SolverOptions* restrict options)
{
  double rho=0.0;


  /* initial rho */
  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] ==
      SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN)
  {
    rho = options->dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }
  else if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] ==
          SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF)
  {
    double norm_1_M =   NM_norm_1(M);
    double norm_1_H =   NM_norm_1(H);
    if((fabs(norm_1_H) > DBL_EPSILON) && (fabs(norm_1_M) > DBL_EPSILON))
      rho = norm_1_M/norm_1_H;
    else
      rho = options->dparam[SICONOS_FRICTION_3D_ADMM_RHO];
  }
  else  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] ==
           SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES)
  {
    double lambda_max =  NM_iterated_power_method(M, 1e-08, 100);
    double lambda_min =  1.0/NM_iterated_power_method(NM_LU_inv(M), 1e-08, 100);

    numerics_printf_verbose(1,"---- GFC3D - ADMM - largest eigenvalue of M = %g ",lambda_max);
    numerics_printf_verbose(1,"---- GFC3D - ADMM - smallest eigenvalue of M = %g ",lambda_min);

    rho = sqrt(lambda_max*lambda_min);
  }
  else
  {
    numerics_error("gfc3d_admm_select_rho", "options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] unknow value");
  }

  /* rho adaptive from initial rho */
  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
      SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING||
      options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
      SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
  {
    *is_rho_variable = 1 ;
  }
  else
    *is_rho_variable = 0 ;
  return rho;
}


static inline void compute_sliding_direction(int nc, double * u, double * sliding_direction)
{
  for(int contact =0; contact < nc; contact ++)
  {
    int pos = contact*3;
    double normUT = sqrt(u[pos + 1] * u[pos + 1] +
                         u[pos + 2] * u[pos + 2]);
    if(normUT > DBL_EPSILON)
    {
      sliding_direction[contact*2]   =  u[pos+1]/normUT;
      sliding_direction[contact*2+1] =  u[pos+2]/normUT;
    }
    else
    {
      sliding_direction[contact*2]   =  0.0;
      sliding_direction[contact*2+1] =  0.0;
    }
  }
}

static inline void gfc3d_ADMM_compute_full_b(
  int nc, double * u,
  double * mu,
  double * b,
  double * b_s,
  SolverOptions *options,
  int update_b)

{
  if(update_b)
  {
    cblas_dcopy(3*nc, b, 1, b_s,1);
    if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]==
        SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
    {
      for(int contact = 0 ; contact < nc ; ++contact)
      {
        int pos = contact * 3;
        b_s[pos] +=  mu[contact]* sqrt(u[pos + 1] * u[pos + 1] + u[pos + 2] * u[pos + 2]);
      }
    }
  }
}

static inline void gfc3d_ADMM_compute_full_H(int nc, double * u,
    double * mu,
    NumericsMatrix * H,
    NumericsMatrix * H_full)
{

  NumericsMatrix* H_correction = NM_create(NM_SPARSE,3*nc,3*nc);
  NM_triplet_alloc(H_correction, 3*nc);
  H_correction->matrix2->origin = NSM_TRIPLET;

  for(int contact =0; contact < nc; contact ++)
  {
    int pos = contact * 3;
    double normUT = sqrt(u[pos + 1] * u[pos + 1] +
                         u[pos + 2] * u[pos + 2]);

    if(normUT > DBL_EPSILON)
    {
      NM_entry(H_correction,pos,pos+1, mu[contact]* u[pos+1]/normUT);
      NM_entry(H_correction,pos,pos+2, mu[contact]* u[pos+2]/normUT);
    }
    else
    {
      /* This really set entry to 0: note that with old NM_zentry
       there were a bug here as 0.0 was not added. NM_entry may be
       improved to remove entry from sparse matrices, but would it be
       more efficient? */
      NM_entry(H_correction,pos,pos+1, 0.0);
      NM_entry(H_correction,pos,pos+2, 0.0);
    }
    NM_entry(H_correction,pos,pos, 1.0);
    NM_entry(H_correction,pos+1,pos+1, 1.0);
    NM_entry(H_correction,pos+2,pos+2, 1.0);

  }
  DEBUG_EXPR(NM_display(H_correction););

  //DEBUG_EXPR(NM_display(H));
  NM_gemm(1.0, H, H_correction, 0.0, H_full);
  NM_clear(H_correction);

  DEBUG_EXPR(NM_display(H_full));
  DEBUG_END("gfc3d_ADMM_compute_H_correction(...)\n");
  /* getchar();  */
}


static void gfc3d_print_problem_info(GlobalFrictionContactProblem* restrict problem, SolverOptions* restrict options)
{
  NumericsMatrix *M = problem->M;
  NumericsMatrix *H = problem->H;
  size_t nc = problem->numberOfContacts;
  size_t n = problem->M->size0;
  size_t m = 3 * nc;
  double* q = problem->q;
  double* b = problem->b;
  double norm_q = cblas_dnrm2(n, q, 1);
  double norm_b = cblas_dnrm2(m, b, 1);

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO] ==
     SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES)
  {
    numerics_printf_verbose(1,"---- GFC3D - ADMM - Problem information");
    numerics_printf_verbose(1,"---- GFC3D - ADMM - 1-norm of M = %g norm of q = %g ", NM_norm_1(M), norm_q);
    numerics_printf_verbose(1,"---- GFC3D - ADMM - inf-norm of M = %g ", NM_norm_inf(M));
    double eig_max = NM_iterated_power_method(M, 1e-08, 100);
    double eig_min =  1./NM_iterated_power_method(NM_LU_inv(M), 1e-08, 100);
    numerics_printf_verbose(1,"---- GFC3D - ADMM - largest eigenvalue of M = %g ", eig_max);
    numerics_printf_verbose(1,"---- GFC3D - ADMM - smallest eigenvalue of M = %g ", eig_min);
    numerics_printf_verbose(1,"---- GFC3D - ADMM - conditioning of M = %g ", eig_max/eig_min);

    numerics_printf_verbose(1,"---- GFC3D - ADMM - 1-norm of H = %g norm of b = %g ", NM_norm_1(H), norm_b);
    numerics_printf_verbose(1,"---- GFC3D - ADMM - inf-norm of H = %g ", NM_norm_inf(H));
    /* W = NM_create(NM_SPARSE,n,n); */
    /* NM_triplet_alloc(W, n); */
    /* W->matrix2->origin = NSM_TRIPLET; */
    /* NM_gemm(1.0, H, Htrans, 0.0, W); */
    /* numerics_printf_verbose(1,"---- GFC3D - ADMM - 1-norm of HH^T = %g  ", NM_norm_1(W)); */
    /* numerics_printf_verbose(1,"---- GFC3D - ADMM - largest eigenvalue of HH^T = %g ", NM_iterated_power_method(W, 1e-08, 100)); */
    /* numerics_printf_verbose(1,"---- GFC3D - ADMM - smallest eigenvalue of HH^T = %g ", NM_iterated_power_method(NM_inv(W), 1e-08, 100)); */
    if(NM_is_symmetric(M))
    {
      numerics_printf_verbose(1,"---- GFC3D - ADMM -  M is symmetric");
    }
    else
    {
      numerics_printf_verbose(1,"---- GFC3D - ADMM -  M is not symmetric");
    }
    //NM_clear(W);
  }
}

void gfc3d_ADMM(GlobalFrictionContactProblem* restrict problem_original, double* restrict reaction,
                double* restrict velocity, double* restrict globalVelocity,
                int* restrict info, SolverOptions* restrict options)
{
  /* verbose=2; */

  int* iparam = options->iparam;
  double* dparam = options->dparam;
  size_t nc = problem_original->numberOfContacts;
  size_t n = problem_original->M->size0;
  size_t m = 3 * nc;

  /**************************************************************************/
  /* Balancing                        ***************************************/
  /**************************************************************************/

  GlobalFrictionContactProblem* problem = gfc3d_balancing_problem(problem_original,options);
  gfc3d_balancing_go_to_balanced_variables(problem, options,
                                           reaction, velocity, globalVelocity);
  GlobalFrictionContactProblem* original_problem=NULL;
  if(iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
  {
    GlobalFrictionContactProblem_balancing_data  *data = (GlobalFrictionContactProblem_balancing_data * ) problem->env;
    original_problem = data->original_problem;
    assert(original_problem);
    original_problem->norm_b = cblas_dnrm2(m, original_problem->b, 1);
    original_problem->norm_q = cblas_dnrm2(n, original_problem->q, 1);
  }



  double* q = problem->q;
  double* b = problem->b;
  double* mu = problem->mu;


  NumericsMatrix* M_original=problem->M;
  NumericsMatrix* H_original=problem->H;

  /**************************************************************************/
  /* Change storage                        **********************************/
  /**************************************************************************/


  if(iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE)
  {
    if (problem->M->storageType == NM_SPARSE_BLOCK)
    {
      DEBUG_PRINT("Force a copy to sparse storage type\n");
      problem->M = NM_create(NM_SPARSE,  problem->M->size0,  problem->M->size1);
      NM_copy_to_sparse(M_original, problem->M, DBL_EPSILON);
    }

    if(problem->H->storageType == NM_SPARSE_BLOCK)
    {
      DEBUG_PRINT("Force a copy to sparse storage type\n");
      problem->H = NM_create(NM_SPARSE,  problem->H->size0,  problem->H->size1);
      NM_copy_to_sparse(H_original, problem->H, DBL_EPSILON);
    }
  }


  /**************************************************************************/
  /***************** Strategies for dealing with symmetry *******************/
  /**************************************************************************/

  DEBUG_PRINT("Strategies for dealing with symmetry \n");
  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = SICONOS_FRICTION_3D_ADMM_SYMMETRIZE;
  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY;

  LinearSolverPtr linear_solver;
  NumericsMatrix *Msym = NULL;

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY ||
     options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY)
  {
    numerics_error("gfc3d_admm", "iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = %i is not implemented\n only the forced symmetry is implemented",
                   options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY]);
  }
  else if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY)
  {
    /* The symmetric version of the algorithm is used even if
     *  the system is not symmetric using the LU solver */
    if(verbose >= 1)
    {
      if(!(NM_is_symmetric(problem->M)))
      {
        double d= NM_symmetry_discrepancy(problem->M);
        numerics_printf_verbose(1,"gfc3d_admm ---- GFC3D - ADMM - M is not symmetric (%e) but gfc3d_admm_symmetric  \nis called with LU solver",d);
      }
    }
    linear_solver = &NM_LU_solve;
  }
  else if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_SYMMETRIZE)
  {
    /* The symmetric version of the algorithm is used and the matrix
     *is systematically symmetrized*/
    NumericsMatrix *MT = NM_transpose(problem->M);
    Msym = NM_add(1/2., problem->M, 1/2., MT );
    //NM_display(Msym);
    //getchar();
    problem->M = Msym;
    NM_clear(MT);
    linear_solver = & NM_Cholesky_solve;
  }
  else if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_ASSUME_SYMMETRY)
  {
    /* The symmetric version of the algorithm is used and we assume
     *  that the data are symmetric */
    linear_solver = & NM_Cholesky_solve;
  }
  else
    numerics_error("gfc3d_admm", "iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = %i is not implemented", options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY]);



  assert((int)problem->H->size1 == problem->numberOfContacts * problem->dimension);

  /**************************************************************************/
  /***************** Rescling cone    ******************* *******************/
  /**************************************************************************/



  /* /\* rescale problem on cone *\/ */
  double cone_scaling=1.0;
  int rescaling_cone=0;
  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE]==SICONOS_FRICTION_3D_RESCALING_CONE_YES)
  {
    rescaling_cone=1;
    cone_scaling=0.1;
    numerics_printf_verbose(2, "The second order cone is rescaled such that mu = %f",cone_scaling);
    NumericsMatrix * P = NM_create(NM_SPARSE,m,m);
    NM_triplet_alloc(P, m);
    P->matrix2->origin = NSM_TRIPLET;
    mu = (double *)malloc(nc*sizeof(double));
    b = (double *)malloc(m*sizeof(double));
    cblas_dcopy(m, problem->b, 1, b, 1);

    for(size_t contact = 0 ; contact < nc ; ++contact)
    {
      int pos = contact*3;
      NM_entry(P,pos,pos, cone_scaling/problem->mu[contact]);
      NM_entry(P,pos+1,pos+1, 1.0);
      NM_entry(P,pos+2,pos+2, 1.0);
      b[pos]=cone_scaling/problem->mu[contact]*b[pos];
      mu[contact]=cone_scaling;
    }
    NumericsMatrix *Htrans =  NM_transpose(problem->H);
    NumericsMatrix *H = problem->H;
    H = NM_create(NM_SPARSE,n,m);
    NM_triplet_alloc(H, n);
    H->matrix2->origin = NSM_TRIPLET;
    NM_copy(problem->H, H);
    NM_gemm(1.0, P, Htrans, 0.0, Htrans);
    NM_gemm(1.0, H, P, 0.0, H);
  }


  /* storage for W = M + rho H H^T */

  NumericsMatrix *W = NM_create(NM_SPARSE,n,n);
  NM_triplet_alloc(W, n);
  W->matrix2->origin = NSM_TRIPLET;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  /* Check for trivial case */
  *info = gfc3d_checkTrivialCaseGlobal(n, q, velocity, reaction, globalVelocity, options);

  if(*info == 0)
    return;


  double norm_q = cblas_dnrm2(n, q, 1);
  problem->norm_q=norm_q;
  double norm_b = cblas_dnrm2(m, b, 1);
  problem->norm_b=norm_b;

  int internal_allocation=0;
  if(!options->dWork || options->dWorkSize != 2*m+n)
  {
    gfc3d_ADMM_init(problem, options);
    internal_allocation = 1;
  }

  NumericsMatrix *H_full = NM_create(NM_SPARSE,n,m);
  NM_triplet_alloc(H_full, n);
  H_full->matrix2->origin = NSM_TRIPLET;

  //getchar();
  int is_rho_variable=0;
  double rho = gfc3d_admm_select_rho(problem->M, problem->H,  &is_rho_variable, options);
  if(rho <= DBL_EPSILON)
    numerics_error("gfc3d_ADMM", "dparam[SICONOS_FRICTION_3D_ADMM_RHO] must be nonzero");

  double eta = dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA];
  double br_tau = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU];
  double br_phi = dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI];

  Gfc3d_ADDM_data * data = (Gfc3d_ADDM_data *)options->solverData;


  double * v = globalVelocity;

  double * u = data->u;
  double * u_k = data->u_k;
  double * u_hat =  data->u_hat;

  double * reaction_k =  data->reaction_k;
  double * reaction_hat = data->reaction_hat;



  double * tmp_m =  options->dWork;
  double * tmp_n =  &options->dWork[m];

  cblas_dscal(m, 1.0/rho, reaction, 1);

  cblas_dcopy(m, reaction, 1, reaction_k, 1);
  cblas_dcopy(m, u, 1, u_k, 1);

  cblas_dcopy(m, reaction, 1, reaction_hat, 1);
  cblas_dcopy(m, u, 1, u_hat, 1);

  double rho_k=0.0, rho_ratio=0.0;
  double e_k = INFINITY, e, alpha, r, s, residual, r_scaled, s_scaled;
  double norm_b_full=0.0;
  double tau, tau_k = 1.0;

  rho_k=rho;
  int has_rho_changed = 1;

  int with_full_Jacobian = options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H];
  int has_full_H_changed = 1;

  double * b_full = data->b_full;

  double * sliding_direction = data->sliding_direction;
  double * sliding_direction_old =  data->sliding_direction_old;
  double * u_old =  data->u_old;
  if(with_full_Jacobian)
    cblas_dcopy(m,u,1,u_old,1);



  NumericsMatrix* M = problem->M;
  NumericsMatrix* H = problem->H;
  NumericsMatrix* Htrans = NM_transpose(problem->H);;


  int update_b =1;
  ComputeErrorGlobalPtr computeError = NULL;

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]==
      SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
  {
    computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;
  }
  else
  {
    computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error_convex;
  }

  /* double * normUT_old  = (double *) malloc(nc*sizeof(double)); */
  /* double * normUT_current  = (double *) malloc(nc*sizeof(double)); */
  /* double delta_normUT =0.0; */
  /* projection. loop through the contact points */
  /* for (int contact = 0 ; contact < nc ; ++contact) */
  /* { */
  /*   problem->mu[contact]=0.6; */
  /* } */

  gfc3d_print_problem_info(problem, options);



  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  assert(rho>0);
  /**************************************************************************/
  /***********************************************  ADMM iterations *********/
  /**************************************************************************/
  while((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /********************/
    /*  0 - Compute b   */
    /********************/
    gfc3d_ADMM_compute_full_b(nc, u, mu, b, b_full, options, update_b);

    /********************/
    /*  1 - Compute v */
    /********************/
    if(with_full_Jacobian || has_rho_changed)
    {

      if(with_full_Jacobian)
      {

        /* for (int contact = 0 ; contact < nc ; ++contact) */
        /* { */
        /*   int pos = contact * 3; */
        /*   normUT_old[contact] = sqrt(u_old[pos + 1] * u_old[pos + 1] + u_old[pos + 2] * u_old[pos + 2]); */
        /*   normUT_current[contact] = sqrt(u[pos + 1] * u[pos + 1] + u[pos + 2] * u[pos + 2]); */
        /* } */
        /* cblas_daxpy(nc, -1.0, normUT_current, 1, normUT_old, 1); */
        /* delta_normUT = cblas_dnrm2(nc, normUT_old, 1); */
        /* printf("delta_normUT = %8.4e\n", delta_normUT); */
        /* printf("normUT = %8.4e\n", cblas_dnrm2(nc, normUT_current, 1)); */




        /* test if we need to update H_full */
        compute_sliding_direction(nc, u, sliding_direction);
        compute_sliding_direction(nc, u_old, sliding_direction_old);
        cblas_daxpy(2*nc, -1.0, sliding_direction, 1, sliding_direction_old, 1);
        double  delta_sliding_direction= cblas_dnrm2(2*nc, sliding_direction_old, 1);
        /* printf("delta_sliding_direction = %8.4e\n", delta_sliding_direction); */
        /* printf("norm_sliding_direction = %8.4e\n", cblas_dnrm2(2*nc, sliding_direction, 1)); */

        double criteria = 0.0;
        if(cblas_dnrm2(2*nc, sliding_direction, 1) > DBL_EPSILON)
          criteria= fabs(delta_sliding_direction) / cblas_dnrm2(2*nc, sliding_direction, 1);

        if(criteria > 1e-01)
        {
          /* printf("Recompute H_full\n"); */
          gfc3d_ADMM_compute_full_H(nc, u, mu, H, H_full);
          cblas_dcopy(3*nc,u,1,u_old,1);
          has_full_H_changed =1;
        }
        else
        {
          has_full_H_changed =0;
        }
        if(has_full_H_changed || has_rho_changed)
        {
          NM_copy(M, W);
          NM_unpreserve(W);
          NM_gemm(rho, H_full, Htrans, 1.0, W);
        }
      }
      else
      {
        NM_copy(M, W);
        NM_unpreserve(W); /* if not unpreserve, the follwing operations are 
                             not made on the destructible pointer */
        NM_gemm(rho, H, Htrans, 1.0, W);
      }
      DEBUG_PRINT("M + rho H H^T: ");
      DEBUG_EXPR(NM_display(W));
    }

    /* compute the rhs */
    /* q --> v */
    cblas_dcopy(n, q, 1, v, 1);

    /* q +  rho H*( u -b + reaction_k) --> v */
    cblas_dcopy(m, u_hat, 1, tmp_m, 1);
    gfc3d_ADMM_compute_full_b(nc, u_hat, mu, b, b_full, options, update_b);
    cblas_daxpy(m, -1.0, b_full, 1, tmp_m, 1);

    if(with_full_Jacobian)
    {
      NM_gemv(rho, H_full, tmp_m, 1.0, v);
      NM_gemv(rho, H, reaction_hat, 1.0, v);
    }
    else
    {
      cblas_daxpy(m, 1.0, reaction_hat, 1, tmp_m, 1);
      NM_gemv(rho, H, tmp_m, 1.0, v);
    }

    /* Linear system solver */
    /* cblas_dcopy(n , w_k , 1 , v, 1); */
    if(with_full_Jacobian)
    {
      /* W destroyed */
      // NM_gesv_expert(W,v,NM_KEEP_FACTORS);
      NM_LU_solve(W, v, 1);
    }

    else
    {
      linear_solver(W,v,1);
    }

    /********************/
    /*  2 - Compute u */
    /********************/

    /* H^T v_k - reaction_k + b */
    gfc3d_ADMM_compute_full_b(nc, u, mu, b,  b_full, options, 0);

    cblas_dcopy(m, b_full, 1, u, 1);
    cblas_daxpy(m, -1.0, reaction_hat, 1, u, 1);
    NM_gemv(1.0, Htrans, v, 1.0, u);

    DEBUG_PRINT("before projection");
    DEBUG_EXPR(NV_display(u,m));

    /* projection. loop through the contact points */
    for(size_t contact = 0 ; contact < nc ; ++contact)
    {
      int pos = contact * 3;
      projectionOnDualCone(&u[pos], mu[contact]);
    }

    double norm_u =  cblas_dnrm2(m, u, 1);
    DEBUG_EXPR(NV_display(u,m));

    /*************************/
    /*  3 - Compute reaction */
    /*************************/


    /* - H^T v_k + u_k -b_full ->  reaction (We use reaction for storing the residual for a while) */
    /* cblas_dscal(m, 0.0, reaction, 1);  */
    NM_gemv(-1.0, Htrans, v, 0.0, reaction);
    double norm_HTv = cblas_dnrm2(m, reaction, 1);

    gfc3d_ADMM_compute_full_b(nc, u, mu, b, b_full, options, 0);
    norm_b_full =  cblas_dnrm2(m, b_full, 1);

    cblas_daxpy(m, -1.0, b_full, 1, reaction, 1);
    cblas_daxpy(m, 1.0, u, 1, reaction, 1);

    /* compute primal residual */
    r = cblas_dnrm2(m, reaction, 1);   /* norm of the residual */

    /* compute the actual reaction */
    /* reaction_hat -  A v_k + u_k -b_full ->  reaction */
    cblas_daxpy(m, 1.0, reaction_hat, 1, reaction, 1);

    /* cblas_dscal(n, 0.0, tmp_n, 1); */
    NM_gemv(1.0*rho, H, reaction, 0.0, tmp_n);
    double norm_rhoHr = cblas_dnrm2(n, tmp_n, 1);

    /*********************************/
    /*  3 - Acceleration and restart */
    /*********************************/

    DEBUG_EXPR(NV_display(u_hat,m));
    DEBUG_EXPR(NV_display(u,m));

    /* Compute dual residual */
    cblas_dcopy(m, u_hat, 1, tmp_m, 1);
    cblas_daxpy(m, -1.0, u, 1, tmp_m, 1);

    /* cblas_dscal(n, 0.0, tmp_n, 1); */
    double s_restart =  rho * cblas_dnrm2(m, tmp_m, 1);

    NM_gemv(1.0*rho, H, tmp_m, 0.0, tmp_n);
    s = cblas_dnrm2(n, tmp_n, 1);


    /* cblas_dcopy(m , u_k , 1 , tmp_m, 1); */
    /* cblas_daxpy(m, -1.0, u, 1, tmp_m , 1); */
    /* cblas_dscal(n, 0.0, tmp_n, 1); */
    /* NM_gemv(1.0*rho, H, tmp_m, 1.0, tmp_n); */
    /* double s_k = cblas_dnrm2(n , tmp_n , 1); */

    /*  Compute full residual for restart */
    e =r*r+s*s;
    e =r*r+s_restart*s_restart;

    DEBUG_PRINTF("residual e = %e \n", e);
    DEBUG_PRINTF("residual r = %e \n", r);
    DEBUG_PRINTF("residual s = %e \n", s);
    DEBUG_PRINTF("residual e_k = %e \n", e_k);
    /* DEBUG_PRINTF("residual s_k = %e \n", s_k); */
    DEBUG_PRINTF("eta  = %e \n", eta);
    if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION ||
        options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART)
    {
      if(e <  eta * e_k)
      {
        tau  = 0.5 *(1 +sqrt(1.0+4.0*tau_k*tau_k));
        alpha = (tau_k-1.0)/tau;

        cblas_dcopy(m, u, 1, u_hat, 1);
        cblas_dscal(m, 1+alpha, u_hat,1);
        cblas_daxpy(m, -alpha, u_k, 1, u_hat, 1);
        DEBUG_EXPR(NV_display(u_hat,m));

        cblas_dcopy(m, reaction, 1, reaction_hat, 1);
        cblas_dscal(m, 1+alpha, reaction_hat,1);
        cblas_daxpy(m, -alpha, reaction_k, 1, reaction_hat, 1);
        DEBUG_EXPR(NV_display(reaction_hat,m));
        DEBUG_PRINTF("Accelerate tau  = %e, \t tau_k  = %e \t alpha  = %e   \n", tau, tau_k, alpha);
        numerics_printf_verbose(2, "Accelerate :tau  = %e, \t tau_k  = %e, \t alpha  = %e ", tau, tau_k, alpha);
        tau_k=tau;
        e_k=e;
      }
      else
      {
        tau_k=1.0;
        e_k = e_k /eta;
        /* e_k=e; */
        DEBUG_PRINTF("Restart tau_k  = %e \n", tau_k);
        numerics_printf_verbose(2,"Restart tau_k  = %e, e = %e\t, e_k = %e, e/e_k = %e\t", tau_k, e, e_k, e/e_k);
        cblas_dcopy(m, reaction_k, 1, reaction_hat, 1);
        cblas_dcopy(m, u_k, 1, u_hat, 1);
      }
    }
    else  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] == SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION)
    {
      tau_k=1.0;
      e_k = e_k /eta;
      numerics_printf_verbose(2,"No acceleration tau_k  = %e  \n", tau_k);
      cblas_dcopy(m, reaction_k, 1, reaction_hat, 1);
      cblas_dcopy(m, u_k, 1, u_hat, 1);
    }
    else
    {
      numerics_error("gfc3d_admm", " options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] value is not recognized");
    }


    rho_k = rho ;
    /* numerics_printf_verbose(2, "gfc3d_admm. residuals : r  = %e, \t  s = %e \t s_k = %e", r, s, s_k); */
    numerics_printf_verbose(2, "gfc3d_admm. residuals : r  = %e, \t  s = %e ", r, s);
    numerics_printf_verbose(2, "gfc3d_admm. scaling   : norm_u  = %e, \t norm_HTv  = %e, \t norm_b = %e, norm_b_full = %e, norm_rhoHr = %e\t", norm_u,  norm_HTv, norm_b, norm_b_full, norm_rhoHr);
    if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] ==
        SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
    {
      r_scaled = r / fmax(norm_u,(fmax(norm_HTv, norm_b_full)));
      s_scaled = s / (norm_rhoHr);
    }
    else
    {
      r_scaled = r;
      s_scaled = s;
    }
    numerics_printf_verbose(2, "gfc3d_admm. scaled residuals : r_scaled  = %e, \t  s_scaled = %e", r_scaled, s_scaled);

    if(is_rho_variable)
    {
      if(r_scaled > br_phi * s_scaled)
      {
        rho = br_tau* rho_k;
        has_rho_changed = 1;
      }
      else if(s_scaled > br_phi * r_scaled)
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
    numerics_printf_verbose(2, "gfc3d_admm. rho = %5.2e\t, rho_k = %5.2e\t ", rho, rho_k);
    rho_ratio = rho_k/rho;

    DEBUG_PRINTF("rho =%e\t,rho_k =%e \n", rho, rho_k);

    cblas_dscal(m, rho_ratio, reaction,1);
    cblas_dscal(m, rho_ratio, reaction_hat,1);

    /* Next step */
    cblas_dcopy(m, reaction, 1, reaction_k, 1);
    cblas_dcopy(m, u, 1, u_k, 1);

    /*********************************/
    /*  4 - Stopping criterium       */
    /*********************************/

    int stopping_criterion =0;

    /* old version */
    residual = sqrt(e);
    /* if (fabs(norm_q) > DBL_EPSILON) */
    /*   residual /= norm_q; */
    /* if (residual < tolerance) */
    /*   stopping_criterion =1; */
    /* printf("norm_b_full =%e, fmax(norm_u,(fmax(norm_HTv, norm_b_full)))=%e\n", norm_b_full, fmax(norm_u,(fmax(norm_HTv, norm_b_full)))); */
    /* printf(" =%e, fmax(norm_u,(fmax(norm_HTv, norm_b_full)))=%e\n", norm_b_full, fmax(norm_u,(fmax(norm_HTv, norm_b_full)))); */

    double scaling_error_primal = fmax(norm_u,(fmax(norm_HTv, norm_b_full))) +  sqrt(m);
    double epsilon_primal = tolerance * scaling_error_primal  ;
    double scaling_error_dual = norm_rhoHr + sqrt(n);
    double epsilon_dual =  tolerance * scaling_error_dual;

    if(r <= epsilon_primal && s <= epsilon_dual)
      stopping_criterion =1;

    numerics_printf_verbose(1,"---- GFC3D - ADMM  - Iteration %i rho = %14.7e, full residual (e) = %14.7e, tol = %14.7e", iter, rho, residual, tolerance);
    numerics_printf_verbose(1,"---- GFC3D - ADMM  -                            primal residual = %14.7e, epsilon_primal = %14.7e", r,  epsilon_primal);
    numerics_printf_verbose(1,"---- GFC3D - ADMM  -                            dual residual = %14.7e, epsilon_dual = %14.7e", s,  epsilon_dual);


    int adapt_tolerance = 1;
    /* We check the full criterion at a given frequency */
    /* by forcing the stopping criterion at 1 but keeping the tolerance */
    if((options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY]>0) &&
        (iter % options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] == 0))
    {
      stopping_criterion =1;
      adapt_tolerance = 0;
    }

    if(stopping_criterion)
    {
      /* check the full criterion */
      cblas_dscal(m, rho, reaction, 1);

      if(rescaling_cone)
      {
        for(size_t contact = 0 ; contact < nc ; ++contact)
        {
          int pos = contact*3;
          reaction[pos] = reaction[pos] * cone_scaling / problem->mu[contact];
        }
      }

      //double error_original=0.0;
      if(iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
      {
        (*computeError)(problem,  reaction, velocity, v,  tolerance, options,
                        norm_q, norm_b,  &error);
        /* printf("############ error  = %g\n", error); */

        gfc3d_balancing_back_to_original_variables(problem,
                                                   options,
                                                   reaction, velocity, v);
        (*computeError)(original_problem,  reaction, velocity, v,  tolerance, options,
                        original_problem->norm_q, original_problem->norm_b,  &error);
        //error_original = error;
        /* printf("############ error original = %g\n", error_original); */
      }
      else
      {
        (*computeError)(problem,  reaction, velocity, v,  tolerance, options,
                        norm_q, norm_b,  &error);
        /* printf("############ error  = %g\n", error); */
      }





      numerics_printf_verbose(1,"---- GFC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);


      if(error < dparam[SICONOS_DPARAM_TOL])
      {
        hasNotConverged = 0;
        numerics_printf_verbose(1,"---- GFC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
      }
      else
      {
        numerics_printf_verbose(1,"---- GFC3D - ADMM  - The tolerance on the  residual is not sufficient to reach accuracy (error =  %14.7e)", error);
        if(adapt_tolerance)
        {
          tolerance = tolerance * fmax(epsilon_dual/scaling_error_dual,epsilon_primal/scaling_error_primal)/error;
          numerics_printf_verbose(1,"---- GFC3D - ADMM  - We reduce the tolerance on the residual to %14.7e", tolerance);
        }
        else
        {
          numerics_printf_verbose(1,"---- GFC3D - ADMM  - We keep the tolerance on the residual to %14.7e", tolerance);
        }

        if(rescaling_cone)
        {
          for(size_t contact = 0 ; contact < nc ; ++contact)
          {
            int pos = contact*3;
            reaction[pos] = reaction[pos] *problem->mu[contact]/cone_scaling;
          }
        }
        cblas_dscal(m, 1.0/rho, reaction, 1);

        if(rescaling_cone)
        {
          for(size_t contact = 0 ; contact < nc ; ++contact)
          {
            int pos = contact*3;
            reaction[pos] = reaction[pos]/ problem->mu[contact];
          }
        }
        if(iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
        {
          gfc3d_balancing_go_to_balanced_variables(problem,
                                                   options,
                                                   reaction, velocity, v);
        }

        
      }
      //getchar();
    }
    *info = hasNotConverged;

  }

  if(iter==itermax)
  {
    cblas_dscal(m, rho, reaction, 1);
    if(iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
    {
      gfc3d_balancing_back_to_original_variables(problem,
                                                 options,
                                                 reaction, velocity, v);
      (*computeError)(original_problem,  reaction, velocity, v,  tolerance, options,
                      original_problem->norm_q, original_problem->norm_b,  &error);
    }
    else
    {
      (*computeError)(problem,  reaction, velocity, v,  tolerance, options,
                      norm_q, norm_b,  &error);
    }
    if(error < dparam[SICONOS_DPARAM_TOL])
    {
      *info = 0;
    }
    numerics_printf_verbose(1,"---- GFC3D - ADMM  - Iteration %i rho = %14.7e \t full error = %14.7e", iter, rho, error);
  }

  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /***** Free memory *****/
  problem = gfc3d_balancing_free(problem, options);
  NM_clear(W);
  NM_clear(Htrans);

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] == SICONOS_FRICTION_3D_ADMM_SYMMETRIZE)
  {
    NM_clear(Msym);
  }

  if(internal_allocation)
  {
    gfc3d_ADMM_free(problem,options);
  }
}



void gfc3d_admm_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;

  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 50;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION] =
    SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART;
  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE] =  SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] =
    SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN;
  /* SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF; */
  /* SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES; */

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY] =
    SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT;
  //SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO] =
    SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H] =
    SICONOS_FRICTION_3D_ADMM_FULL_H_NO;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]=
    SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options->dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 0.1;
  options->dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA] = 0.999;
  options->dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU]=2.0;
  options->dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI]=10.0;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY;

  options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]=SICONOS_FRICTION_3D_RESCALING_NO;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE]=SICONOS_FRICTION_3D_RESCALING_CONE_NO;
}
