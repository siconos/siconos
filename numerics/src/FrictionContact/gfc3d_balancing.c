/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "debug.h"                         // for DEBUG_EXPR
#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
#endif


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
  double* q = problem->q;
  double* b = problem->b;
  double* mu = problem->mu;

  NumericsMatrix *M = problem->M;
  NumericsMatrix *H = problem->H;

 
  data->original_problem = problem;

  double alpha_r=0.0, beta_r=0.0;
  BalancingMatrices * B_for_M = NULL;
  BalancingMatrices * B_for_H = NULL;

  NumericsMatrix *Htrans =  NM_transpose(H);
  /* Compute M + rho H H^T (storage in W)*/
  NumericsMatrix *W = NM_create(NM_SPARSE,n,n);
  NM_triplet_alloc(W, n);
  W->matrix2->origin = NSM_TRIPLET;

  
  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_SCALAR)
  {
    alpha_r = NM_norm_inf(M);
    beta_r = NM_norm_inf(H);
    numerics_printf_verbose(1,"---- GFC3D - ADMM - Scalar rescaling of the problem");
    numerics_printf_verbose(1,"---- GFC3D - ADMM - alpha_r = %e\t beta_r= %e\n", alpha_r, beta_r);

    globalFrictionContact_rescaling(rescaled_problem, 1./alpha_r, 1.0/beta_r, 1.0);

    data->alpha= 1.0/alpha_r;
    data->beta=  1.0/beta_r;
    data->gamma= 1.0 ;
  }
  else if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_M)
  {
    numerics_printf_verbose(1,"---- GFC3D - ADMM - Rescaling of the problem by balancing M");
    data->B_for_M  = NM_BalancingMatrices_new(problem->M);
    globalFrictionContact_balancing_M(rescaled_problem, data->B_for_M);
  }
  /* else if (options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]==SICONOS_FRICTION_3D_RESCALING_BALANCING_H) */
  /* { */
  /*   numerics_printf_verbose(1,"---- GFC3D - ADMM - Rescaling of the problem by balancing H"); */
  /*   data->B_for_H  = NM_BalancingMatrices_new(problem->H); */
  /*   globalFrictionContact_balancing_H(rescaled_problem, data->B_for_H); */
  /* } */
  else
  {
    numerics_printf_verbose(1,"---- GFC3D - ADMM - No rescaling of the problem");
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
    else
      numerics_printf_verbose(1,"---- GFC3D - ADMM - rescaling type is not implemented");
      
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
    else
      numerics_printf_verbose(1,"---- GFC3D - ADMM - rescaling type is not implemented");
      
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
