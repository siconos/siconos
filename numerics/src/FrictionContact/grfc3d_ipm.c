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

#include "CSparseMatrix_internal.h"
#include "global_rolling_fc_Solvers.h"  // for GRFCProb, SolverOpt, Friction_cst, grfc3d_...
#include "gfc3d_compute_error.h"
#include "SiconosLapack.h"

#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "float.h"
#include "JordanAlgebra.h"              // for JA functions

#include "NumericsSparseMatrix.h"
#include "NumericsMatrix.h"

#include "projectionOnCone.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

const char* const SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM_STR = "GRFC3D_IPM"; // for printing the name of solver

typedef struct
{
  double * globalVelocity;
  double * reaction;
  double * velocity;
}
  IPM_point;

typedef struct
{
  NumericsMatrix* mat;
  NumericsMatrix* inv_mat;
}
  IPM_change_of_variable;

typedef struct
{
  double alpha_primal; // primal step length
  double alpha_dual;   // dual step length
  double sigma;        // centering parameter
  double barr_param;   // barrier parameter
}
  IPM_internal_params;

typedef struct
{
  /* initial interior points */
  IPM_point* starting_point;
  IPM_point* tmp_point;

  /* change of variable matrix */
  IPM_change_of_variable* P_mu;

  /* initial internal solver parameters */
  IPM_internal_params* internal_params;

  double **tmp_vault_nd;
  double **tmp_vault_m;
}
  Grfc3d_IPM_data;






/* optimization method */
void grfc3d_IPM(GlobalRollingFrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, double* restrict globalVelocity,
               int* restrict info, SolverOptions* restrict options)
{
  printf("\n#################### grfc3d_ipm.c OK ####################\n");
}


/* initialize solver (allocate memory) */
void grfc3d_IPM_init(GlobalRollingFrictionContactProblem* problem, SolverOptions* options)
{
//   unsigned int m = problem->M->size0;
//   unsigned int nd = problem->H->size1;
//   unsigned int d = problem->dimension;

//   if(!options->dWork || options->dWorkSize != (size_t)(m + nd + nd))
//   {
//     options->dWork = (double*)calloc(m + nd + nd, sizeof(double));
//     options->dWorkSize = m + nd + nd;
//   }


//   /* ------------- initialize starting point ------------- */
//   options->solverData=(Grfc3d_IPM_data *)malloc(sizeof(Grfc3d_IPM_data));
//   Grfc3d_IPM_data * data = (Grfc3d_IPM_data *)options->solverData;

//   /* --------- allocate memory for tmp point ----------- */
//   data->tmp_point = (IPM_point*)malloc(sizeof(IPM_point));
//   data->tmp_point->globalVelocity = (double*)calloc(m, sizeof(double));
//   data->tmp_point->velocity = (double*)calloc(nd, sizeof(double));
//   data->tmp_point->reaction = (double*)calloc(nd, sizeof(double));

//   /* 1. v */
//   data->starting_point = (IPM_point*)malloc(sizeof(IPM_point));
//   data->starting_point->globalVelocity = (double*)calloc(m, sizeof(double));
//   for(unsigned int i = 0; i < m; ++ i)
//     data->starting_point->globalVelocity[i] = 0.01;

//   /* 2. u */
//   data->starting_point->velocity = (double*)calloc(nd, sizeof(double));
//   for(unsigned int i = 0; i < nd; ++ i)
//   {
//     data->starting_point->velocity[i] = 0.001;
//     if(i % d == 0)
//       data->starting_point->velocity[i] = 3.0;
//   }

//   /* 3. r */
//   data->starting_point->reaction = (double*)calloc(nd, sizeof(double));
//   for(unsigned int i = 0; i < nd; ++ i)
//   {
//     data->starting_point->reaction[i] = 0.04;
//     if(i % d == 0)
//       data->starting_point->reaction[i] = 0.5;
//   }

//   /* ------ initialize the change of variable matrix P_mu ------- */
//   data->P_mu = (IPM_change_of_variable*)malloc(sizeof(IPM_change_of_variable));
//   data->P_mu->mat = NM_create(NM_SPARSE, nd, nd);
//   NM_triplet_alloc(data->P_mu->mat, nd);
//   data->P_mu->mat->matrix2->origin = NSM_TRIPLET;
//   for(unsigned int i = 0; i < nd; ++i)
//     if(i % d == 0)
//       /* NM_entry(data->P_mu->mat, i, i, 1. / problem->mu[(int)(i/d)]); */
//       NM_entry(data->P_mu->mat, i, i, 1.);
//     else
//       /* NM_entry(data->P_mu->mat, i, i, 1.); */
//       NM_entry(data->P_mu->mat, i, i, problem->mu[(int)(i/d)]);

// // here, need to change the indices of mu
// // then add mu_r






//   /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
//   data->P_mu->inv_mat = NM_create(NM_SPARSE, nd, nd);
//   NM_triplet_alloc(data->P_mu->inv_mat, nd);
//   data->P_mu->inv_mat->matrix2->origin = NSM_TRIPLET;
//   for(unsigned int i = 0; i < nd; ++i)
//     if(i % d == 0)
//       /* NM_entry(data->P_mu->inv_mat, i, i, problem->mu[(int)(i/d)]); */
//       NM_entry(data->P_mu->inv_mat, i, i, 1.);
//     else
//       /* NM_entry(data->P_mu->inv_mat, i, i, 1.); */
//       NM_entry(data->P_mu->inv_mat, i, i, 1.0/problem->mu[(int)(i/d)]);
//   /* ------ initial parameters initialization ---------- */
//   data->internal_params = (IPM_internal_params*)malloc(sizeof(IPM_internal_params));
//   data->internal_params->alpha_primal = 1.0;
//   data->internal_params->alpha_dual = 1.0;
//   data->internal_params->sigma = 0.1;
//   data->internal_params->barr_param = 1.0;


//   /* ----- temporary vaults initialization ------- */
//   data->tmp_vault_nd = (double**)malloc(17 * sizeof(double*));
//   for(unsigned int i = 0; i < 17; ++i)
//     data->tmp_vault_nd[i] = (double*)calloc(nd, sizeof(double));

//   data->tmp_vault_m = (double**)malloc(2 * sizeof(double*));
//   for(unsigned int i = 0; i < 2; ++i)
//     data->tmp_vault_m[i] = (double*)calloc(m, sizeof(double));

}



/* deallocate memory */
void grfc3d_IPM_free(GlobalRollingFrictionContactProblem* problem, SolverOptions* options)
{

}


/* setup default solver parameters */
void grfc3d_IPM_set_default(SolverOptions* options)
{

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] = SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  /* 0: convex case;  1: non-smooth case */
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 1;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] = 1;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING] = 0;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.095; //0.095

}
