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


#include "NMS.h"

#include <float.h>
#include <assert.h>

#include "SiconosBlas.h"
#include "Newton_methods.h"
#include "SiconosSets.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

int NMS(NMS_data* data_NMS, void* data, functions_LSA* functions, double* restrict z, double* restrict z_N, int force_watchdog_step, int force_d_step_merit_check, double check_ratio)
{
  double theta_iter;
  int n = data_NMS->size;

  /* H could be sparse ... */
  NumericsMatrix* Htmp = data_NMS->H;
  double* H = Htmp->matrix0;
  double* F = NMS_get_F(data_NMS->workspace, n);
  double* F_merit = NMS_get_F_merit(data_NMS->workspace, n);
  double* JacThetaF_merit = NMS_get_JacTheta_F_merit(data_NMS->workspace, n);
  double* workV1 = NMS_get_dir(data_NMS->workspace, n);
  double* workV2 = NMS_get_generic_workV(data_NMS->workspace, n);

  double dotprod, ref_merit, preRHS;
  double alpha_projected_gradient, alpha_watchdog;

  /* first see whether the path search was successful */
  if (force_watchdog_step) goto watchdog_step;

  /* see if we can do a d_step */
  if (data_NMS->n < data_NMS->n_max)
  {
    /* compute the norm of the increment */
    cblas_dcopy(n, z_N, 1, workV1, 1);
    cblas_daxpy(n, -1.0, z, 1, workV1, 1);
    double dist = cblas_dnrm2(n, workV1, 1);

    if (force_d_step_merit_check || (dist <= data_NMS->delta*(1.0 + check_ratio))) /* try a d_step */
    {
      /* see if the nonmonotone criterion is fulfilled */
      /* F should be up-to-date, H and F_merit have to be recomputed */
      functions->compute_H(data, z, F, workV1, workV2, Htmp);
      functions->compute_F_merit(data, z, F, F_merit);
      cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, H, n, F_merit, 1, 0.0, JacThetaF_merit, 1);

      /* Compute the variation of the merit function */
      functions->compute_F(data, z_N, F);
      functions->compute_F_merit(data, z_N, F, F_merit);
      theta_iter = .5 * cblas_ddot(n, F_merit, 1, F_merit, 1);

      /* check if the merit value does not increase too much */
      if (force_d_step_merit_check || theta_iter < data_NMS->merit_incr*data_NMS->ref_merit*(1.0 + check_ratio))
      {
        DEBUG_PRINTF("NMS d_step accepted theta_iter = %2.2e < tol = %2.2e\n", theta_iter, data_NMS->merit_incr*data_NMS->ref_merit);
        /* we accept the step */
        data_NMS->n++;
        data_NMS->delta *= data_NMS->delta_var;
        double dotprod = cblas_ddot(n, JacThetaF_merit, 1, workV1, 1);

        int criterion = check_nmd_criterion(theta_iter, data_NMS->ref_merit, data_NMS->sigma, dotprod);
        if (criterion == 0)
        {
          goto update_checkpoint;
        }
        else
        {
          goto accept_z_N;
        }
      }
      /* else watchdog step */
      DEBUG_PRINTF("NMS d_step rejected theta_iter = %2.2e; tol = %2.2e\n", theta_iter, data_NMS->merit_incr*data_NMS->ref_merit);
    }
    else
    {
      DEBUG_PRINTF("NMS d_step rejected dist = %2.2e; delta = %2.2e\n", dist, data_NMS->delta*(1.0+ check_ratio));
      goto m_step;
    }
  }
  else /* try the m_step */
  {
m_step:
    /* compute the gradient JacThetaF_merit */
    functions->compute_H(data, z, F, workV1, workV2, Htmp);
    functions->compute_F_merit(data, z, F, F_merit);
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, H, n, F_merit, 1, 0.0, JacThetaF_merit, 1);
    /* compute d = z - z_N */
    cblas_dcopy(n, z, 1, workV1, 1);
    cblas_daxpy(n, -1.0, z_N, 1, workV1, 1);

    /* compute the new merit value */
    functions->compute_F(data, z_N, F);
    functions->compute_F_merit(data, z_N, F, F_merit);
    theta_iter = .5 * cblas_ddot(n, F_merit, 1, F_merit, 1);

    double dotprod = cblas_ddot(n, JacThetaF_merit, 1, workV1, 1);

    int criterion = check_nmd_criterion(theta_iter, data_NMS->ref_merit, data_NMS->sigma, dotprod);
    if (criterion == 0)
    {
      DEBUG_PRINT("NMS m_step succeeded !\n");
      goto update_checkpoint;
    }
    DEBUG_PRINT("NMS m_step unsuccessful !\n");
    /* now we have to do the watchdog */
  } /* end try m_step */

watchdog_step:
  /* "worst-case" scenario: we need to perform a pathsearch. There are
   * multiple options here: line search, arc search or backward pathsearch,
   * which starts from the newton point generated by the modified Lemke
   * procedure and we go back in the sequency of t to find the a point yielding
   * a sufficient descent*/

  switch(data_NMS->watchdog_search_type)
  {
    case LINESEARCH:
      /* compute the gradient JacThetaF_merit at the checkpoint z_c(0)*/
      functions->compute_F(data, NMS_checkpoint_0(data_NMS, n), F);
      functions->compute_H(data, NMS_checkpoint_0(data_NMS, n), F, workV1, workV2, Htmp);
      functions->compute_F_merit(data, NMS_checkpoint_0(data_NMS, n), F, F_merit);
      cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, H, n, F_merit, 1, 0.0, JacThetaF_merit, 1);
      /* compute d = z_c(T_k) - z_c(0) */
      cblas_dcopy(n, NMS_checkpoint_T(data_NMS, n), 1, workV1, 1);
      cblas_daxpy(n, -1.0, NMS_checkpoint_0(data_NMS, n), 1, workV1, 1);

      dotprod = cblas_ddot(n, JacThetaF_merit, 1, workV1, 1);
      ref_merit = data_NMS->ref_merit;
      if (dotprod < 0.0) /* check condition on <JacThetaF_merit, z_c(0) - z_c(T_k)>*/
      {
//        preRHS = -data_NMS->sigma*dotprod; /* we expect a plus in the LS function */
        preRHS = -dotprod; /* we expect a plus in the LS function */
      }
      else
      {
//        preRHS = -data_NMS->sigma*ref_merit;
        preRHS = -ref_merit;
      }

      data_NMS->ls_data->z = NMS_checkpoint_0(data_NMS, n);
      data_NMS->ls_data->alpha0 = .95;
      data_NMS->ls_data->searchtype = LINESEARCH;
      alpha_watchdog = search_Armijo_standalone(n, &ref_merit, preRHS, data_NMS->ls_data);
      if (alpha_watchdog > data_NMS->alpha_min_watchdog)
      {
        /* successful line search */
        /* XXX this is not needed. we could get it from ls_data->zc */
        cblas_dcopy(n, NMS_checkpoint_0(data_NMS, n), 1, z_N, 1);
        cblas_daxpy(n, alpha_watchdog, workV1, 1, z_N, 1); //  z_N + tau*d --> z_N
        theta_iter = ref_merit; /* on output ref_merit contains the current value of theta*/
        DEBUG_PRINTF("NMS :: successful watchdog_step: alpha = %2.2e; theta_iter = %2.2e\n", alpha_watchdog, theta_iter);
        goto update_checkpoint;
      }
      /* else, we have to do a projected gradient step */
        DEBUG_PRINTF("NMS :: unsuccessful watchdog_step: alpha = %2.2e\n", alpha_watchdog);
      break;
    case ARCSEARCH:
      printf("watchdog_step: arc search not implemented yet !\n");
      exit(EXIT_FAILURE);
      //break;
    case BACKWARD_PATHSEARCH:
      printf("watchdog_step: path search not implemented yet !\n");
      exit(EXIT_FAILURE);
      //break;
    default:
      printf("watchdog_step: unknown search type : %d\n", data_NMS->watchdog_search_type);
      exit(EXIT_FAILURE);
  }

  /* compute the gradient JacThetaF_merit at the bestpoint z_b(0)*/
  functions->compute_F(data, NMS_bestpoint(data_NMS, n), F);
  functions->compute_H(data, NMS_bestpoint(data_NMS, n), F, workV1, workV2, Htmp);
  functions->compute_F_merit(data, NMS_bestpoint(data_NMS, n), F, F_merit);
  cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, H, n, F_merit, 1, 0.0, JacThetaF_merit, 1);

  switch(data_NMS->projected_gradient_search_type)
  {
    case LINESEARCH:
      /* workV1 is the same vector as the desc_dir used in the search */
      /* compute d = z_b(0) - JacThetaF_merit and project it onto the set*/
      cblas_dcopy(n, NMS_bestpoint(data_NMS, n), 1, workV1, 1);
      cblas_daxpy(n, -1.0, JacThetaF_merit, 1, workV1, 1);

      /* compute d_B = proj_B(d) */
      project_on_set(n, workV1, data_NMS->set);
      /* compute d_B - z_b(0) */
      cblas_daxpy(n, -1.0, NMS_bestpoint(data_NMS, n), 1, workV1, 1);
      /* compute <JacThetaF_merit, d_B - z_b(0)> */
      dotprod = cblas_ddot(n, JacThetaF_merit, 1, workV1, 1);

      /* check is reversed because we compute the <JacThetaF_merit, d_B - z_b(0)> */
      if (dotprod > 0.0) /* check condition on <JacThetaF_merit, z_b(0) - z_b(T_k)>*/
      {
//        preRHS = data_NMS->sigma*dotprod; /* we expect a plus in the LS function */
        preRHS = dotprod; /* we expect a plus in the LS function */
      }
      else
      {
//        preRHS = -data_NMS->sigma*data_NMS->merit_bestpoint; /* we expect a plus in the LS function */
        preRHS = -data_NMS->merit_bestpoint; /* we expect a plus in the LS function */
      }

      data_NMS->ls_data->z = NMS_bestpoint(data_NMS, n);
      data_NMS->ls_data->alpha0 = 1.0;
      data_NMS->ls_data->searchtype = LINESEARCH;
      alpha_projected_gradient = search_Armijo_standalone(n, &data_NMS->merit_bestpoint, preRHS, data_NMS->ls_data);
      if (alpha_projected_gradient > data_NMS->alpha_min_pgrad)
      {
        /* successful line search */
        theta_iter = data_NMS->merit_bestpoint; /* on output ref_merit contains the current value of theta*/
        DEBUG_PRINTF("NMS :: projected_gradient_step succeed: alpha = %2.2e; theta_iter = %2.2e\n", alpha_projected_gradient, theta_iter);
        /** \todo check sign here */
        cblas_dcopy(n, NMS_bestpoint(data_NMS, n), 1, z_N, 1);
        cblas_daxpy(n, alpha_projected_gradient, workV1, 1, z_N, 1); //  z_N + tau*d --> z_N
        goto update_checkpoint;
      }
      else
      {
        DEBUG_PRINTF("NMS :: projected_gradient_step unsuccessful: alpha = %2.2e\n", alpha_projected_gradient);
        printf("NMS: projected gradient unsuccessful ! Don't know what to do ... \n");
        return 1;
      }
    case ARCSEARCH:
      /* workV1 is the same vector as the desc_dir used in the search */
      /* desc_dir = -JacThetaF_merit */
      cblas_dcopy(n, JacThetaF_merit, 1, workV1, 1);
      cblas_dscal(n, -1.0, workV1, 1);
      /* preRHS is not useful here */
      data_NMS->ls_data->z = NMS_bestpoint(data_NMS, n);
      data_NMS->ls_data->alpha0 = 1.0;
      data_NMS->ls_data->searchtype = ARCSEARCH;
      alpha_projected_gradient = search_Armijo_standalone(n, &data_NMS->merit_bestpoint, preRHS, data_NMS->ls_data);
      if (alpha_projected_gradient > data_NMS->alpha_min_pgrad)
      {
        /* successful line search */
        theta_iter = data_NMS->merit_bestpoint; /* on output ref_merit contains the current value of theta*/
        DEBUG_PRINTF("NMS :: projected_gradient_step succeed: alpha = %2.2e; theta_iter = %2.2e\n", alpha_projected_gradient, theta_iter);
        /** \todo check sign here */
        cblas_dcopy(n, NMS_bestpoint(data_NMS, n), 1, z_N, 1);
        cblas_daxpy(n, alpha_projected_gradient, workV1, 1, z_N, 1); //  z_N + tau*d --> z_N
        goto update_checkpoint;
      }
      else
      {
        DEBUG_PRINTF("NMS :: projected_gradient_step unsuccessful: alpha = %2.2e\n", alpha_projected_gradient);
        printf("NMS: projected gradient unsuccessful ! Don't know what to do ... \n");
        return 1;
      }
    default:
      printf("projected gradient: unknown search type : %d\n", data_NMS->projected_gradient_search_type);
      exit(EXIT_FAILURE);
  }


/* This code is executed if there a m_step or a watchdog_step have been done*/
update_checkpoint:
  /* save the new checkpoint */
  DEBUG_PRINTF("NMS :: updating checkpoint : theta_iter = %e; ref_merit = %e\n", theta_iter, data_NMS->ref_merit);
  cblas_dcopy(n, z, 1, NMS_checkpoint_0(data_NMS, n), 1);
  cblas_dcopy(n, z_N, 1, NMS_checkpoint_T(data_NMS, n), 1);

  /* reinit some NMS  variables */
  data_NMS->n = 0;
  data_NMS->delta_checkpoint = data_NMS->delta;

  /* compute the new reference value associated with the merit function*/
  update_non_monotone_ref(data_NMS->ref_merit_data, theta_iter);
  data_NMS->ref_merit = theta_iter;
  get_non_monotone_ref(data_NMS->ref_merit_data, &data_NMS->ref_merit);

  /* save current checkpoint as bestpoint if there is a decrease in the merit
   * value */
  if (theta_iter <= data_NMS->merit_bestpoint)
  {
    /* save z_c(0) as bestpoint */
    /** \todo check if z_c(T) is really not needed */
    cblas_dcopy(n, NMS_checkpoint_0(data_NMS, n), 1, NMS_bestpoint(data_NMS, n), 1);
    /* update merit value and delta*/
    DEBUG_PRINTF("NMS :: updating bestpoint : theta_iter = %e; merit_bestpoint = %e\n", theta_iter, data_NMS->merit_bestpoint);
    data_NMS->delta_bestpoint = data_NMS->delta;
    data_NMS->merit_bestpoint = theta_iter;
  }

accept_z_N:
  DEBUG_PRINT("z_N z_0\n");
  DEBUG_EXPR_WE(for (unsigned i = 0; i < n; ++i)
      { DEBUG_PRINTF("%e %e\n", z_N[i], z[i]) });

  cblas_dcopy(n, z_N, 1, z, 1);

  return 0;
}


NMS_data* create_NMS_data(unsigned size, int matrix_type, int* restrict iparam, double* restrict dparam)
{
  NMS_data* data = (NMS_data*) malloc(sizeof(NMS_data));
  assert(data);
  data->size = size;
  data->n = 0;

  /* we have first to properly implement backward path search */
  data->path_data = NULL;

  data->checkpoint = (double*)malloc(2*size*sizeof(double));
  data->bestpoint = (double*)malloc(size*sizeof(double));
  data->workspace = (double*)malloc(5*size*sizeof(double));
  data->H = NM_create(matrix_type, size, size);
  data->ls_data = (search_data *)malloc(sizeof(search_data));

  /* set to NULL stuff that we can't create here */
  data->set = NULL;
  data->path_data = NULL;

  switch (iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS])
  {
    case NM_LS_MAX:
    case NM_LS_MEAN:
      data->ref_merit_data = malloc(sizeof(nm_ref_struct));
      fill_nm_data((nm_ref_struct*)data->ref_merit_data, iparam);
      break;
    default:
      printf("create_NMS_data :: unknown search type %d\n", iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS]);
      exit(EXIT_FAILURE);
  }

  /* set some parameters for the NMS; default values should be set */
  data->watchdog_search_type = iparam[SICONOS_IPARAM_NMS_WATCHDOG_TYPE];
  data->projected_gradient_search_type = iparam[SICONOS_IPARAM_NMS_PROJECTED_GRADIENT_TYPE];
  data->n_max = iparam[SICONOS_IPARAM_NMS_N_MAX];

  data->delta = dparam[SICONOS_DPARAM_NMS_DELTA];
  data->delta_var = dparam[SICONOS_DPARAM_NMS_DELTA_VAR];
  data->sigma = dparam[SICONOS_DPARAM_NMS_SIGMA];
  data->alpha_min_watchdog = dparam[SICONOS_DPARAM_NMS_ALPHA_MIN_WATCHDOG];
  data->alpha_min_pgrad = dparam[SICONOS_DPARAM_NMS_ALPHA_MIN_PGRAD];
  data->merit_incr = dparam[SICONOS_DPARAM_NMS_MERIT_INCR];

  /* init some values */
  data->delta_checkpoint = data->delta;
  data->delta_bestpoint = data->delta;

  return data;
}

void free_NMS_data(NMS_data* data)
{
  assert(data);
  if (data->ref_merit_data)
  {
    /* XXX fix this */
    /** \todo sort the nm data thing */
    free_nm_data((nm_ref_struct *)data->ref_merit_data);
    free(data->ref_merit_data);
  }

  free(data->checkpoint);
  free(data->bestpoint);
  free(data->workspace);
  NM_free(data->H);
  free(data->H);
  free_siconos_set(data->set);
  free(data->set);
  free_ls_data(data->ls_data);
  free(data->ls_data);

  if (data->path_data)
    free(data->path_data);

  free(data);
}
