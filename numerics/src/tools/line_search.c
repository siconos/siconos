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

#include <assert.h>
#include "line_search.h"
#include <math.h>
#include "SiconosCompat.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
using namespace std;
#endif

double line_search_generic(int n, double theta, double preRHS, search_data* ls_data, unsigned searchtype, sn_ls_fn ls_fn)
{
  double theta_ref = theta;

  if (ls_data->nm_ref_data)
    get_non_monotone_ref(ls_data->nm_ref_data, &theta_ref);

  ls_data->searchtype = searchtype;
  assert(searchtype != ARCSEARCH || ls_data->set);
  double alpha = (*ls_fn)(n, &theta_ref, preRHS, ls_data);
 
  if (ls_data->nm_ref_data)
  {
    if (isfinite(alpha)) update_non_monotone_ref(ls_data->nm_ref_data, theta_ref);
    else zero_nm_data((nm_ref_struct*)ls_data->nm_ref_data);
  }

  return alpha;
}

void update_non_monotone_ref(void* nm_ref_data, double cur_merit)
{
  assert(nm_ref_data);
  nm_ref_struct* data = (nm_ref_struct*) nm_ref_data;
  if (data->m < data->M)
  {
    data->previous_thetas[data->m] = cur_merit;
    data->m++;
  }
  else if (data->M > 0)
  {
    for (int i = 0; i < data->M-1; ++i) data->previous_thetas[i] = data->previous_thetas[i+1];
    data->previous_thetas[data->M-1] = cur_merit;
  }
}

void get_non_monotone_ref(void* nm_ref_data, double* theta_ref)
{
  assert(nm_ref_data);
  assert(theta_ref);
  double local_theta_ref = 0.;
  nm_ref_struct* data_max;
  nm_ref_struct* data_mean;

  switch (get_nonmonotone_type(nm_ref_data))
  {
    case NM_LS_MAX: // classical nonmonotone theta_ref = max theta_j
      data_max = (nm_ref_struct*) nm_ref_data;
      local_theta_ref = *theta_ref;
      for (size_t i = 0; i < data_max->m; ++i)
      {
        if (data_max->previous_thetas[i] > local_theta_ref)
        {
          local_theta_ref = data_max->previous_thetas[i];
        }
      }
      *theta_ref = local_theta_ref;
      break;

    case NM_LS_MEAN: // mean like value : theta_ref = max { theta, mean(theta) }
      data_mean = (nm_ref_struct*)nm_ref_data;
      if (data_mean->m > 0)
      {
        for (size_t i = 0; i < data_mean->m; ++i)
        {
          local_theta_ref += data_mean->previous_thetas[i];
        }
        local_theta_ref /= (double)(data_mean->m);
        if (local_theta_ref < *theta_ref)
        {
          *theta_ref = local_theta_ref;
        }
      }
      break;

    case NM_LS_DISABLE:
      break;
    case NM_LS_ZHANG_HAGER:
    default:
      printf("update_non_monotone_ref :: not implemented for type %d\n",
          get_nonmonotone_type(nm_ref_data));
      assert(0 && "unknown nonmonotone strategy!");
      exit(EXIT_FAILURE);
  }
}


void fill_nm_data(nm_ref_struct* nm_ref_data, int* restrict iparam)
{
  assert(nm_ref_data);
  assert(iparam);

  nm_ref_data->type = iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS];
  nm_ref_data->M = iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M];
  nm_ref_data->m = 0;
  if (nm_ref_data->M > 0)
    nm_ref_data->previous_thetas = (double*)calloc(nm_ref_data->M, sizeof(double));
}

void zero_nm_data(nm_ref_struct* nm_ref_data)
{
  assert(nm_ref_data);

  nm_ref_data->m = 0;
}

void free_nm_data(nm_ref_struct* nm_ref_data)
{
  assert(nm_ref_data);
  if (nm_ref_data->M > 0)
  {
    free(nm_ref_data->previous_thetas);
    nm_ref_data->previous_thetas = NULL;
  }
}

void free_ls_data(search_data* ls_data)
{
  if (ls_data->nm_ref_data)
  {
    free_nm_data((nm_ref_struct*)ls_data->nm_ref_data);
  }
  if (ls_data->extra_params)
  {
    free(ls_data->extra_params);
    ls_data->extra_params = NULL;
  }
}
