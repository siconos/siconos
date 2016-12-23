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


//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"
#include "float.h"
#include "numerics_verbose.h"
#include <assert.h>

#include "SiconosBlas.h"
#include "ArmijoSearch.h"
#include "SiconosSets.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

double search_Armijo_standalone(int n, double* theta, double preRHS, search_data* ls_data)
{
  assert(ls_data->alpha0 > 0.0);
  assert(ls_data->alpha0 > ls_data->alpha_min);
  double alpha = ls_data->alpha0;
  double theta_iter = *theta, theta_ref = *theta;
  double* z = ls_data->z;
  double* zc = ls_data->zc;
  double* F = ls_data->F;
  double* F_merit = ls_data->F_merit;
  double* desc_dir = ls_data->desc_dir;
  void* data = ls_data->data;
  bool arcsearch = ls_data->searchtype == ARCSEARCH;
  void* set = ls_data->set;
  double RHS;

  armijo_extra_params* aep = (armijo_extra_params*) ls_data->extra_params;
  assert(aep);
  preRHS *= aep->gamma;

  while (alpha >= ls_data->alpha_min)
  {
     // desc_dir contains the direction d
     cblas_dcopy(n, z, 1, zc, 1);
     cblas_daxpy(n, alpha, desc_dir, 1, zc, 1);     //  z + alpha*d --> z
     if (arcsearch)
     {
       project_on_set(n, zc, set);
       /* we use F as a work vector here */
       cblas_dcopy(n, z, 1, F, 1);
       cblas_daxpy(n, -1.0, zc, 1, F, 1); /* F = z(0) - z(alpha) */
       /* warning desc_dir = -JacMerit !*/
       double dotprod = cblas_ddot(n, desc_dir, 1, F, 1);
       if (dotprod > 0.0)
         RHS = ls_data->sigma*dotprod;
       else
         RHS = -alpha*ls_data->sigma*theta_ref;
     }
     else
     {
       RHS = alpha*preRHS;
     }

     // compute new F_merit
     ls_data->compute_F(data, zc, F);
     ls_data->compute_F_merit(data, zc, F, F_merit);

     DEBUG_PRINT("z ");
     DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
         { DEBUG_PRINTF("% 2.2e ", zc[i]) }
         DEBUG_PRINT("\n"));
 
     DEBUG_PRINT("F ");
     DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
         { DEBUG_PRINTF("% 2.2e ", F[i]) }
         DEBUG_PRINT("\n"));
 
     DEBUG_PRINT("F_merit ");
     DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
         { DEBUG_PRINTF("% 2.2e ", F_merit[i]) }
         DEBUG_PRINT("\n"));
 
     theta_iter = 0.5 * cblas_ddot(n, F_merit, 1, F_merit, 1);
 
     DEBUG_PRINTF("search_Armijo :: alpha %g\n", alpha);
     DEBUG_PRINTF("search_Armijo :: theta_iter %.*e ; theta_ref %.*e  \n", DECIMAL_DIG, theta_iter, DECIMAL_DIG, theta_ref);
 
     // acceptance test
     if (theta_iter <= theta_ref + RHS)
     {
       if (verbose > 1)
         printf("search_Armijo :: alpha %g\n", alpha);
       break;
     }
     else
     {
       // alpha too large, decrease it
       alpha /= 2.0;
     }
  }
  *theta = theta_iter;
  return alpha;
}

void search_Armijo_params_init(armijo_extra_params* p)
{
  p->gamma = 1e-4;
}
