/* Siconos-Numerics, Copyright INRIA 2005-2014.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */


//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"
#include "SiconosBlas.h"
#include "LineSearch.h"

double linesearch_Armijo2(int n, double theta, double preRHS, linesearch_data* ls_data)
{
  double alpha = 2.0;
  double theta_iter;
  double* z = ls_data->z;
  double* zc = ls_data-> zc;
  double* F = ls_data->F;
  double* F_merit = ls_data->F_merit;
  double* desc_dir = ls_data->desc_dir;
  void* data = ls_data->data;

  double theta_ref = theta;

  switch (ls_data->nonmonotone)
  {
    case 1: // classical nonmonotone theta_ref = max theta_j
      for (int i = 0; i < ls_data->m; ++i)
      {
        if (ls_data->previous_theta[i] > theta_ref)
        {
          theta_ref = ls_data->previous_theta[i];
        }
      }
      break;

    case 2: // mean like value : theta_ref = max { theta, mean(theta) }
      for (int i = 0; i < ls_data->m; ++i)
      {
        theta_ref += ls_data->previous_theta[i];
      }
      theta_ref /= (double)(ls_data->m+1);
      if (theta_ref < theta)
      {
        theta_ref = theta;
      }
      break;

    default:
      ;
  }

  while (1)
  {
     // workV1 contains the direction d
     cblas_dcopy(n, z, 1, zc, 1);
     cblas_daxpy(n, alpha, desc_dir, 1, zc, 1);     //  z + alpha*d --> z
 
     // compute new F_merit
     ls_data->compute_F(data, zc, F);
     ls_data->compute_F_merit(data, zc, F, F_merit);
 
     DEBUG_PRINT("z ");
     DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
         { DEBUG_PRINTF("% 2.2e ", workV[i]) }
         DEBUG_PRINT("\n"));
 
     DEBUG_PRINT("F ");
     DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
         { DEBUG_PRINTF("% 2.2e ", F[i]) }
         DEBUG_PRINT("\n"));
 
     DEBUG_PRINT("F_merit ");
     DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
         { DEBUG_PRINTF("% 2.2e ", F_merit[i]) }
         DEBUG_PRINT("\n"));
 
     theta_iter = cblas_dnrm2(n, F_merit, 1);
     theta_iter = 0.5 * theta_iter * theta_iter;
 
     DEBUG_PRINTF("newton_FBLSA :: alpha %g\n", alpha);
     DEBUG_PRINTF("newton_FBLSA :: theta_iter %.*e ; theta %.*e  \n", DECIMAL_DIG, theta_iter, DECIMAL_DIG, theta);
 
     // acceptance test
     if (theta_iter <= theta_ref + alpha*preRHS)
     {
       break;
       if (verbose > 1)
         printf("newton_FBLSA :: alpha %g\n", alpha);
     }
     else
     {
       // alpha too large, decrease it
       alpha /= 2.0;
     }
  }

  if (ls_data->m < ls_data->M)
  {
    ls_data->previous_theta[ls_data->m] = theta_iter;
    ls_data->m++;
  }
  else if (ls_data->M > 0)
  {
    for (int i = 0; i < ls_data->M-1; ++i) ls_data->previous_theta[i] = ls_data->previous_theta[i+1];
    ls_data->previous_theta[ls_data->m] = theta_iter;
  }
  return alpha;
}


