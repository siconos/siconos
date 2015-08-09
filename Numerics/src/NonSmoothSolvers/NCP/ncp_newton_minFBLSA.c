/* Siconos-Numerics, Copyright INRIA 2005-2014
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


#include <stdio.h>
#include <math.h>
#include <float.h>

#include "NCP_Solvers.h"
#include "SiconosLapack.h"
#include "Newton_Methods.h"
#include "FischerBurmeister.h"
#include "min_merit.h"
#include "ncp_newton_FBLSA.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

static void ncp_min(void* data_opaque, double* z, double* F, double* Fmin)
{
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem *)data_opaque;

  F_min(0, data->n, z, F, Fmin);
}

static void min_compute_H_ncp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, NumericsMatrix* H)
{
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem *)data_opaque;

  data->compute_nabla_F(data->env, data->n, z, data->nabla_F);

  Jac_F_min(0, data->n, z, F, data->nabla_F, H);
}

void ncp_newton_minFBLSA(NonlinearComplementarityProblem* problem, double *z, double* F, int *info , SolverOptions* options)
{
  functions_LSA functions_minFBLSA_ncp;
  init_lsa_functions(&functions_minFBLSA_ncp, &FB_compute_F_ncp, &ncp_FB);
  functions_minFBLSA_ncp.compute_H = &FB_compute_H_ncp;
  functions_minFBLSA_ncp.compute_error = &FB_compute_error_ncp;
  functions_minFBLSA_ncp.compute_RHS_desc = &ncp_min;
  functions_minFBLSA_ncp.compute_H_desc = &min_compute_H_ncp;

  set_lsa_params_data(options, problem->nabla_F);
  newton_LSA(problem->n, z, F, info, (void *)problem, options, &functions_minFBLSA_ncp);
}
