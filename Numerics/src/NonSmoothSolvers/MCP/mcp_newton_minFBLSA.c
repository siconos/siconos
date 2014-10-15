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


#include <stdio.h>
#include <math.h>
#include <float.h>

#include "MCP_Solvers.h"
#include "SiconosLapack.h"
#include "Newton_Methods.h"
#include "FischerBurmeister.h"
#include "min_merit.h"
#include "mcp_newton_FBLSA.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void mcp_min(void* data_opaque, double* z, double* F, double* Fmin)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;

  F_min(data->n1, data->n2, z, F, Fmin);
}

void min_compute_H_mcp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, double* H)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;
  assert(data->nabla_Fmcp);

  data->compute_nabla_Fmcp(data->env, data->n1, data->n2, z, data->nabla_Fmcp);

  Jac_F_min(data->n1, data->n2, z, F, data->nabla_Fmcp, H);
}

void mcp_newton_minFBLSA(MixedComplementarityProblem2* problem, double *z, double* Fmcp, int *info , SolverOptions* options)
{
  functions_LSA functions_minFBLSA_mcp;
  functions_minFBLSA_mcp.compute_F = &FB_compute_F_mcp;
  functions_minFBLSA_mcp.compute_F_merit = &mcp_FB;
  functions_minFBLSA_mcp.compute_H = &FB_compute_H_mcp;
  functions_minFBLSA_mcp.compute_error = &FB_compute_error_mcp;
  functions_minFBLSA_mcp.compute_RHS_desc = &mcp_min;
  functions_minFBLSA_mcp.compute_H_desc = &min_compute_H_mcp;

 newton_LSA(problem->n1 + problem->n2, z, Fmcp, info, (void *)problem, options, &functions_minFBLSA_mcp);
}
