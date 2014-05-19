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
#include "MCP_cst.h"
#include "SiconosLapack.h"
#include "Newton_Methods.h"
#include "FischerBurmeister.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void FB_compute_F_mcp(void* data_opaque, double* z, double* Fmcp)
{
  // Computation of the new value F(z)
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;
  data->compute_Fmcp(data->env_compute_Fmcp, data->n1, data->n2, z, Fmcp);
}

void FB_compute_H_mcp(void* data_opaque, double* z, double* Fmcp, double* workV1, double* workV2, double* H)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;

  data->compute_nabla_Fmcp(data->env_compute_nabla_Fmcp, data->n1, data->n2, z, data->nabla_Fmcp);

  Jac_F_FB(data->n1, data->n2, z, Fmcp, workV1, workV2, data->nabla_Fmcp, H);
}

void FB_compute_error_mcp(void* data_opaque, double* z, double* w, double* Jac_F_merit, double tol, double* err)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;
  unsigned int n = data->n1 + data->n2;
  err[0] = cblas_dnrm2(n, Jac_F_merit, 1);
}

void mcp_FB(void* data_opaque, double* z, double* F, double* F_FB)
{
  MixedComplementarityProblem2* data = (MixedComplementarityProblem2 *)data_opaque;
  phi_Mixed_FB(data->n1, data->n2, z, F, F_FB);
}

void mcp_newton_FBLSA(MixedComplementarityProblem2* problem, double *z, double* Fmcp, int *info , SolverOptions* options)
{
  functions_FBLSA functions_FBLSA_lcp;
  functions_FBLSA_lcp.compute_F = &FB_compute_F_mcp;
  functions_FBLSA_lcp.compute_F_merit = &mcp_FB;
  functions_FBLSA_lcp.compute_H = &FB_compute_H_mcp;
  functions_FBLSA_lcp.compute_error = &FB_compute_error_mcp;
  functions_FBLSA_lcp.compute_H_desc = NULL;
  functions_FBLSA_lcp.compute_F_desc = NULL;

  newton_FBLSA(problem->n1 + problem->n2, z, Fmcp, info, (void *)problem, options, &functions_FBLSA_lcp);
}

int mixedComplementarity_newton_FBLSA_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the Newton based FBLSA MCP Solver\n");
  }

  options->solverId = SICONOS_MCP_NEWTON_FBLSA;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;

  options->iparam[0] = 1000;
  options->dparam[0] = 1e-10;

  return 0;
}

