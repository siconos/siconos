/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include "LA.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"
#include "NonSmoothNewton.h"
#include "FischerBurmeister.h"
#include "MCP_solvers.h"


/* Static object which contains the MCP problem description.
Ugly but required to deal with function pointer connection
in  FischerFunc_MCP and its jacobian.
*/
static MixedComplementarityProblem * localProblem = 0;


void mcp_FischerBurmeister_init(MixedComplementarityProblem * problem, SolverOptions* options)
{
  localProblem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));
  /* Connect local static problem with the "real" MCP */
  localProblem->sizeEqualities = problem->sizeEqualities ;
  localProblem->sizeInequalities = problem->sizeInequalities ;
  localProblem->computeFmcp = problem->computeFmcp ;
  localProblem->computeNablaFmcp = problem->computeNablaFmcp ;

  int fullSize = localProblem->sizeEqualities + localProblem->sizeInequalities ;
  // Memory allocation for working vectors
  int lwork = 3 * (problem->sizeEqualities + problem->sizeInequalities) ;
  options->dWork = (double *) malloc(lwork * sizeof(double));
  localProblem->nablaFmcp = &(options->dWork[fullSize]) ;
  localProblem->Fmcp = options->dWork ;
}

void mcp_FischerBurmeister_reset(MixedComplementarityProblem * problem, SolverOptions* options)
{
  localProblem->Fmcp = NULL;
  localProblem->nablaFmcp = NULL;
  freeMixedComplementarityProblem(localProblem);
}

// Must corresponds to a NewtonFunctionPtr
void FischerFunc_MCP(int size, double* z, double* phi, int dummy)
{
  // This function uses a user-defined function, set in the problem, to compute
  // the Fisher function

  int sizeEq = localProblem->sizeEqualities;
  int sizeIneq = localProblem->sizeInequalities;
  /* First call user-defined function to compute Fmcp function, */
  localProblem->computeFmcp(sizeEq + sizeIneq, z, localProblem->Fmcp) ;
  /* and compute the corresponding Fischer function */
  phi_Mixed_FB(sizeEq, sizeIneq, z, localProblem->Fmcp, phi) ;
}

// Must corresponds to a NewtonFunctionPtr
void nablaFischerFunc_MCP(int size, double* z, double* nablaPhi, int dummy)
{
  int sizeEq = localProblem->sizeEqualities;
  int sizeIneq = localProblem->sizeInequalities;
  /* First call user-defined function to compute Fmcp function, */
  localProblem->computeNablaFmcp(sizeEq + sizeIneq, z, localProblem->nablaFmcp) ;
  /* and compute the corresponding jacobian of the Fischer function */
  jacobianPhi_Mixed_FB(sizeEq, sizeIneq, z, localProblem->Fmcp, localProblem->nablaFmcp, nablaPhi) ;
}

void mcp_FischerBurmeister(MixedComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  *info = 1;
  int fullSize = problem->sizeEqualities + problem->sizeInequalities ;

  // Set links to Fisher functions and its jacobian
  NewtonFunctionPtr phi = &FischerFunc_MCP ;
  NewtonFunctionPtr nablaPhi = &nablaFischerFunc_MCP ;

  // Call semi-smooth Newton solver
  *info = nonSmoothNewton(fullSize, z, &phi, &nablaPhi, options->iparam, options->dparam);

  // todo : compute error function

  // Check output
  if (*info > 0)
    fprintf(stderr, "Numerics, mcp_FB failed, reached max. number of iterations without convergence. Error = %f\n", options->dparam[1]);

  return;
}

int mixedComplementarity_FB_setDefaultSolverOptions(MixedComplementarityProblem* problem, SolverOptions* pSolver)
{
  mixedComplementarity_default_setDefaultSolverOptions(problem, pSolver);
  return 0;
}


