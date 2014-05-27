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

//#include "SolverOptions.h"
//#include "MixedComplementarityProblem.h"
#include <assert.h>

#include "MCP_Solvers.h"
#include "MCP_cst.h"
#include "MCP_FischerBurmeister.h"

char SICONOS_MCP_FB_STR[] = "NewtonFB";
char SICONOS_MCP_NEWTON_FBLSA_STR[] = "Newton FBLSA";
char SICONOS_MCP_NEWTON_MINFBLSA_STR[] = "Newton minFBLSA";

int mcp_driver2(MixedComplementarityProblem2* problem, double *z , double *Fmcp, SolverOptions* options, NumericsOptions* global_options)
{
  assert(options != NULL);

  /* Set global options */
  if (global_options)
    setNumericsOptions(global_options);

  /* Checks inputs */
  assert(problem != NULL);
  assert(z != NULL);
  assert(Fmcp != NULL);
  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  switch (options->solverId)
  {
  case SICONOS_MCP_NEWTON_FBLSA: // Fischer-Burmeister/Newton -- new version
    mcp_newton_FBLSA(problem, z, Fmcp, &info, options);
    break;

  case SICONOS_MCP_NEWTON_MINFBLSA: // Fischer-Burmeister/Newton + min descent direction
    mcp_newton_minFBLSA(problem, z, Fmcp, &info, options);
    break;

  default:
    fprintf(stderr, "mcp_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

  return info;
}

int mcp_driver(MixedComplementarityProblem* problem, double *z , double *w, SolverOptions* options,  NumericsOptions* global_options)
{
  if (options == NULL)
    numericsError("mcp_driver ", "null input for solver options.\n");

  /* Set global options */
  if (global_options)
    setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("mcp_driver", "null input for MixedComplementarityProblem and/or unknowns (z,w)");
  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  switch (options->solverId)
  {
  case SICONOS_MCP_FB: // Fischer-Burmeister/Newton
    mcp_FischerBurmeister(problem, z, w, &info, options);
    break;

  default:
    fprintf(stderr, "mcp_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);

  }

  return info;
}

void mcp_driver_init(MixedComplementarityProblem* problem, SolverOptions* options)
{
  switch (options->solverId)
  {
  case SICONOS_MCP_FB :
    mcp_FischerBurmeister_init(problem, options) ;
    break ;
  default :
    fprintf(stderr, "mcp_driver_init error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

}

void mcp_driver_reset(MixedComplementarityProblem* problem, SolverOptions* options)
{
  switch (options->solverId)
  {
  case SICONOS_MCP_FB :
    mcp_FischerBurmeister_reset(problem, options) ;
    break ;
  default :
    fprintf(stderr, "mcp_driver_init error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

}
