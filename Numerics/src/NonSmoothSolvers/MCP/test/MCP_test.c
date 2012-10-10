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

#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "mcp_cst.h"



void testF(int size, double *z, double * F)
{
  printf("call to MCP function F(z) ...\n");
}

void testNablaF(int size, double *z, double *F)
{
  printf("call to MCP function nablaF(z) ...\n");
}

int main(void)
{
  printf(" Start tests for MCP solvers.\n");

  int info = 0 ;

  /* Set solver options */
  SolverOptions options;

  /* FB solver */
  options.solverId = SICONOS_MCP_FB;

  /* Create a MixedComplementarityProblem */
  MixedComplementarityProblem* problem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));

  problem->sizeEqualities = 2;
  problem->sizeInequalities = 3;
  problem->computeFmcp = &testF ;
  problem->computeNablaFmcp = &testNablaF ;

  int size = 5;
  double z[4];
  double F[4];
  double nablaF[16];
  problem->computeFmcp(size, z, F);
  problem->computeNablaFmcp(size, z, nablaF);

  /* Initialize the solver */
  mcp_driver_init(problem, &options) ;

  /// TODO : write a real test ... ////

  printf("End of MCP solvers test. \n");

  return info;
}
