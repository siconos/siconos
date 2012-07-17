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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "mcp_cst.h"

using namespace std;

void testF(double *z)
{
  cout << "call to MCP function F(z) ..." << endl;
}

void testNablaF(double *z)
{
  cout << "call to MCP function nablaF(z) ..." << endl;
}

int main(void)
{
  cout << " Start tests for MCP solvers." << endl;

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


  double z[4];
  problem->computeFmcp(z);
  problem->computeNablaFmcp(z);

  /* Initialize the solver */
  mcp_driver_init(problem, &options) ;

  /// TODO : write a real test ... ////

  cout << "End of MCP solvers test. " << endl;

  return info;
}
