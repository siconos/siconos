/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#include "frictionContact_test_function.h"


int main(void)
{
  NumericsOptions NO;
  setDefaultNumericsOptions(&NO);
  NO.verboseMode = 1; // turn verbose mode to off by default

  int total_info = 0;

  double q[] = { -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[] = {0.1, 0.1, 0.1};

  double Wdata[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  NumericsMatrix* tmpM = createNumericsMatrixFromData(NM_DENSE, 9, 9, Wdata);
  NumericsMatrix* W = createNumericsMatrix(NM_SPARSE, 9, 9);
  NM_copy_to_sparse(tmpM, W);

  int solvers_to_test[] = {SICONOS_FRICTION_3D_NSGS,
                           SICONOS_FRICTION_3D_NSN_AC,
                           SICONOS_FRICTION_3D_NSN_FB,
                           SICONOS_FRICTION_3D_NSN_NM,
                           SICONOS_FRICTION_3D_SOCLCP,
                           SICONOS_FRICTION_3D_PROX};

  for (size_t s = 0; s < sizeof(solvers_to_test); ++s)
  {
    int solver_id = solvers_to_test[s];

    FrictionContactProblem* FC = frictionContactProblem_new(3, 3, W, q, mu);
    double r[9] = {0.};
    double u[9] = {0.};

    SolverOptions SO;;
    fc3d_setDefaultSolverOptions(&SO, solver_id);
    int info = fc3d_driver(FC, r, u, &SO, &NO);

    if (info)
    {
      fprintf(stderr, "Solver %s failed with error %d\n", idToName(solver_id), info);
      total_info = 1;
    }
    FC->M = NULL;
    FC->q = NULL;
    FC->mu = NULL;
    deleteSolverOptions(&SO);
    freeFrictionContactProblem(FC);
    free(FC);
  }

  freeNumericsMatrix(W);
  tmpM->matrix0 = NULL;
  freeNumericsMatrix(tmpM);
  free(W);
  free(tmpM);

  return total_info;
}
