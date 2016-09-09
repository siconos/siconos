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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"


int main(void)
{
  int total_info = 0;

  double q[] = { -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[] = {0.1, 0.1, 0.1};

  double Wdata[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  NumericsMatrix* tmpM = createNumericsMatrixFromData(NM_DENSE, 9, 9, Wdata);
  NumericsMatrix* W = createNumericsMatrix(NM_SPARSE, 9, 9);
  NM_copy_to_sparse(tmpM, W);

  int solvers_to_test[] = {SICONOS_FRICTION_3D_NSGS};
  /* Note computeSparseAWpB only for sparse block storage at the moment */
  /*                         SICONOS_FRICTION_3D_NSN_AC,
                             SICONOS_FRICTION_3D_NSN_FB,
                             SICONOS_FRICTION_3D_NSN_NM,
                             SICONOS_FRICTION_3D_SOCLCP,
                             SICONOS_FRICTION_3D_PROX};*/

  for (size_t s = 0; s < 1; ++s)
  {
    int solver_id = solvers_to_test[s];

    FrictionContactProblem* FC = frictionContactProblem_new(3, 3, W, q, mu);
    double r[9] = {0.};
    double u[9] = {0.};

    SolverOptions SO;;
    fc3d_setDefaultSolverOptions(&SO, solver_id);
    int info = fc3d_driver(FC, r, u, &SO);

    if (info)
    {
      fprintf(stderr, "Solver %s failed with error %d\n", solver_options_id_to_name(solver_id), info);
      total_info = 1;
    }
    FC->M = NULL;
    FC->q = NULL;
    FC->mu = NULL;
    solver_options_delete(&SO);
    freeFrictionContactProblem(FC);
  }

  freeNumericsMatrix(W);
  tmpM->matrix0 = NULL;
  freeNumericsMatrix(tmpM);
  free(W);
  free(tmpM);

  return total_info;
}
