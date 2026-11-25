/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "SolverOptions.h"
#include "numerics_verbose.h"

int main(void) {
  int total_info = 0;

  //  numerics_set_verbose(1);

  double q[] = {-1, 1, 3, -1, 1, 3, -1, 1, 3};
  double theta[] = {0.1, 0.1, 0.1};
  double eta[] = {0.1, 0.1, 0.1};

  double Wdata[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  NumericsMatrix* tmpM = NM_create_with_data(NM_DENSE, 9, 9, Wdata);
  NumericsMatrix* W = NM_copy_to_sparse(tmpM, 1e-16);

  PlasticityProblem* PLASTICITY_2D = plasticity2DProblem_new_with_data(3, 3, W, q, eta, theta);
  double stress[9] = {0.};
  double plastic_strain_rate[9] = {0.};
  plasticity2D_display(PLASTICITY_2D);
  plasticity2D_printInFilename(PLASTICITY_2D, "plasticity_2d_example1.dat");

  PLASTICITY_2D->M = NULL;
  PLASTICITY_2D->q = NULL;
  PLASTICITY_2D->model.drucker_prager->eta = NULL;
  PLASTICITY_2D->model.drucker_prager->theta = NULL;
  plasticity2DProblem_free(PLASTICITY_2D);

  PlasticityProblem* PLASTICITY_2D_r =
      plasticity2D_new_from_filename("plasticity_2d_example1.dat");
  plasticity2D_display(PLASTICITY_2D_r);
  plasticity2DProblem_free(PLASTICITY_2D_r);

  NM_clear(W);
  tmpM->dense = NULL;
  NM_clear(tmpM);
  free(W);
  free(tmpM);

  return total_info;
}
