/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FrictionContactProblem.h"
#include "Friction_cst.h"
#include "MohrCoulomb2DProblem.h"
#include "NumericsMatrix.h"
#include "NonSmoothDrivers.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#include "frictionContact_test_utils.h"
#include "numerics_verbose.h"


int main(void) {
  int total_info = 0;

  //  numerics_set_verbose(1);

  double q[] = {-1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[] = {0.1, 0.1, 0.1};
  double eta[] = {0.1, 0.1, 0.1};

  double Wdata[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  NumericsMatrix* tmpM = NM_create_from_data(NM_DENSE, 9, 9, Wdata);
  NumericsMatrix* W = NM_create(NM_SPARSE, 9, 9);
  NM_copy_to_sparse(tmpM, W, 1e-16);


  MohrCoulomb2DProblem* MC2D = mohrCoulomb2DProblem_new_with_data(3, 3, W, q, eta, mu);
  double stress[9] = {0.};
  double plastic_strain_rate[9] = {0.};
  mohrCoulomb2D_display(MC2D);
  mohrCoulomb2D_printInFilename(MC2D, "mc2d_example1.dat");
  
  MC2D->M = NULL;
  MC2D->q = NULL;
  MC2D->eta = NULL;
  MC2D->mu = NULL;
  mohrCoulomb2DProblem_free(MC2D);
  
  MohrCoulomb2DProblem* MC2D_r= mohrCoulomb2D_new_from_filename("mc2d_example1.dat");
  mohrCoulomb2D_display(MC2D_r);
  mohrCoulomb2DProblem_free(MC2D_r);
  
  FrictionContactProblem* FC = frictionContactProblem_new_with_data(3, 3, W, q, mu);
  double r[9] = {0.};
  double u[9] = {0.};

  FC->M = NULL;
  FC->q = NULL;
  FC->mu = NULL;
  frictionContactProblem_free(FC);


  FrictionContactProblem* FCdense = frictionContactProblem_new_with_data(3, 3, tmpM, q, mu);
  double rdense[9] = {0.};
  double udense[9] = {0.};
    
  FCdense->M = NULL;
  FCdense->q = NULL;
  FCdense->mu = NULL;
  frictionContactProblem_free(FCdense);
  
  NM_clear(W);
  tmpM->matrix0 = NULL;
  NM_clear(tmpM);
  free(W);
  free(tmpM);

  return total_info;
}
