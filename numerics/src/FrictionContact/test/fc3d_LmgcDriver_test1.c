/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdlib.h>            // for free, malloc
#include "Friction_cst.h"      // for SICONOS_FRICTION_3D_NSGS
#include "NonSmoothDrivers.h"  // for fc3d_LmgcDriver

int main(void)
{
  int info = 0 ;
  int nc = 4;//Number of contacts
  int nb = 4;//Number of blocks

  double q[12] = { -1, 1, 3, -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[4] = {0.1, 0.1, 0.1, 0.1};

  unsigned int row[4] = {1, 2, 3, 4};
  unsigned int column[4] = {1, 2, 3, 4};
  double W[36] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1};


  double *reaction = (double*)malloc(3 * nc * sizeof(double));
  double *velocity = (double*)malloc(3 * nc * sizeof(double));
  for (int i = 0; i < 3 * nc; i++)
  {
    reaction[i] = 0.0;
    velocity[i] = 0.0;
  }

  int solver_id = SICONOS_FRICTION_3D_NSGS; // 500
  double tolerance = 1e-16;
  int itermax = 100;

  info = fc3d_LmgcDriver(reaction,
                                      velocity,
                                      q,
                                      mu,
                                      W,
                                      row,
                                      column,
                                      nc,
                                      nb,
                                      solver_id,
                                      tolerance,
                                      itermax,
                                      0, 0, 0, 0);

  free(reaction);
  free(velocity);

  return info;
}
