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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"
#include "numerics_verbose.h"
#include "Friction_cst.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  int info = 0 ;


  int nc = 4 ;
  int nb = 16 ;
  double mu[4] =
  {
    7.500000e-01,   7.500000e-01,   7.500000e-01,   7.500000e-01
  };
  unsigned int row[16] =
  {
    1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,   4
  };
  unsigned int column[16] =
  {
    1,  2,  3,  4,  1,  2,  3,  4,  1,  2,  3,  4,  1,  2,  3,   4
  };
  double q[12] =
  {
    -9.810000e-03,  0.000000e+00, -6.006893e-19,  -9.810000e-03,  0.000000e+00, -6.006893e-19,  -9.810000e-03,  0.000000e+00, -6.006893e-19,  -9.810000e-03,  0.000000e+00,  -6.006893e-19
  };
  double W[144] =
  {
    9.906380e+00,   -2.190000e+00,  -2.166000e+00,  -2.190000e+00,  5.643519e+00,   -2.470594e+00,  -2.166000e+00,  -2.470594e+00,  5.697969e+00,
    2.087120e+00,   -2.190000e+00,  2.166000e+00,   -2.190000e+00,  7.564814e-01,   2.470594e+00,   -2.166000e+00,  -2.470594e+00,  5.697969e+00,
    -5.818780e+00,  2.142000e+00,   2.166000e+00,   -2.190000e+00,  7.564814e-01,   2.470594e+00,   -2.166000e+00,  2.416444e+00,   7.567814e-01,
    2.000480e+00,   2.142000e+00,   -2.166000e+00,  -2.190000e+00,  5.643519e+00,   -2.470594e+00,  -2.166000e+00,  2.416444e+00,   7.567814e-01,
    2.087120e+00,   -2.190000e+00,  -2.166000e+00,  -2.190000e+00,  7.564814e-01,   -2.470594e+00,  2.166000e+00,   2.470594e+00,   5.697969e+00,
    9.906380e+00,   -2.190000e+00,  2.166000e+00,   -2.190000e+00,  5.643519e+00,   2.470594e+00,   2.166000e+00,   2.470594e+00,   5.697969e+00,
    2.000480e+00,   2.142000e+00,   2.166000e+00,   -2.190000e+00,  5.643519e+00,   2.470594e+00,   2.166000e+00,   -2.416444e+00,  7.567814e-01,
    -5.818780e+00,  2.142000e+00,   -2.166000e+00,  -2.190000e+00,  7.564814e-01,   -2.470594e+00,  2.166000e+00,   -2.416444e+00,  7.567814e-01,
    -5.818780e+00,  -2.190000e+00,  -2.166000e+00,  2.142000e+00,   7.564814e-01,   2.416444e+00,   2.166000e+00,   2.470594e+00,   7.567814e-01,
    2.000480e+00,   -2.190000e+00,  2.166000e+00,   2.142000e+00,   5.643519e+00,   -2.416444e+00,  2.166000e+00,   2.470594e+00,   7.567814e-01,
    9.733100e+00,   2.142000e+00,   2.166000e+00,   2.142000e+00,   5.643519e+00,   -2.416444e+00,  2.166000e+00,   -2.416444e+00,  5.589669e+00,
    1.913840e+00,   2.142000e+00,   -2.166000e+00,  2.142000e+00,   7.564814e-01,   2.416444e+00,   2.166000e+00,   -2.416444e+00,  5.589669e+00,
    2.000480e+00,   -2.190000e+00,  -2.166000e+00,  2.142000e+00,   5.643519e+00,   2.416444e+00,   -2.166000e+00,  -2.470594e+00,  7.567814e-01,
    -5.818780e+00,  -2.190000e+00,  2.166000e+00,   2.142000e+00,   7.564814e-01,   -2.416444e+00,  -2.166000e+00,  -2.470594e+00,  7.567814e-01,
    1.913840e+00,   2.142000e+00,   2.166000e+00,   2.142000e+00,   7.564814e-01,   -2.416444e+00,  -2.166000e+00,  2.416444e+00,   5.589669e+00,
    9.733100e+00,   2.142000e+00,   -2.166000e+00,  2.142000e+00,   5.643519e+00,   2.416444e+00,   -2.166000e+00,  2.416444e+00,   5.589669e+00
  };

  double *reaction = (double*)malloc(3 * nc * sizeof(double));
  double *velocity = (double*)malloc(3 * nc * sizeof(double));
  for (int i = 0; i < 3 * nc; i++)
  {
    reaction[i] = 0.0;
    velocity[i] = 0.0;
  }



  int solver_id = SICONOS_FRICTION_3D_NSN_AC; // 500
  double tolerance = 1e-10;
  int itermax = 500;

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

  printf("info: %d\n", info);

  return info;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

}
