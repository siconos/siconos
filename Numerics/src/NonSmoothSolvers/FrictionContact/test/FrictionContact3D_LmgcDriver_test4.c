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
#include "frictionContact_test_function.h"


int main(void)
{
  int info = 0 ;


  int nc = 4 ;
  int nb = 16 ;
  double mu[4] =
  {
    7.500000e-01,   7.500000e-01,   7.500000e-01,   7.500000e-01
  };
  int row[16] =
  {
    1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,   4
  };
  int column[16] =
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



  int solver_id = SICONOS_FRICTION_3D_GLOBALAC; // 500
  double tolerance = 1e-10;
  int itermax = 500;

  info = frictionContact3D_LmgcDriver(reaction,
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
                                      0, 0);
  printf("reaction:");
  printm(1, 3 * nc, reaction);

  printf("velocity:");
  printm(1, 3 * nc, velocity);

  free(reaction);
  free(velocity);

  printf("info: %d\n", info);

  return info;
}
