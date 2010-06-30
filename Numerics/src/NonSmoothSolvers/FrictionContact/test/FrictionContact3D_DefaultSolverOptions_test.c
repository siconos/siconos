/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
  printf("\n Start of test on Default SolverOptions\n");
  int info = 0 ;
  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));



  info = frictionContact3D_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_NSGS);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = frictionContact3D_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_NSGSV);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = frictionContact3D_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_PROX);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = frictionContact3D_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_TFP);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = frictionContact3D_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_DSFP);
  printSolverOptions(options);
  deleteSolverOptions(options);

  free(options);

  printf("\n End of test on Default SolverOptions\n");
  return info;
}
