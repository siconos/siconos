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
#include <string.h>
#include <time.h>
#include <float.h>
#include "LA.h"
#include "NumericsOptions.h"
#include "FrictionContact2D_Solvers.h"
#include "NonSmoothDrivers.h"

int frictionContact2D_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  int info = -1;
  switch (solverId)
  {
  case SICONOS_FRICTION_2D_NSGS:
  {
    info =    frictionContact2D_sparse_nsgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_2D_PGS:
  case SICONOS_FRICTION_2D_NLGS:
  {
    info =    frictionContact2D_nlgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_2D_CPG:
  {
    info =    frictionContact2D_cpg_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_FRICTION_2D_LATIN:
  {
    info =    frictionContact2D_latin_setDefaultSolverOptions(options);
    break;
  }

  default:
  {
    numericsError("frictionContact2D_setDefaultSolverOptions", "Unknown Solver");

  }
  }

  return info;
}
