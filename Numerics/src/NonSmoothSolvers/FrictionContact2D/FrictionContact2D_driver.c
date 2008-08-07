/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "LA.h"
#include "Numerics_Options.h"
#include "FrictionContact2D_Solvers.h"
#include "NonSmoothDrivers.h"

int frictionContact2D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options)
{

  /* Set global options */
  setNumericsOptions(global_options);

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the Solver_Options structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(3, options);

  if (verbose > 0)
    printSolverOptions(options);

  /* Solver name */
  char * name = options->solverName;


  int info = -1 ;

  /* Non Smooth Gauss Seidel (NSGS) */
  if (strcmp(name, "NSGS") == 0)
  {
    if (problem->M->storageType == 1)
    {
      if (verbose)
        printf(" ======================= Call Sparse NSGS solver for Friction-Contact 2D problem ======================\n");
      frictionContact2D_sparse_nsgs(problem, reaction , velocity , &info , options);
    }
    else
    {
      if (verbose)
        printf(" ========================== Call NSGS solver for Friction-Contact 2D problem ==========================\n");

    };

    if (strcmp(name, "NSGS") == 0)
    {
      if (verbose)
        printf(" ========================== Call NSGS solver for Friction-Contact 2D problem ==========================\n");
      /*frictionContact2D_nsgs(problem, reaction , velocity , &info , options);*/
    }

    return info;


  }
}
