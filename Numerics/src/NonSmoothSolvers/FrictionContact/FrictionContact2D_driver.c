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
#include <string.h>
#include <time.h>
#include "LA.h"
#include "NumericsOptions.h"
#include "FrictionContact2D_Solvers.h"
#include "NonSmoothDrivers.h"
char *  SICONOS_FRICTION_2D_NSGS_STR  = "F2D_NSGS";
char *  SICONOS_FRICTION_2D_PGS_STR  = "F2D_PGS";
char *  SICONOS_FRICTION_2D_CPG_STR  = "F2D_CPG";
char *  SICONOS_FRICTION_2D_LATIN_STR  = "F2D_LATIN";
#ifdef DUMP_PROBLEM
static int fccounter = 0;
#endif

int frictionContact2D_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options, NumericsOptions* global_options)
{

#ifdef DUMP_PROBLEM
  char fname[256];
  sprintf(fname, "FrictionContactProblem%.5d.dat", fccounter++);

  FILE * foutput  =  fopen(fname, "w");
  frictionContact_printInFile(problem, foutput);
  fclose(foutput);
#endif

  if (options == NULL || global_options == NULL)
    numericsError("FrictionContact2D_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || reaction == NULL || velocity == NULL)
    numericsError("FrictionContact2D_driver", "null input for FrictionContactProblem and/or unknowns (reaction,velocity)");

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the SolverOptions structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(3, options);

  if (verbose > 0)
    printSolverOptions(options);





  /* Solver name */
  /*char * name = options->solverName;*/


  int info = -1 ;

  if (problem->dimension != 2)
    numericsError("FrictionContact2D_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");


  /* Non Smooth Gauss Seidel (NSGS) */

  if (problem->M->storageType == 1)
  {

    if (options->solverId == SICONOS_FRICTION_2D_NSGS)
    {
      if (verbose)
        printf(" ======================= Call Sparse NSGS solver for Friction-Contact 2D problem ======================\n");
      frictionContact2D_sparse_nsgs(problem, reaction , velocity , &info , options);
    }
    else
    {
      fprintf(stderr, "FrictionContact2D_driver error: unknown solver named: %s\n", idToName(options->solverId));
      exit(EXIT_FAILURE);
    }
  }
  else if (problem->M->storageType == 0)
  {

    switch (options->solverId)
    {
      /****** NLGS algorithm ******/
    case SICONOS_FRICTION_2D_PGS:
    case SICONOS_FRICTION_2D_NSGS:
    {
      if (verbose)
        printf(" ========================== Call NLGS solver for Friction-Contact 2D problem ==========================\n");
      FrictionContact2D_nsgs(problem, reaction, velocity, &info, options);
      break;
    }
    /****** CPG algorithm ******/
    case SICONOS_FRICTION_2D_CPG:
    {
      printf(" ========================== Call CPG solver for Friction-Contact 2D problem ==========================\n");
      FrictionContact2D_cpg(problem, reaction, velocity, &info, options);
      break;
    }
    /****** Latin algorithm ******/
    case SICONOS_FRICTION_2D_LATIN:
    {
      printf(" ========================== Call Latin solver for Friction-Contact 2D problem ==========================\n");
      FrictionContact2D_latin(problem, reaction, velocity, &info, options);
      break;
    }
    /*error */
    default:
    {
      fprintf(stderr, "FrictionContact2D_driver error: unknown solver named: %s\n", idToName(options->solverId));
      exit(EXIT_FAILURE);
    }
    }
  }
  else
  {
    numericsError("FrictionContact2D_driver", " error: unknown storagetype named");
    exit(EXIT_FAILURE);
  }

  return info;

}
