/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include "NumericsOptions.h"
#include "VariationalInequality_Solvers.h"
#include "NonSmoothDrivers.h"

char *  SICONOS_VI_EG_STR = "VI_EG";

void snPrintf(int level, SolverOptions* opts, const char *fmt, ...);

int variationalInequality_driver(VariationalInequality* problem, 
                                 double *x, double *fx, 
                                 SolverOptions* options, 
                             NumericsOptions* global_options)
{
  if (options == NULL)
    numericsError("variationalInequality_driver", "null input for solver and/or global options");
  /* Set global options */
  if (global_options)
  {
    setNumericsOptions(global_options);
    options->numericsOptions = global_options;
  }

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the SolverOptions structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(3, options);

  if (verbose > 0)
    printSolverOptions(options);

  /* Solver name */
  /*char * name = options->solverName;*/

  int info = -1 ;

  /* Check for trivial case */
  info = checkTrivialCase_vi(problem, x, fx, options);

  if (info == 0)
    return info;


  switch (options->solverId)
  {
    /* Non Smooth Gauss Seidel (NSGS) */
  /* Extra Gradient algorithm */
  case SICONOS_VI_EG:
  {
    snPrintf(1, options, 
             " ========================== Call ExtraGradient (EG) solver for VI problem ==========================\n");
    variationalInequality_ExtraGradient(problem, x , fx , &info , options);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, variationalInequality_driver failed. Unknown solver.\n");
    exit(EXIT_FAILURE);

  }
  }

  return info;

}

int checkTrivialCase_vi(VariationalInequality* problem, double* x, 
                     double* fx, SolverOptions* options)
{
  
  if (verbose == 1)
    printf("variationalInequality driver, take off, trivial solution reaction = 0, velocity = q.\n");
  return 0;
}
