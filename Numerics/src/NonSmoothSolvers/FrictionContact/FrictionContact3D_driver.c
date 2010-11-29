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
#include "FrictionContact3D_Solvers.h"
#include "NonSmoothDrivers.h"

char  SICONOS_FRICTION_3D_NSGS_STR [] = "F3D_NSGS";
char  SICONOS_FRICTION_3D_NSGSV_STR [] = "F3D_NSGSV";
char  SICONOS_FRICTION_3D_TFP_STR [] = "F3D_TFP";
char  SICONOS_FRICTION_3D_GLOBALAC_STR [] = "F3D_GLOBALAC";
char  SICONOS_FRICTION_3D_DSFP_STR[] = "F3D_DSFP";
char  SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint_STR[] = "F3D_NCPGlockerFBFixedPoint";
char  SICONOS_FRICTION_3D_AlartCurnierNewton_STR[] = "F3D_AlartCurnierNewton";
char  SICONOS_FRICTION_3D_DampedAlartCurnierNewton_STR[] = "F3D_DampedAlartCurnierNewton";
char  SICONOS_FRICTION_3D_NCPGlockerFBNewton_STR[] = "F3D_NCPGlockerFBNewton";
char SICONOS_FRICTION_3D_ProjectionOnConeWithDiagonalization_STR[] = "F3D_ProjectionOnConeWithDiagonalization";
char SICONOS_FRICTION_3D_ProjectionOnCone_STR[] = "F3D_ProjectionOnCone";
char SICONOS_FRICTION_3D_ProjectionOnConeWithLocalIteration_STR[] = "F3D_ProjectionOnConeWithLocalIteration";
char SICONOS_FRICTION_3D_projectionOnConeWithRegularization_STR[] = "F3D_projectionOnConeWithRegularization";
char SICONOS_FRICTION_3D_NCPGlockerFBPATH_STR[] = "F3D_NCPGlockerFBPATH";
char SICONOS_FRICTION_3D_projectionOnCylinder_STR[] = "F3D_projectionOnCylinder";
char SICONOS_FRICTION_3D_ProjectionOnCone_velocity_STR[] = "F3D_ProjectionOnCone_velocity";
char SICONOS_FRICTION_3D_PGoC_STR[] = "F3D_PGoC";
char SICONOS_FRICTION_3D_DeSaxceFixedPoint_STR[] = "F3D_DeSaxceFixedPoint";
char SICONOS_FRICTION_3D_PROX_STR[] = "F3D_PROX";
char SICONOS_FRICTION_3D_QUARTIC_STR[] = "F3D_QUARTIC";
char SICONOS_FRICTION_3D_QUARTIC_NU_STR[] = "F3D_QUARTIC_NU";

int frictionContact3D_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options, NumericsOptions* global_options)
{
  if (options == NULL)
    numericsError("FrictionContact3D_driver", "null input for solver and/or global options");
  /* Set global options */
  if (global_options)
    setNumericsOptions(global_options);

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the SolverOptions structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(3, options);

  if (verbose > 0)
    printSolverOptions(options);

  /* Solver name */
  /*char * name = options->solverName;*/

  int info = -1 ;

  if (problem->dimension != 3)
    numericsError("FrictionContact3D_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Check for trivial case */
  info = checkTrivialCase(problem, velocity, reaction, iparam, dparam);


  if (info == 0)
    return info;


  switch (options->solverId)
  {
    /* Non Smooth Gauss Seidel (NSGS) */
  case SICONOS_FRICTION_3D_NSGS:
  {
    if (verbose == 1)
      printf(" ========================== Call NSGS solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_nsgs(problem, reaction , velocity , &info , options);
    break;
  }
  case SICONOS_FRICTION_3D_NSGSV:
  {
    if (verbose == 1)
      printf(" ========================== Call NSGSV solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_nsgs_velocity(problem, reaction , velocity , &info , options);
    break;
  }
  /* Proximal point algorithm */
  case SICONOS_FRICTION_3D_PROX:
  {
    if (verbose == 1)
      printf(" ========================== Call PROX (Proximal Point) solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_proximal(problem, reaction , velocity , &info , options);
    break;
  }
  /* Tresca Fixed point algorithm */
  case SICONOS_FRICTION_3D_TFP:
  {
    if (verbose == 1)
      printf(" ========================== Call TFP (Tresca Fixed Point) solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_TrescaFixedPoint(problem, reaction , velocity , &info , options);
    break;
  }
  /* Projected Gradient algorithm */
  case SICONOS_FRICTION_3D_DSFP:
  {
    if (verbose == 1)
      printf(" ========================== Call DeSaxce Fized Point (DSFP) solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_DeSaxceFixedPoint(problem, reaction , velocity , &info , options);
    break;
  }
  /* Global Alart Curnier */
  case SICONOS_FRICTION_3D_GLOBALAC:
  {
    if (verbose == 1)
      printf(" ========================== Call Global Alart Curnier solver for Friction-Contact 3D problem ==========================\n");
    if (problem->M->matrix0)
    {
      frictionContact3D_globalAlartCurnier(problem, reaction , velocity , &info , options);
    }
    else
    {
      frictionContact3D_sparseGlobalAlartCurnier(problem, reaction , velocity , &info , options);
    }
    break;
  }
  case SICONOS_FRICTION_3D_QUARTIC_NU:
  case SICONOS_FRICTION_3D_QUARTIC:
  {
    if (verbose == 1)
      printf(" ========================== Call Quartic solver for Friction-Contact 3D problem ==========================\n");
    frictionContact3D_unitary_enumerative(problem, reaction , velocity , &info , options);
    break;
  }
  case SICONOS_FRICTION_3D_AlartCurnierNewton:
  case SICONOS_FRICTION_3D_DampedAlartCurnierNewton:
  {
    if (verbose == 1)
      printf(" ========================== Call Quartic solver for Friction-Contact 3D problem ==========================\n");
    info = frictionContact3D_Newton_solve(problem, reaction , options);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, FrictionContact3D_driver failed. Unknown solver.\n");
    exit(EXIT_FAILURE);

  }
  }
  return info;

}

int checkTrivialCase(FrictionContactProblem* problem, double* velocity, double* reaction, int* iparam, double* dparam)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  /* Dimension of the problem */
  int n = 3 * nc;
  int i = 0;
  /*take off? R=0 ?*/
  for (i = 0; i < nc; i++)
  {
    if (q[3 * i] < -DBL_EPSILON)
      return -1;
  }
  for (i = 0 ; i < n ; ++i)
  {
    velocity[i] = q[i];
    reaction[i] = 0.;
  }
  iparam[2] = 0;
  iparam[4] = 0;
  dparam[1] = 0.0;
  dparam[3] = 0.0;
  if (verbose == 1)
    printf("FrictionContact3D driver, take off, trivial solution reaction = 0, velocity = q.\n");
  return 0;
}

