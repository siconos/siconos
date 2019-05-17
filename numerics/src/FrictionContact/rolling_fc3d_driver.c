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
#include <string.h>
#include <time.h>
#include <float.h>
#include <assert.h>

#include "rolling_fc3d_Solvers.h"
#include "NonSmoothDrivers.h"
#include "numerics_verbose.h"

const char* const   SICONOS_ROLLING_FRICTION_3D_NSGS_STR = "RFC3D_NSGS";

const char* const  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration_STR = "RFC3D_ProjectionOnConeWithLocalIteration";

const char* const  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone_STR = "RFC3D_ProjectionOnCone";

int rolling_fc3d_driver(RollingFrictionContactProblem* problem,
                        double *reaction, double *velocity,
                        SolverOptions* options)
{
  /* verbose=3; */
  /* rollingFrictionContact_display(problem); */
  /* rollingFrictionContact_printInFilename(problem, "rfc3d_sphere_1.dat"); */


  if (options == NULL)
    numerics_error("rolling_fc3d_driver", "null input for solver options");

  assert(options->isSet); /* true(1) if the SolverOptions structure has been filled in else false(0) */

  if (verbose > 1)
    solver_options_print(options);

  int info = -1 ;

  if (problem->dimension != 5)
    numerics_error("rolling_fc3d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* Check for trivial case */
  info = rolling_fc3d_checkTrivialCase(problem, velocity, reaction, options);
  if (info == 0)
  {
    /* If a trivial solution is found, we set the number of iterations to 0
       and the reached acuracy to 0.0 .
    */
    options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
    goto exit;
  }


  switch (options->solverId)
  {
    /* Non Smooth Gauss Seidel (NSGS) */
  case SICONOS_ROLLING_FRICTION_3D_NSGS:
  {
    numerics_printf(" ========================== Call NSGS solver for Rolling Friction-Contact 3D problem ==========================\n");
    rolling_fc3d_nsgs(problem, reaction , velocity , &info , options);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, rolling_fc3d_driver failed. Unknown solver.\n");
    exit(EXIT_FAILURE);

  }
  }

exit:

  return info;

}

int rolling_fc3d_checkTrivialCase(RollingFrictionContactProblem* problem, double* velocity,
                                  double* reaction, SolverOptions* options)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  /* Dimension of the problem */
  int n = 5 * nc;
  int i = 0;
  /*take off? R=0 ?*/
  for (i = 0; i < nc; i++)
  {
    if (q[5 * i] < -DBL_EPSILON)
      return -1;
  }
  for (i = 0 ; i < n ; ++i)
  {
    velocity[i] = q[i];
    reaction[i] = 0.;
  }

  numerics_printf("fc3d fc3d_checkTrivialCase, take off, trivial solution reaction = 0, velocity = q.\n");
  return 0;
}
