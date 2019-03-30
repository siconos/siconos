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

#include "rolling_fc3d_Solvers.h"
#include "Friction_cst.h"
#include "numerics_verbose.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>



void rolling_fc3d_set_internalsolver_tolerance(RollingFrictionContactProblem* problem,
                                               SolverOptions* options,
                                               SolverOptions* internalsolver_options,
                                               double error)
{
  int* iparam = options->iparam;
  if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] == SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE )
  {
    internalsolver_options->dparam[0] = fmax(error/options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO], options->dparam[0]/problem->numberOfContacts);
    numerics_printf_verbose(2,"fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver tolerance is set to %e\n", internalsolver_options->dparam[0] );
  }
  else if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] == SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT )
  {
    internalsolver_options->dparam[0] = error/(options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO]*problem->numberOfContacts);
    numerics_printf_verbose(2,"fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver tolerance is set to %e", internalsolver_options->dparam[0] );
  }
  else if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] == SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE)
  {
    // We use the user value for the error of the local solver
    numerics_printf_verbose(2,"fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver tolerance is set to %e", internalsolver_options->dparam[0] );
  }
  else
  {
    numerics_error("rolling_fc3d__set_internalsolver_tolerance","Unknown strategy for driving the tolerance");
  }


}
