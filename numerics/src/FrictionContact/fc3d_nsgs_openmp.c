/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_Path.h"
#include "fc3d_NCPGlockerFixedPoint.h"
#include "fc3d_projection.h"
#include "fc3d_unitary_enumerative.h"
#include "fc3d_compute_error.h"
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <alloca.h>

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

#include <SiconosConfig.h>
#if defined(WITH_OPENMP) && defined(_OPENMP)
#define USE_OPENMP 1
#include <omp.h>
#endif

void fc3d_nsgs_openmp_redblack(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options);

void fc3d_nsgs_openmp_for(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options);


void fc3d_nsgs_openmp(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;


  if (iparam[11] == 0)
  {
    fc3d_nsgs_openmp_for(problem, reaction, velocity, info, options) ;
  }
  else if (iparam[11] == 1)
  {
    fc3d_nsgs_openmp_redblack(problem, reaction, velocity, info, options) ;
  }
  else
  {
    numericsError("fc3d_nsgs_openmp", "The method defined by iparam[11] is not recognized");
  }

}



int fc3d_nsgs_openmp_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS OPENMP Solver\n");
  }

  /*  strcpy(options->solverName,"NSGS");*/
  options->solverId = SICONOS_FRICTION_3D_NSGS_OPENMP;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 15;
  options->dSize = 15;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  null_SolverOptions(options);
  for (i = 0; i < 15; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_onecontact_nonsmooth_Newtow_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
