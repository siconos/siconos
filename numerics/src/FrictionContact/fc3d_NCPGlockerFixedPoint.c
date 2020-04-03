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

#include <stdio.h>                      // for NULL, fprintf, stderr
#include <stdlib.h>                     // for exit, EXIT_FAILURE
#include "Friction_cst.h"               // for SICONOS_FRICTION_3D_NCPGlocke...
#include "NCP_FixedP.h"                 // for Fixe
#include "NumericsFwd.h"                // for SolverOptions, FrictionContac...
#include "SolverOptions.h"              // for SolverOptions
#include "fc3d_2NCP_Glocker.h"          // for NCPGlocker_initialize, comput...
#include "fc3d_NCPGlockerFixedPoint.h"  // for F_GlockerFixedP, fc3d_FixedP_...
#include "fc3d_Solvers.h"               // for FreeSolverPtr, PostSolverPtr
#include "SiconosBlas.h"                      // for cblas_dcopy

/* Pointer to function used to update the solver, to formalize the local problem for example. */
typedef void (*UpdateSolverPtr)(int, double*);

static UpdateSolverPtr updateSolver = NULL;
static PostSolverPtr postSolver = NULL;
static FreeSolverPtr freeSolver = NULL;

/* size of a block */
static int Fsize;

/** writes \f$ F(z) \f$ using Glocker formulation
 */
void F_GlockerFixedP(int sizeF, double* reaction, double* FVector, int up2Date)
{
  /* Glocker formulation */
  double* FGlocker = NULL;
  computeFGlocker(&FGlocker, up2Date);
  /* Note that FGlocker is a static var. in fc3d2NCP_Glocker and thus there is no memory allocation in
     the present file.
  */

  /* TMP COPY: review memory management for FGlocker ...*/
  cblas_dcopy(sizeF, FGlocker, 1, FVector, 1);
  FGlocker = NULL;
}

/** writes \f$ \nabla_z F(z) \f$  using Glocker formulation and the Fischer-Burmeister function.
 */


void fc3d_FixedP_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * localsolver_options)
{

  /*
     Initialize solver (Compute F) according to the chosen formulation.
  */

  /* Glocker formulation */
  if(localsolver_options->solverId == SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint)
  {
    Fsize = 5;
    NCPGlocker_initialize(problem, localproblem);
    /*     updateSolver = &NCPGlocker_update; */
    postSolver = &NCPGlocker_post;
    freeSolver = &NCPGlocker_free;
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

int fc3d_FixedP_solve(FrictionContactProblem * localproblem, double* reaction, SolverOptions * options)
{
  int * iparam = options->iparam;
  double * dparam = options->dparam;

  double * reactionBlock = reaction;

  int info = Fixe(Fsize, reactionBlock, iparam, dparam);

  if(info > 0)
  {
    fprintf(stderr, "Numerics, fc3d_FixedP failed, reached max. number of iterations without convergence. Residual = %f\n", dparam[SICONOS_DPARAM_RESIDU]);
    exit(EXIT_FAILURE);
  }
  return info;

  /*   (*postSolver)(contact,reaction); */
}

void fc3d_FixedP_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions * localsolver_option)
{
  updateSolver = NULL;
  postSolver = NULL;
  (*freeSolver)();
}

/*
double fc3d_FixedP_computeError(int contact, int dimReaction, double* reaction, double * error)
{
  return 0.0;
}
*/
