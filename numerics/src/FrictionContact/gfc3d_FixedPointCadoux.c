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

#include <assert.h>

#include "gfc3d_Gams.h"
#include "GlobalFrictionContactProblem.h"

#include "gfc3d_Solvers.h"

#include "SolverOptions.h"


/* Alart & Curnier solver for sparse global problem */
void gfc3d_FixedPointCadoux(
  GlobalFrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  double *globalVelocity,
  int *info,
  SolverOptions* options)
{

  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);
  assert(problem->H);

  assert(!problem->M->matrix0);

  /* M is square */
  assert(problem->M->size0 == problem->M->size1);

  assert(problem->M->size0 == problem->H->size0);

  unsigned int erritermax = options->iparam[7];

  if (erritermax == 0)
  {
    /* output a warning here */
    erritermax = 1;
  }

  assert(options->iparam[0] > 0);
  assert(options->iparam[3] > 0);

  assert(options->dparam[0] > 0);

  /* sparse triplet storage */
//  NM_setup(problem->M);
//  NM_setup(problem->H);

  //Mlu = (double*)
  /* save multiple data M, H */



}
