
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "FrictionContactProblem_as_ConvexQP.h"

#include "ConvexQP.h"                // for ConvexQP
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "SolverOptions.h"           // for SolverOptions
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "SiconosBlas.h"           // for cblas_dcopy
#include "projectionOnCylinder.h"  // for projectionOnCylinder
#include "projectionOnDisk.h"      // for projectionOnDisk
#include "siconos_debug.h"         // for DEBUG_PRINT

void Projection_ConvexQP_FC3D_Cylinder(void *cqpIn, double *x, double *PX) {
  DEBUG_PRINT("Projection_ConvexQP_FC3D_Cylinder(void *cqpIn, double *x, double *PX)\n")

  ConvexQP *cqp = (ConvexQP *)cqpIn;
  FrictionContactProblem_as_ConvexQP *pb = (FrictionContactProblem_as_ConvexQP *)cqp->env;
  FrictionContactProblem *fc3d = pb->fc3d;
  SolverOptions *options = pb->options;
  // frictionContact_display(fc3d);

  int contact = 0;
  int nLocal = fc3d->dimension;
  int n = fc3d->numberOfContacts * nLocal;
  cblas_dcopy(n, x, 1, PX, 1);
  for (contact = 0; contact < fc3d->numberOfContacts; ++contact) {
    projectionOnCylinder(&PX[contact * nLocal], options->dWork[contact]);
  }
}

void Projection_ConvexQP_FC3D_Disk(void *cqpIn, double *x, double *PX) {
  DEBUG_PRINT("Projection_ConvexQP_FC3D_Cylinder(void *cqpIn, double *x, double *PX)\n")

  ConvexQP *cqp = (ConvexQP *)cqpIn;
  FrictionContactProblem_as_ConvexQP *pb = (FrictionContactProblem_as_ConvexQP *)cqp->env;
  FrictionContactProblem *fc3d = pb->fc3d;
  SolverOptions *options = pb->options;
  // frictionContact_display(fc3d);

  int nLocal = 2;
  int n = fc3d->numberOfContacts * nLocal;

  cblas_dcopy(n, x, 1, PX, 1);

  for (int contact = 0; contact < fc3d->numberOfContacts; ++contact) {
    projectionOnDisk(&PX[contact * nLocal], options->dWork[contact]);
  }
}
