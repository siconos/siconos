// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "FrictionContactProblem_as_ConvexQP.h"

#include "ConvexQP.h"                // for ConvexQP
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "SolverOptions.h"           // for SolverOptions
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "SiconosBlas.h"  // for cblas_dcopy
#include "numerics_errors.h"
#include "projectionOnCylinder.h"  // for projectionOnCylinder
#include "projectionOnDisk.h"      // for projectionOnDisk
#include "siconos_debug.h"         // for DEBUG_PRINT
void Projection_ConvexQP_FC3D_Cylinder(void* cqpIn, double* x, double* PX) {
  DEBUG_PRINT("Projection_ConvexQP_FC3D_Cylinder(void *cqpIn, double *x, double *PX)\n")

  ConvexQP* cqp = (ConvexQP*)cqpIn;
  FrictionContactProblem_as_ConvexQP* pb = (FrictionContactProblem_as_ConvexQP*)cqp->env;
  FrictionContactProblem* fc3d = pb->fc3d;
  SolverOptions* options = pb->options;
  // frictionContactProblem_display(fc3d);

  int contact = 0;
  int nLocal = fc3d->dimension;
  int n = fc3d->numberOfContacts * nLocal;
  cblas_dcopy(n, x, 1, PX, 1);
  for (contact = 0; contact < fc3d->numberOfContacts; ++contact) {
    projectionOnCylinder(&PX[contact * nLocal], options->dWork[contact]);
  }
}

void Projection_ConvexQP_FC3D_Disk(void* cqpIn, double* x, double* PX) {
  DEBUG_PRINT("Projection_ConvexQP_FC3D_Cylinder(void *cqpIn, double *x, double *PX)\n")

  ConvexQP* cqp = (ConvexQP*)cqpIn;
  FrictionContactProblem_as_ConvexQP* pb = (FrictionContactProblem_as_ConvexQP*)cqp->env;
  FrictionContactProblem* fc3d = pb->fc3d;
  SolverOptions* options = pb->options;
  // frictionContactProblem_display(fc3d);

  int nLocal = 2;
  int n = fc3d->numberOfContacts * nLocal;

  cblas_dcopy(n, x, 1, PX, 1);

  for (int contact = 0; contact < fc3d->numberOfContacts; ++contact) {
    projectionOnDisk(&PX[contact * nLocal], options->dWork[contact]);
  }
}
