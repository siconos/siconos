
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

#include "GlobalFrictionContactProblem_as_ConvexQP.h"
#include "ConvexQP.h"                      // for ConvexQP
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "projectionOnCone.h"              // for projectionOnCone
#include "SiconosBlas.h"                         // for cblas_dcopy
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"


void Projection_ConvexQP_GFC3D_DualCone(void *cqpIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_ConvexQP_FC3D_Cylinder(void *cqpIn, double *x, double *PX)\n")

  ConvexQP * cqp = (ConvexQP *) cqpIn;
  GlobalFrictionContactProblem_as_ConvexQP* pb = (GlobalFrictionContactProblem_as_ConvexQP*)cqp->env;
  GlobalFrictionContactProblem * gfc3d = pb->gfc3d;

  //globalFrictionContact_display(gfc3d);

  int contact =0;
  int nLocal =  gfc3d->dimension;
  int n = gfc3d->numberOfContacts* nLocal;

  cblas_dcopy(n, x, 1, PX, 1);
  for(contact = 0 ; contact < gfc3d->numberOfContacts  ; ++contact)
  {
    projectionOnCone(&PX[ contact * nLocal ], 1.0/gfc3d->mu[contact]);
  }
}

