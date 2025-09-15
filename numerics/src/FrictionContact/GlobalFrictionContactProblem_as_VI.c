
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
#include "GlobalFrictionContactProblem_as_VI.h"

#include <math.h>  // for sqrt

#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NM_gemv, NM_tgemv, Numeric...
#include "SiconosBlas.h"                   // for cblas_dcopy
#include "VariationalInequality.h"         // for VariationalInequality
#include "projectionOnCone.h"              // for projectionOnCone
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"  // for DEBUG_EXPR, DEBUG_BEGIN

void Function_VI_GFC3D(void *self, int n_notused, double *x, double *F) {
  DEBUG_BEGIN("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality *vi = (VariationalInequality *)self;
  GlobalFrictionContactProblem_as_VI *pb = (GlobalFrictionContactProblem_as_VI *)vi->env;
  GlobalFrictionContactProblem *gfc3d = pb->gfc3d;

  // DEBUG_EXPR(globalFrictionContact_display(gfc3d););
  int nLocal = gfc3d->dimension;

  int m = gfc3d->numberOfContacts * gfc3d->dimension;
  int n = gfc3d->M->size0;
  DEBUG_EXPR(NM_vector_display(x, n + m));

  double *globalVelocity = &x[0];
  double *reaction = &x[n];

  cblas_dcopy(n, gfc3d->q, 1, F, 1);
  for (int i = 0; i < n; i++) F[i] *= -1.0; /* F= -q*/

  NM_gemv(1.0, gfc3d->M, globalVelocity, 1.0, F); /* F= M v -q */
  NM_gemv(-1.0, gfc3d->H, reaction, 1.0, F);      /* F= M v -q - Hr  */

  /* cblas_dcopy(n , gfc3d->q , 1 , F, 1); */

  /* NM_gemv(-1.0, gfc3d->M, globalVelocity, 1.0, F); /\* F= M v -q *\/ */
  /* NM_gemv(1.0, gfc3d->H, reaction, 1.0, F); /\* F= M v -q - Hr  *\/ */

  double *localvelocity = &F[n];
  cblas_dcopy(m, gfc3d->b, 1, localvelocity, 1);             /* localvelocity = b */
  NM_tgemv(1., gfc3d->H, globalVelocity, 1., localvelocity); /* localvelocity = b + H^T V*/

  for (int contact = 0; contact < gfc3d->numberOfContacts; ++contact) {
    double normUT =
        sqrt(localvelocity[contact * nLocal + 1] * localvelocity[contact * nLocal + 1] +
             localvelocity[contact * nLocal + 2] * localvelocity[contact * nLocal + 2]);
    localvelocity[contact * nLocal] += (gfc3d->mu[contact] * normUT);
  }
  // frictionContact_display(gfc3d);

  DEBUG_EXPR(NM_vector_display(F, n + m));
  DEBUG_END("Function_VI_FC3D(void * self, double *x, double *F)\n")
}

void Projection_VI_GFC3D(void *viIn, double *x, double *PX) {
  DEBUG_BEGIN("Projection_VI_FC3D(void *viIn, double *x, double *PX)\n");

  VariationalInequality *vi = (VariationalInequality *)viIn;
  GlobalFrictionContactProblem_as_VI *pb = (GlobalFrictionContactProblem_as_VI *)vi->env;
  GlobalFrictionContactProblem *gfc3d = pb->gfc3d;
  // frictionContact_display(fc3d);

  int nLocal = gfc3d->dimension;
  int m = gfc3d->numberOfContacts * nLocal;
  int n = gfc3d->M->size0;
  DEBUG_EXPR(NM_vector_display(x, n + m));
  cblas_dcopy(n + m, x, 1, PX, 1);

  double *reaction = &PX[n];

  for (int contact = 0; contact < gfc3d->numberOfContacts; ++contact) {
    projectionOnCone(&reaction[contact * nLocal], gfc3d->mu[contact]);
  }
  DEBUG_EXPR(NM_vector_display(PX, n + m));
  DEBUG_END("Projection_VI_FC3D(void *viIn, double *x, double *PX)\n");
}
