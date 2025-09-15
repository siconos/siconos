
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
#include "ConvexQP_as_VI.h"

#include "ConvexQP.h"               // for ConvexQP
#include "NumericsMatrix.h"         // for NM_gemv
#include "VariationalInequality.h"  // for VariationalInequality
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "SiconosBlas.h"    // for cblas_dcopy
#include "siconos_debug.h"  // for DEBUG_PRINT

void Function_VI_CQP(void *self, int n_notused, double *x, double *F) {
  DEBUG_PRINT("Function_VI_CQP(void * self, double *x, double *F)\n")
  VariationalInequality *vi = (VariationalInequality *)self;
  ConvexQP_as_VI *pb = (ConvexQP_as_VI *)vi->env;
  ConvexQP *cqp = pb->cqp;

  int n = cqp->size;

  cblas_dcopy(n, cqp->q, 1, F, 1);
  NM_gemv(1.0, cqp->M, x, 1.0, F);
}

void Projection_VI_CQP(void *viIn, double *x, double *PX) {
  DEBUG_PRINT("Projection_VI_CQP(void *viIn, double *x, double *PX)\n")

  VariationalInequality *vi = (VariationalInequality *)viIn;
  ConvexQP_as_VI *pb = (ConvexQP_as_VI *)vi->env;
  ConvexQP *cqp = pb->cqp;

  cqp->ProjectionOnC(cqp, x, PX);
}
