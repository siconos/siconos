
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "LinearComplementarityProblem_as_ConvexQP.h"
#include <math.h>                          // for fmax
#include "ConvexQP.h"                      // for ConvexQP
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"                         // for DEBUG_PRINT

void Projection_ConvexQP_LCP(void *cqpIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_ConvexQP_LCP(void *cqpIn, double *x, double *PX)\n")

  ConvexQP * cqp = (ConvexQP *) cqpIn;
  LinearComplementarityProblem_as_ConvexQP* pb = (LinearComplementarityProblem_as_ConvexQP*)cqp->env;
  LinearComplementarityProblem * lcp = pb->lcp;

  int n = lcp->size;
  for(int i = 0 ; i <n ; ++i)
  {
    PX[i] = fmax(0, x[i]);
  }
}


