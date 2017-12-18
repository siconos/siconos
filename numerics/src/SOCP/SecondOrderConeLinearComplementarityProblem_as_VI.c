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
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "SecondOrderConeLinearComplementarityProblem_as_VI.h"

#include "projectionOnCone.h"
#include "numerics_verbose.h"
#include "SiconosBlas.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

void Function_VI_SOCLCP(void * self, int n_notused, double *x, double *F)
{
  DEBUG_PRINT("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality * vi = (VariationalInequality *) self;
  SecondOrderConeLinearComplementarityProblem_as_VI* pb = (SecondOrderConeLinearComplementarityProblem_as_VI*)vi->env;
  SecondOrderConeLinearComplementarityProblem * soclcp = pb->soclcp;
  //frictionContact_display(fc3d);

  int n =   soclcp->n;

  cblas_dcopy(n , soclcp->q , 1 , F, 1);
  NM_gemv(1.0, soclcp->M, x, 1.0, F);
}


void Projection_VI_SOCLCP(void *viIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_VI_SOCLCP(void *viIn, double *x, double *PX)\n")

  VariationalInequality * vi = (VariationalInequality *) viIn;
  SecondOrderConeLinearComplementarityProblem_as_VI* pb = (SecondOrderConeLinearComplementarityProblem_as_VI*)vi->env;
  SecondOrderConeLinearComplementarityProblem * soclcp = pb->soclcp;
  //SecondOrderConeLinearComplementarityProblem_display(soclcp);

  int cone =0;
  int n = soclcp->n;
  cblas_dcopy(n , x , 1 , PX, 1);
  int dim;
  for (cone = 0 ; cone < soclcp->nc  ; ++cone)
  {
    dim=soclcp->coneIndex[cone+1]-soclcp->coneIndex[cone];
    projectionOnSecondOrderCone(&PX[soclcp->coneIndex[cone]], soclcp->tau[cone], dim);
  }
}
