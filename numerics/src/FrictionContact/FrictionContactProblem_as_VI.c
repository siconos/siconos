
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
#include "FrictionContactProblem_as_VI.h"
#include "FrictionContactProblem.h"
#include "VariationalInequality.h"
#include "projectionOnCone.h"
#include "projectionOnCylinder.h"
#include "projectionOnDisk.h"
#include "numerics_verbose.h"
#include "SiconosBlas.h"
#include "NumericsMatrix.h"
#include "SolverOptions.h"
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"




void Function_VI_FC3D(void * self, int n_notused, double *x, double *F)
{
  DEBUG_PRINT("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality * vi = (VariationalInequality *) self;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int nLocal =  fc3d->dimension;
  int n = fc3d->numberOfContacts *  fc3d->dimension;

  cblas_dcopy(n , fc3d->q , 1 , F, 1);
  NM_gemv(1.0, fc3d->M, x, 1.0, F);

  int contact =0;

  for (contact = 0 ; contact <  fc3d->numberOfContacts ; ++contact)
  {
    double  normUT = sqrt(F[contact * nLocal + 1] * F[contact * nLocal + 1]
                            + F[contact * nLocal + 2] * F[contact * nLocal + 2]);
    F[contact * nLocal] +=  (fc3d->mu[contact] * normUT);
  }

}


void Projection_VI_FC3D(void *viIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_VI_FC3D(void *viIn, double *x, double *PX)\n")

  VariationalInequality * vi = (VariationalInequality *) viIn;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int contact =0;
  int nLocal =  fc3d->dimension;
  int n = fc3d->numberOfContacts* nLocal;
  cblas_dcopy(n , x , 1 , PX, 1);
  for (contact = 0 ; contact < fc3d->numberOfContacts  ; ++contact)
  {
    projectionOnCone(&PX[ contact * nLocal ], fc3d->mu[contact]);
  }
}

void Function_VI_FC3D_Cylinder(void * self, int n_notused, double *x, double *F)
{
  DEBUG_PRINT("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality * vi = (VariationalInequality *) self;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int n = fc3d->numberOfContacts *  fc3d->dimension;

  cblas_dcopy(n , fc3d->q , 1 , F, 1);
  NM_gemv(1.0, fc3d->M, x, 1.0, F);
}

void Projection_VI_FC3D_Cylinder(void *viIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_VI_FC3D_Cylinder(void *viIn, double *x, double *PX)\n")

  VariationalInequality * vi = (VariationalInequality *) viIn;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  SolverOptions * options = pb->options;
  //frictionContact_display(fc3d);

  int contact =0;
  int nLocal =  fc3d->dimension;
  int n = fc3d->numberOfContacts* nLocal;
  cblas_dcopy(n , x , 1 , PX, 1);
  for (contact = 0 ; contact < fc3d->numberOfContacts  ; ++contact)
  {
    projectionOnCylinder(&PX[ contact * nLocal ], options->dWork[contact]);
  }
}



void Function_VI_FC3D_Disk(void * self, int n_notused, double *x, double *F)
{
  DEBUG_PRINT("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality * vi = (VariationalInequality *) self;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  SolverOptions * options = pb->options;
  SplittedFrictionContactProblem * splitted_problem = (SplittedFrictionContactProblem *)options->solverData;

  int n = 2 * fc3d->numberOfContacts ;

  cblas_dcopy(n , splitted_problem->q_t , 1 , F, 1);
  NM_gemv(1.0, splitted_problem->M_tt, x, 1.0, F);
}

void Projection_VI_FC3D_Disk(void *viIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_VI_FC3D_Cylinder(void *viIn, double *x, double *PX)\n")

  VariationalInequality * vi = (VariationalInequality *) viIn;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  SolverOptions * options = pb->options;
  //frictionContact_display(fc3d);

  int nLocal  =  2;
  int n = fc3d->numberOfContacts* nLocal;

  cblas_dcopy(n , x , 1 , PX, 1);

  for (int contact = 0 ; contact < fc3d->numberOfContacts  ; ++contact)
  {
    projectionOnDisk(&PX[ contact * nLocal ], options->dWork[contact]);
  }
}
