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
#ifndef FRICTIONCONTACTPROBLEM_AS_VI_H
#define FRICTIONCONTACTPROBLEM_AS_VI_H

/*! \page fcProblemVI Friction-contact problems (2D or 3D) as VI
 *
 *
 */


/*!\file FrictionContactProblem_as_VI.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
*/

#include "NumericsFwd.h"
#include "SiconosConfig.h"

/** \struct FrictionContactProblem_as_VI FrictionContactProblem_as_VI.h
 *
 */
struct FrictionContactProblem_as_VI
{
  /* the VI associated with the FC3D problem */
  VariationalInequality * vi;
  /* the FC3D associated with the VI  */
  FrictionContactProblem * fc3d;
  /* the SolverOptions that might be used to pass some numerical parameters */
  SolverOptions * options;
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void Function_VI_FC3D(void * self, int n, double *x, double *F);

  void Projection_VI_FC3D(void *viIn, double *x, double *PX);

  void Function_VI_FC3D_Cylinder(void * self, int n, double *x, double *F);

  void Projection_VI_FC3D_Cylinder(void *viIn, double *x, double *PX);

  void Function_VI_FC3D_Disk(void * self, int n, double *x, double *F);

  void Projection_VI_FC3D_Disk(void *viIn, double *x, double *PX);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
