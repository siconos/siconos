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
#ifndef FRICTIONCONTACT3D_nonsmooth_Newton_NaturalMap_H
#define FRICTIONCONTACT3D_nonsmooth_Newton_NaturalMap_H

/*!\file fc3d_nonsmooth_Newton_natural_map.h

  \brief Typedef and functions declarations related to natural map
  formulation for 3 dimension frictional contact problems in Local Coordinates.

  Subroutines used when the friction-contact 3D problem is written
  using natural map

  [...]

  where M is an n by n matrix, q an n-dimensional vector, reaction an
  n-dimensional vector and velocity an n-dimensional vector.

  We consider a "global" (ie for several contacts) problem, used to
  initialize the static global variables.


 */
#include "NumericsFwd.h"  // for SolverOptions, FrictionContactProblem
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** The natural map function signature for a 3x3 block.
   */
  typedef void (*NaturalMapFun3x3Ptr)(double* reaction,
                                             double* velocity,
                                             double mu,
                                             double* rho,
                                             double* F,
                                             double* A,
                                             double* B);

  /** Nonsmooth Newton solver based on the Natural--Map function for the
   * local (reduced) frictional contact problem in the dense form
   * \param problem the problem to solve in dense or sparse block form
   * \param reaction solution and initial guess for reaction
   * \param velocity solution and initial guess for velocity
   * \param info returned info
   * \param options  the solver options
   */
  void fc3d_nonsmooth_Newton_NaturalMap(
    FrictionContactProblem* problem,
    double *reaction,
    double *velocity,
    int *info,
    SolverOptions *options);


  /** The natural map function for several contacts.
      On each contact, the specified natural map function in iparam[9] is called.
      \param problemSize the number of contacts.
      \param computeACFun3x3 the block 3x3 natural map function.
      \param reaction3D the reactions at each contact (size: 3 x problemSize)
      \param velocity3D the velocities at each contact (size: 3 x problemSize)
      \param mu the mu parameter (size : problemSize)
      \param rho3D the rho parameters (size : 3 x problemSize)
      \param output_blocklist3 the computed natural map function (size : 3 x problemSize)
      \param output_blocklist3x3_1 the computed A part of gradient (size : 9 x problemSize)
      \param output_blocklist3x3_2 the computed B param of gradient (size : 9 x problemSize)
  */
  void fc3d_NaturalMapFunction(
    unsigned int problemSize,
    NaturalMapFun3x3Ptr computeACFun3x3,
    double *reaction3D,
    double *velocity3D,
    double *mu,
    double *rho3D,
    double *output_blocklist3,
    double *output_blocklist3x3_1,
    double *output_blocklist3x3_2);

  int fc3d_nonsmooth_Newton_NaturalMap_compute_error(
    FrictionContactProblem* problem,
    double *z , double *w, double tolerance,
    SolverOptions * options, double * error);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
