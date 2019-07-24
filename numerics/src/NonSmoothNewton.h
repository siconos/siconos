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
#ifndef NONSMOOTHNEWTON_H
#define NONSMOOTHNEWTON_H

/*!\file NonSmoothNewton.h
  Typedef and functions declarations related to non-smooth Newton solver

  Solve \f$ \phi(z) = 0 \f$ using a Newton method.

  The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
  for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

 */

/* Pointer to function that corresponds to the function \f$ \phi \f$ */
typedef void (*NewtonFunctionPtr)(int, double*, double*, int);

#include "SiconosConfig.h"
#include "SolverOptions.h"

enum NONSMOOTH_NEWTON_SOLVER
{
  SICONOS_NONSMOOTH_NEWTON_LSA = 11000
};

extern const char* const   SICONOS_NONSMOOTH_NEWTON_LSA_STR ;


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Armijo Linesearch
   *  \param n size of the vector z
   *  \param z unknown vector
   *  \param dir search direction
   *  \param psi_k initial value of the merit function
   *  \param descentCondition descent condition
   *  \param phi pointer to function used to compute phi(z)
   */
  void linesearch_Armijo(int n, double *z, double* dir, double psi_k,
                         double descentCondition, NewtonFunctionPtr* phi);


  /** Newton solver with line Search
  \param n size of the vector z
  \param z unknown vector, in-out argument
  \param phi pointer to \f$ \phi \f$ function
  \param jacobianPhi pointer to \f$ \nabla_z \phi(z) \f$ function
  \param iparam vector of int parameters:
   - [0] : max. number of iterations
   - [1] : number of iterations processed
  \param dparam vector of double parameters:
   - [0]: tolerance
   - [1]: error
  \return int 0 if ok
  */
  int nonSmoothNewton(int n, double* z, NewtonFunctionPtr* phi,
                      NewtonFunctionPtr* jacobianPhi,
                      SolverOptions * options);

  /** Newton solver without line Search
  \param n size of the vector z
  \param z unknown vector, in-out argument
  \param phi pointer to \f$ \phi \f$ function
  \param jacobianPhi pointer to \f$ \nabla_z \phi(z) \f$ function
  \param iparam vector of int parameters:
   - [0] : max. number of iterations
   - [1] : number of iterations processed
  \param dparam vector of double parameters:
   - [0]: tolerance
   - [1]: error
  \return int 0 if ok
  */
  int nonSmoothDirectNewton(int n, double* z, NewtonFunctionPtr* phi,
                            NewtonFunctionPtr* jacobianPhi,
                            SolverOptions * options);
  
  void nonSmoothNewton_setDefaultSolverOptions(SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
