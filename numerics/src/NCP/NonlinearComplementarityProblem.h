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
#ifndef NCP_PROBLEM_H
#define NCP_PROBLEM_H

#include "NumericsFwd.h"
#include "SiconosConfig.h"


/*!\file NonlinearComplementarityProblem.h
 * \brief data structure to formalize a Nonlinear Complementarity Problem (NCP)
 *
 * \author Olivier Huber
*/

/** type for user defined function used to compute F and its jacobian.
 */
typedef void (*ptrFunctionNCP)(void* env, int n, double* z, double* F);
typedef void (*ptrFunctionJacNCP)(void* env, int n, double* z, NumericsMatrix* jacF);

/** \struct  NonlinearComplementarityProblem NonlinearComplementarityProblem.h
 * The structure that defines a Nonlinear Complementarity Problem (NCP) : Find two vectors \f$(z,w \in {{\mathrm{I\!R}}}^{n})\f$ such that:\n
  \f{align*}{
  w &= F(z) \\
  0 &\le w \perp z \ge 0
  \f}
 */
struct NonlinearComplementarityProblem
{
  unsigned n; /**< size of the problem */
  ptrFunctionNCP compute_F; /**< pointer to the function used to compute \f$F(z)\f$ */
  ptrFunctionJacNCP compute_nabla_F; /**< pointer to the function used to compute \f$\nabla_z F(z)\f$ */
  NumericsMatrix* nabla_F; /**< storage for \f$\nabla_z F\f$*/
  void* env; /**< environment for the compute_Fmcp and compute_nabla_F function.
               When called from Python, it contains an object with compute_F and compute_nabla_F as methods.
               When called from C, it can reference a data struct containing variables needed for the computations.*/
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** free an NCP problem 
   * \param ncp structure to free
   */
  void freeNCP(NonlinearComplementarityProblem* ncp);

  /** create an empty NCP problem
   * \return an MixedComplementarityProblem instance
   */
  NonlinearComplementarityProblem* newNCP(void);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
