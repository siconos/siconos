/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#ifndef NCPPath_H
#define NCPPath_H

#include "Standalone_Path.h"
#include "PathAlgebra.h"
/*!\file NCP_Path.h

  \brief Interface to Path Solver for NCP problems

  Solves the following Nonlinear Complementarity Problem with Path Solver:
  \f{eqnarray*}
  0 \le F(z) \perp z \ge 0 \\
  \f}

  \author Franck Perignon, 21/05/2008
*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Path solver for NCP problem
      \param n size of the vector z
      \param z vector
      \param F pointer to function used to compute \f$ F(z) \f$
      \param jacobianF pointer to function used to compute \f$ \nabla_zF(z) \f$
      \param iparam vector of int parameters (useless at the time)
      \param dparam vector of double parameters (useless at the time)
      \return 0 if successfull
  */
  int NCP_Path(int n, double* z, FuncEvalPtr F, JacEvalPtr jacobianF, int* iparam, double* dparam);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
