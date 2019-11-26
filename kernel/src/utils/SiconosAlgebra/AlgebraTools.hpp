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
#ifndef ALGEBRATOOLS_HPP
#define ALGEBRATOOLS_HPP

/*! \file AlgebraTools.hpp
  Utilities to perform operations on matrices (exponantial ...)
*/

class SiconosMatrix;


namespace Siconos {
  namespace algebra {
    namespace tools {

/** Compute the matrix exponential Exp = exp(A) for general matrices,
    using scaling and Padé approximation. See expm.hpp.
    \param A : input matrix
    \param Exp : result = exp(A)
    \param computeAndAdd : if true, result = result + exp(A)
**/
      void expm(SiconosMatrix& A, SiconosMatrix& Exp, bool computeAndAdd = false);
    
    } // namespace tools
  } // namespace algebra
} // namespace Siconos


#endif
