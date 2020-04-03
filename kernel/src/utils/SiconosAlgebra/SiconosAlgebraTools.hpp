/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*! \file SiconosAlgebraTools.hpp
  \brief Standalone functions used with matrices and vectors.
*/
#ifndef SICONOSALGEBRATOOLS_H
#define SICONOSALGEBRATOOLS_H

class SiconosMatrix;
class BlockVector;

namespace Siconos {
  namespace Algebra {


    /** test if two BlockVectors have the same number of blocks with
        blocks of the same size when at the same position
        \param v1 first vector to compare with
        \param v2 second vector to compare with
    */
    bool isComparableTo(const BlockVector& v1, const BlockVector& v2);

    /** test if two matrices have the same number of blocks with
        blocks of the same dimension when at the same position
        \param v1 first matrix to compare with
        \param v2 second matrix to compare with
    */
    bool isComparableTo(const SiconosMatrix& m1, const SiconosMatrix& m2);
  }

}

#endif
