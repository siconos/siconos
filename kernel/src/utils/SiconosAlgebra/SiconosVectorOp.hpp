/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file SiconosVectorOp.hpp
  \brief Toolbox for operators acting on vectors
  For matrix-vector op, see SiconosMatrixVectorOp.hpp
*/
#ifndef SICOVECOP_H
#define SICOVECOP_H

#include <memory>
#include <vector>

#include "BlockVector.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

/** test if two BlockVectors have the same number of blocks with
    blocks of the same size when at the same position
    \param v1 first vector to compare with
    \param v2 second vecstor to compare with
*/
bool isComparableTo(const BlockVector& v1, const BlockVector& v2);

}  // namespace siconos::algebra

#endif
