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

#ifndef SM_SETBLOCK_HPP
#define SM_SETBLOCK_HPP

#include "SiconosAlgebraTypeDef.hpp"

// NOT A FRIEND
/** Copy a subBlock of MIn into a sub-block of MOut - Dim and positions of the sub-block are given in dim and start.
 * \param MIn a SPC::SiconosMatrix
 * \param[in,out] MOut a SP::SiconosMatrix
 * \param dim an Index, dim[0], dim[1]: number of rows and columns of the sub-block
 * \param start an Index, start[0], start[1]: position (row, column) of the first element of the sub-block in MIn
 *  start[2], start[3]: position (row, column) of the first element of the sub-block in MOut.
 */
void setBlock(SPC::SiconosMatrix MIn, SP::SiconosMatrix MOut, const Index& dim, const Index& start);

#endif
