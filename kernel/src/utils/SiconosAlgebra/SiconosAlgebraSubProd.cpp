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

#include <algorithm>

#include "BlockMatrix.hpp"
#include "BlockVector.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for subprod
#include "SimpleMatrix.hpp"

void siconos::algebra::subprod(const SiconosMatrix &A, const SiconosVector &x,
                               SiconosVector &y,
                               const std::vector<std::size_t> &coord,
                               bool init) {
  // To compute subY = subA * subX in an "optimized" way (in comparison with y =
  // prod(A,x) ) or subY += subA*subX if init = false.

  // coord is [r0A r1A c0A c1A r0x r1x r0y r1y]
  //
  // subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and
  // columns between c0A and c1A-1. The same for x and y with rix and riy.

  // Check dims
  auto rowA = coord[1] - coord[0];
  auto colA = coord[3] - coord[2];
  auto dimX = coord[5] - coord[4];
  auto dimY = coord[7] - coord[6];
  if (colA != dimX) THROW_EXCEPTION("inconsistent sizes between A and x.");

  if (rowA != dimY) THROW_EXCEPTION("inconsistent sizes between A and y.");

  if (dimX > x.size() || dimY > y.size() || rowA > A.size(0) ||
      colA > A.size(1))
    THROW_EXCEPTION("input index too large.");

  if (init) {
    y.segment(coord[6], dimY) =
        A.block(coord[0], coord[2], rowA, colA) * x.segment(coord[4], dimX);
  } else {
    y.segment(coord[6], dimY) +=
        A.block(coord[0], coord[2], rowA, colA) * x.segment(coord[4], dimX);
  }
}

void siconos::algebra::subprod(const SiconosMatrix &A, const BlockVector &x,
                               SiconosVector &y,
                               const std::vector<std::size_t> &coord,
                               bool init) {
  // Number of the subvector of x that handles element at position coord[4]
  auto firstBlockNum = x.getNumVectorAtPos(coord[4]);
  // Number of the subvector of x that handles element at position coord[5]
  auto lastBlockNum = x.getNumVectorAtPos(coord[5]);
  auto subCoord = coord;
  auto tmp = x.vector(firstBlockNum);
  auto subSize = static_cast<size_t>(tmp->size());  // Size of the sub-vector
  auto xTab = x.tabIndex();
  if (firstBlockNum != 0) {
    subCoord[4] -= (*xTab)[firstBlockNum - 1];
    subCoord[5] = std::min(coord[5] - (*xTab)[firstBlockNum - 1], subSize);
  } else
    subCoord[5] = std::min(coord[5], subSize);

  if (firstBlockNum == lastBlockNum) {
    subprod(A, *tmp, y, subCoord, init);
  } else {
    decltype(firstBlockNum) xPos =
        0;  // Position in x of the current sub-vector of x
    bool firstLoop = true;
    subCoord[3] = coord[2] + subCoord[5] - subCoord[4];
    for (const auto &vec : x) {
      if (xPos >= firstBlockNum && xPos <= lastBlockNum) {
        tmp = x.vector(xPos);
        if (firstLoop) {
          subprod(A, *tmp, y, subCoord, init);
          firstLoop = false;
        } else {
          subCoord[2] +=
              subCoord[5] - subCoord[4];  // !! old values for 4 and 5
          subSize = tmp->size();
          subCoord[4] = 0;
          subCoord[5] = std::min(coord[5] - (*xTab)[xPos - 1], subSize);
          subCoord[3] = subCoord[2] + subCoord[5] - subCoord[4];
          subprod(A, *tmp, y, subCoord, false);
        }
      }
      xPos++;
    }
  }
}
