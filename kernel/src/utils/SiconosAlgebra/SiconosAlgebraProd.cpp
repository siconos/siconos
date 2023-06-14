/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "BlockVector.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"        // For prod
#include "SiconosMatrixVectorOp.hpp"  // For prod
#include "SiconosVector.hpp"
#include "SiconosVectorIterator.hpp"
#include "SiconosVectorOp.hpp"
#include "SimpleMatrix.hpp"




void siconos::algebra::prod(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y,
                            bool init) {

  unsigned int startRow = 0;
  // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix
  // of A corresponding to y[i] position.
  //       // private_prod takes into account the fact that x and y[i] may be block vectors.
  for (auto& it : y) {
    // A.private_prod(startRow, x, *it, init);
    if(init)
      it->zero();
    unsigned int startRow = 0;
    it->noalias() += A.block(startRow, 0, it->size(), x.size()) * *it;
    startRow += it->size();
  }
}

void siconos::algebra::prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y,
                            bool init) {
  if (init) y.zero();
  unsigned int startRow = 0;
  unsigned int startCol = 0;
  // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA x[i]),
  // with subA a submatrix of A, starting from position startRow in rows and startCol in
  // columns. private_prod takes also into account the fact that each block of x can also be a
  // block.
  for (auto& it : x) {
    assert(y != *it);
    y += A.block(startRow, startCol, y.size(), it->size()) * *it;
    startCol += it->size();
  }
}


void siconos::algebra::prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y,
                            bool init) {
  // To compute y = A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += A*x if init = false.

  if (A.size(1) != x.size()) THROW_EXCEPTION("inconsistent sizes between A and x.")

  if (A.size(0) != y.size()) THROW_EXCEPTION("inconsistent sizes between A and y.");

  
  // === First case: y is not a block vector ===
  if (init) {
    if (&x != &y)  // if no common memory between x and y.
    {
      assert(y != x);
      y.noalias() = A*x;
    } else  // if x and y are the same object => alias
    {
      y = A*x;
    }
  } else  // += case
  {
    if (&x != &y)  // if no common memory between x and y.
    {
      y.noalias() += A*x;
    } else  // if x and y are the same object => alias
    {
      y += A*x;
    }
  }
}

void siconos::algebra::prod(const SiconosVector& x, const SiconosMatrix& A, SiconosVector& y,
                            bool init) {
  // To compute y = trans(A) * x in an "optimized" way, if init = true
  // (or y = trans(A) * x + y if init = false

  if (A.size(0) != x.size()) THROW_EXCEPTION("inconsistent sizes between A and x.");

  if (A.size(1) != y.size()) THROW_EXCEPTION("inconsistent sizes between A and y.");

  if (init) {
    if (&x != &y)  // if no common memory between x and y.
    {
      y.noalias() = A.transpose() * x;
    } else  // if x and y are the same object => alias
    {
      y = A.transpose() * x;
    }
  } else  // += case
  {
    if (&x != &y)  // if no common memory between x and y.
    {
      y.noalias() += A.transpose() * x;
    } else  // if x and y are the same object => alias
    {
      y += A.transpose() * x;
    }
  }  
}



void siconos::algebra::prod(const SiconosVector& x, const SiconosMatrix& A, BlockVector& y,
                            bool init) {

  if (A.size(0) != x.size()) THROW_EXCEPTION("inconsistent sizes between A and x.");

  if (A.size(1) != y.size()) THROW_EXCEPTION("inconsistent sizes between A and y.");

  unsigned int pos = 0;
  // For Each subvector of y, y[i], computes y[i] = transpose(subA) x, subA being a submatrix
  // of A corresponding to y[i] position. private_prod takes into account the fact that x and
  // y[i] may be block vectors.
  for (auto& it : y) {
    taxpy(x, A, pos, 0, it, init);
    pos += it->size();
  }
}

// ========== Products matrix - vector

siconos::algebra::SiconosVector siconos::algebra::prod(const SiconosMatrix& A,
                                                       const SiconosVector& x) {
  // To compute y = A * x
  if (A.size(1) != x.size()) THROW_EXCEPTION("inconsistent sizes between A and x.");
  
  return A*x;
}


const siconos::algebra::SimpleMatrix siconos::algebra::prod(const SiconosMatrix& A,
                                                            const SiconosMatrix& B) {

  SimpleMatrix C(A.size(0), B.size(1));
  prod(A, B, C);
  return C;
}


void siconos::algebra::prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C,
                            bool init) {
  // To compute C = A * B
  if ((A.size(1) != B.size(0))) THROW_EXCEPTION("inconsistent sizes between A and B");

  if (A.size(0) != C.size(0) || B.size(1) != C.size(1))
    THROW_EXCEPTION("inconsistent sizes between A and C or B and C.");

 // neither A or B is equal to identity or zero.
  if (init) {
    if (&C == &A || &C == &B)  // if common memory between A and C
    {
      C = A*B;      
    }
    else  // if no alias between C and A or B.
    {
      C.noalias() = A*B;
    }
  } else  // += case
  {
    if (&C == &A || &C == &B)  // if common memory between A and C
    {
      C += A*B;
    } else  // if no alias between C and A or B.
    {
      C.noalias() += A*B;
    }
  }
}

void siconos::algebra::prod(double a, const SiconosMatrix& A, const SiconosVector& x,
                            SiconosVector& y, bool init) {
  // To compute y = a*A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += a*A*x if init = false.

  if (A.size(1) != x.size()) THROW_EXCEPTION("inconsistent sizes between A and x.");

  if (A.size(0) != y.size()) THROW_EXCEPTION("inconsistent sizes between A and y.");

  if (init) {
    if (&x != &y)  // if no common memory between x and y.
    {
      y.noalias() = a * A * x;
    } else  // if x and y are the same object => alias
    {
      y = a * A * x;
    }
  } else  // += case
  {
    if (&x != &y)  // if no common memory between x and y.
    {
      y.noalias() += a * A * x;
    } else  // if x and y are the same object => alias
    {
      y += a * A * x;   
    }
  }  
}

void siconos::algebra::taxpy(const SiconosVector& x, const SiconosMatrix& A,
                             unsigned int startRow, unsigned int startCol,
                             std::shared_ptr<SiconosVector> y, bool init) {
  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between
  // el. of A of index (col) startCol and startCol + sizeY
  if (init)  // y = subA * x , else y += subA * x
    y->zero();

  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and
  // between columns startCol and (startCol+sizeX). Then computation of y = subA*x + y.
  auto sizeX = x.size();
  auto sizeY = y->size();

  assert(*y != x); // WARNING : compare addresses ?
  y->noalias() += A.block(startRow, startCol, sizeY, sizeX).transpose() * x;
}
