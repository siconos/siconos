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

#include "SiconosVector.hpp"

#include <boost/numeric/ublas/io.hpp>            // for >>
#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <boost/numeric/ublas/vector_sparse.hpp>

#include "BlockVector.hpp"
#include "SiconosConfig.h"
#include "SiconosException.hpp"
#include "SiconosMatrixVectorOp.hpp"  // For outer_prod decl. ...
#include "SiconosVectorOp.hpp"        // For setBlock decl. ...
#include "boost/numeric/bindings/blas.hpp"
#include "boost/numeric/bindings/std/vector.hpp"
#include "boost/numeric/bindings/ublas/matrix.hpp"
#include "boost/numeric/bindings/ublas/vector.hpp"
// #define DEBUG_MESSAGES
#include "SimpleMatrix.hpp"
#include "io.hpp"
#include "siconos_debug.h"

namespace ublas = boost::numeric::ublas;
namespace bindings_blas = boost::numeric::bindings::blas;

// namespace siconos::algebra {
// // // Do not document
// // /// @cond
// // #include "Question.hpp"
// struct IsDense : public siconos::tools::Question<bool> {
//   using siconos::tools::SiconosVisitor::visit;

//   void visit(const SiconosVector& v) { answer = v._dense; }

//   void visit(const BlockVector& v) { answer = false; }
// };

// struct IsSparse : public siconos::tools::Question<bool> {
//   using siconos::tools::SiconosVisitor::SiconosVisitor::visit;

//   void visit(const SiconosVector& v) { answer = !v._dense; }

//   void visit(const BlockVector& v) { answer = false; }
// };

// struct IsBlock : public siconos::tools::Question<bool> {
//   using siconos::tools::SiconosVisitor::SiconosVisitor::visit;

//   void visit(const SiconosVector& v) { answer = false; }

//   void visit(const BlockVector& v) { answer = true; }
// };

// /// @endcond
// }  // namespace siconos::algebra

// =================================================
//                CONSTRUCTORS
// =================================================

// Default
siconos::algebra::SiconosVector::SiconosVector() {
  _dense = true;
  vect.Dense = new DenseVect(ublas::zero_vector<double>());
}

// parameters: dimension and type.
siconos::algebra::SiconosVector::SiconosVector(unsigned row,
                                               siconos::algebra::UblasType type) {
  if (type == UblasType::SPARSE) {
    _dense = false;
    vect.Sparse = new SparseVect(ublas::zero_vector<double>(row));
  } else if (type == UblasType::DENSE) {
    _dense = true;
    vect.Dense = new DenseVect(ublas::zero_vector<double>(row));
  } else {
    THROW_EXCEPTION("invalid type");
  }
}

// parameters: dimension, default value for all components and type.
siconos::algebra::SiconosVector::SiconosVector(unsigned row, double val,
                                               siconos::algebra::UblasType type) {
  if (type == UblasType::SPARSE) {
    _dense = false;
    vect.Sparse = new SparseVect(row);
    fill(val);
  } else if (type == UblasType::DENSE) {
    _dense = true;
    vect.Dense = new DenseVect(ublas::scalar_vector<double>(row, val));
  } else {
    THROW_EXCEPTION("invalid type");
  }
}

// parameters: a vector (stl) of double and the type.
siconos::algebra::SiconosVector::SiconosVector(const std::vector<double>& v,
                                               siconos::algebra::UblasType typ) {
  if (typ != UblasType::DENSE) THROW_EXCEPTION("invalid type");

  _dense = true;
  vect.Dense = new DenseVect(v.size());
  std::copy(v.begin(), v.end(), (vect.Dense)->begin());
}

// Copy
siconos::algebra::SiconosVector::SiconosVector(const SiconosVector& svect)
    : std::enable_shared_from_this<SiconosVector>() {
  // if (siconos::tools::ask<IsDense>(svect))  // dense
  if (svect.num() == UblasType::DENSE) {
    _dense = true;
    vect.Dense = new DenseVect(svect.size());
    noalias(*vect.Dense) = (*svect.dense());
    // std::copy((vect.Dense)->begin(), (vect.Dense)->end(), (svect.dense())->begin());
  } else  // sparse
  {
    _dense = false;
    vect.Sparse = new SparseVect(svect.size());
    noalias(*vect.Sparse) = (*svect.sparse());
    // std::copy((vect.Sparse)->begin(), (vect.Sparse)->end(), (svect.sparse())->begin());
  }
  // Note FP: using constructor + noalias = (or std::copy) seems to be more
  // efficient than a call to ublas::vector copy constructor, this for
  // large or small vectors.
}

siconos::algebra::SiconosVector::SiconosVector(const DenseVect& m) {
  _dense = true;
  vect.Dense = new DenseVect(m.size());
  noalias(*vect.Dense) = m;
}

siconos::algebra::SiconosVector::SiconosVector(const SparseVect& m) {
  _dense = false;
  vect.Sparse = new SparseVect(m.size());
  noalias(*vect.Sparse) = m;
}

siconos::algebra::SiconosVector::SiconosVector(const std::string& file, bool ascii) {
  _dense = true;
  vect.Dense = new DenseVect();
  if (ascii) {
    siconos::algebra::io::read(file, *this, io::ASCII_IN);
  } else {
    siconos::algebra::io::read(file, *this, io::BINARY_IN);
  }
}

siconos::algebra::SiconosVector::SiconosVector(const SiconosVector& v1,
                                               const SiconosVector& v2) {
  unsigned int size1 = v1.size();
  auto num1 = v1.num();
  auto num2 = v2.num();
  //  if (siconos::tools::ask<IsDense>(v1) && siconos::tools::ask<IsDense>(v2)) {
  if (num1 == UblasType::DENSE && num2 == UblasType::DENSE) {
    _dense = true;
    vect.Dense = new DenseVect(size1 + v2.size());
  }
  //  else if (siconos::tools::ask<IsSparse>(v1) && siconos::tools::ask<IsSparse>(v2)) {
  else if (num1 == UblasType::SPARSE && num2 == UblasType::SPARSE) {
    _dense = false;
    vect.Sparse = new SparseVect(size1 + v2.size());
  } else {
    THROW_EXCEPTION("mixed dense and sparse vector detected");
  }
  setBlock(0, v1);
  setBlock(size1, v2);
}

// Copy a block vector into a SiconosVector
// This is mostly used to handle contiguous memory.
siconos::algebra::SiconosVector::SiconosVector(const BlockVector& input) {
  std::cout << "BLOCK COPY !!!!" << std::endl;
  if (input.isDense()) {
    _dense = true;
    vect.Dense = new DenseVect(input.size());
  } else {
    _dense = false;
    vect.Sparse = new SparseVect(input.size());
  }

  unsigned int pos = 0;
  for (auto it = input.begin(); it != input.end(); ++it) {
    setBlock(pos, **it);
    pos += (*it)->size();
  }
}

siconos::algebra::SiconosVector::~SiconosVector() noexcept {
  if (_dense)
    delete (vect.Dense);
  else
    delete (vect.Sparse);
}

// =================================================
//        get Ublas component (dense or sparse)
// =================================================

siconos::algebra::SparseVect* siconos::algebra::SiconosVector::sparse() const {
  if (_dense) THROW_EXCEPTION("the current vector is not a Sparse vector");

  return vect.Sparse;
}

double* siconos::algebra::SiconosVector::getArray() const {
  assert(
      vect.Dense &&
      "siconos::algebra::SiconosVector::getArray() : not yet implemented for sparse vector.");

  return &(((*vect.Dense).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void siconos::algebra::SiconosVector::zero() {
  if (_dense)
    bindings_blas::scal(0.0, *vect.Dense);

  else {
    assert(vect.Sparse);
    *vect.Sparse *= 0.0;
  }
}

void siconos::algebra::SiconosVector::fill(double value) {
  if (!_dense) {
    for (unsigned int i = 0; i < (vect.Sparse)->size(); ++i)
      (vect.Sparse)->push_back(i, value);
  } else
    bindings_blas::set(value, *vect.Dense);
}

//=======================
// set vector dimension
//=======================

void siconos::algebra::SiconosVector::resize(unsigned int n, bool preserve) {
  if (_dense)
    (vect.Dense)->resize(n, preserve);
  else
    (vect.Sparse)->resize(n, preserve);
}

//=======================
//       get norm
//=======================

double siconos::algebra::SiconosVector::normInf() const {
  if (_dense)
    return norm_inf(*vect.Dense);
  else  // if(num==UblasType::SPARSE)
    return norm_inf(*vect.Sparse);
}

double siconos::algebra::SiconosVector::norm2() const {
  if (_dense)
    return ublas::norm_2(*vect.Dense);
  else  // if(num==UblasType::SPARSE)
    return ublas::norm_2(*vect.Sparse);
}
//======================================
// get sum of all elements of the vector
//=====================================
double siconos::algebra::SiconosVector::vector_sum() const {
  if (_dense)
    return ublas::sum(*vect.Dense);
  else
    return ublas::sum(*vect.Sparse);
}

//=====================
// screen display
//=====================

void siconos::algebra::SiconosVector::display() const {
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);
  if (_dense)
    std::cout << *vect.Dense << std::endl;
  else if (vect.Sparse)
    std::cout << *vect.Sparse << std::endl;
}

//=====================
// convert to an ostream
//=====================

std::ostream& siconos::algebra::operator<<(std::ostream& os, const SiconosVector& sv) {
  if (sv._dense)
    os << *sv.vect.Dense;
  else
    os << *sv.vect.Sparse;
  return os;
}

//=============================
// Elements access (get or set)
//=============================

double siconos::algebra::SiconosVector::getValue(unsigned int row) const {
  assert(row < size() &&
         "siconos::algebra::SiconosVector::getValue(index) : Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row);
}

void siconos::algebra::SiconosVector::setValue(unsigned int row, double value) {
  assert(row < size() &&
         "siconos::algebra::SiconosVector::setValue(index, value) : Index out of range");
  if (_dense)
    (*vect.Dense)(row) = value;
  else
    (*vect.Sparse)(row) = value;
}

double& siconos::algebra::SiconosVector::operator()(unsigned int row) {
  assert(row < size() &&
         "siconos::algebra::SiconosVector::operator ( index ): Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

double siconos::algebra::SiconosVector::operator()(unsigned int row) const {
  assert(row < size() &&
         "siconos::algebra::SiconosVector::operator ( index ): Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return ((*vect.Sparse)(row)).ref();
}

//============================================
// Access (get or set) to blocks of elements
//============================================

void siconos::algebra::SiconosVector::setBlock(unsigned int index, const SiconosVector& vIn) {
  // Set current vector elements, starting from position "index", to the values of vector vIn

  // Exceptions ...
  assert(&vIn != this &&
         "siconos::algebra::SiconosVector::this->setBlock(pos,vIn): vIn = this.");

  assert(index < size() && "siconos::algebra::SiconosVector::setBlock : invalid ranges");

  auto end = vIn.size() + index;
  assert(end <= size() && "siconos::algebra::SiconosVector::setBlock : invalid ranges");

  assert(vIn.num() == num() &&
         "siconos::algebra::SiconosVector::setBlock: inconsistent types.");

  if (_dense)
    noalias(ublas::subrange(*vect.Dense, index, end)) = *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, end)) = *vIn.sparse();
}

void siconos::algebra::SiconosVector::toBlock(SiconosVector& vOut, unsigned int sizeB,
                                              unsigned int startIn,
                                              unsigned int startOut) const {
  // To copy a subBlock of the vector (from position startIn to startIn+sizeB) into vOut (from
  // pos. startOut to startOut+sizeB). Check dim ...
  assert(startIn < size() &&
         "vector toBlock(v1,v2,...): start position in input vector is out of range.");

  assert(startOut < vOut.size() &&
         "vector toBlock(v1,v2,...): start position in output vector is out of range.");

  assert(startIn + sizeB <= size() &&
         "vector toBlock(v1,v2,...): end position in input vector is out of range.");
  assert(startOut + sizeB <= vOut.size() &&
         "vector toBlock(v1,v2,...): end position in output vector is out of range.");

  unsigned int endOut = startOut + sizeB;
  auto numIn = num();
  auto numOut = vOut.num();

  if (numIn == numOut) {
    if (numIn == UblasType::DENSE)
      noalias(ublas::subrange(*vOut.dense(), startOut, endOut)) =
          ublas::subrange(*vect.Dense, startIn, startIn + sizeB);
    else  // if(numIn == UblasType::SPARSE)// vIn / vOut are Sparse
      noalias(ublas::subrange(*vOut.sparse(), startOut, endOut)) =
          ublas::subrange(*vect.Sparse, startIn, startIn + sizeB);
  } else  // vIn and vout of different types ...
  {
    if (numIn == UblasType::DENSE)  // vIn Dense
      noalias(ublas::subrange(*vOut.sparse(), startOut, endOut)) =
          ublas::subrange(*vect.Dense, startIn, startIn + sizeB);
    else  // if(numIn == UblasType::SPARSE)// vIn Sparse
      noalias(ublas::subrange(*vOut.dense(), startOut, endOut)) =
          ublas::subrange(*vect.Sparse, startIn, startIn + sizeB);
  }
}

void siconos::algebra::SiconosVector::addBlock(unsigned int index, const SiconosVector& vIn) {
  // Add vIn to the current vector, starting from position "index".

  if (&vIn == this) THROW_EXCEPTION("try to add a vector to itself.");

  unsigned int end = vIn.size();
  if ((index + end) > size()) THROW_EXCEPTION("invalid ranges");

  auto numVin = vIn.num();

  if (numVin != num()) THROW_EXCEPTION("inconsistent types.");

  if (_dense)
    noalias(ublas::subrange(*vect.Dense, index, index + end)) += *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, index + end)) += *vIn.sparse();
}

void siconos::algebra::SiconosVector::subBlock(unsigned int index, const SiconosVector& vIn) {
  // Add vIn from the current vector, starting from position "index".

  unsigned int end = vIn.size();
  if ((index + end) > size()) THROW_EXCEPTION("invalid ranges");

  auto numVin = vIn.num();
  if (numVin != num()) THROW_EXCEPTION("inconsistent types.");

  if (_dense)
    noalias(ublas::subrange(*vect.Dense, index, index + end)) -= *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, index + end)) -= *vIn.sparse();
}

//===============
//  Assignment
//===============

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator=(
    const SiconosVector& vIn) {
  if (&vIn == this) return *this;  // auto-assignment.

  assert(size() == vIn.size() &&
         "siconos::algebra::SiconosVector::operator = failed: inconsistent sizes.");

  auto vInNum = vIn.num();
  {
    switch (num()) {
      case UblasType::DENSE:
        switch (vInNum) {
          case UblasType::DENSE:
            noalias(*vect.Dense) = *vIn.dense();
            break;
          case UblasType::SPARSE:
            noalias(*vect.Dense) = *vIn.sparse();
            break;
          default:
            THROW_EXCEPTION("invalid type");
            break;
        }
        break;
      case UblasType::SPARSE:
        if (vInNum == UblasType::DENSE)
          noalias(*vect.Sparse) = *vIn.dense();
        else if (vInNum == UblasType::SPARSE)
          noalias(*vect.Sparse) = *vIn.sparse();
        else
          THROW_EXCEPTION("invalid type");
        break;
      default:
        THROW_EXCEPTION("invalid type");
        break;
    }
  }
  return *this;
}

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator=(
    const BlockVector& vIn) {
  unsigned int pos = 0;
  for (auto it = vIn.begin(); it != vIn.end(); ++it) {
    setBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator=(
    const DenseVect& d) {
  if (!_dense) THROW_EXCEPTION("the current vector is not dense.");
  if (d.size() != size()) THROW_EXCEPTION("inconsistent size.");

  bindings_blas::copy(d, *vect.Dense);
  return *this;
}

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator=(
    const SparseVect& sp) {
  if (_dense) THROW_EXCEPTION("current vector is not sparse.");
  if (sp.size() != size()) THROW_EXCEPTION("inconsistent size.");

  noalias(*vect.Sparse) = sp;

  return *this;
}

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator=(const double* d) {
  assert(_dense &&
         "siconos::algebra::SiconosVector::operator = double* : forbidden: the current vector "
         "is not dense.");

  bindings_blas::detail::copy(vect.Dense->size(), d, 1, getArray(), 1);
  return *this;
}

unsigned siconos::algebra::SiconosVector::copyData(double* data) const {
  assert(_dense &&
         "siconos::algebra::SiconosVector::copyData : forbidden: the current vector is not "
         "dense.");

  unsigned size = vect.Dense->size();
  bindings_blas::detail::copy(vect.Dense->size(), getArray(), 1, data, 1);
  return size;
}

//=================================
// Op. and assignment (+=, -= ... )
//=================================

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator+=(
    const SiconosVector& vIn) {
  if (&vIn == this)  // alias
  {
    // Note: using this *= 2.0 is much more time-consuming.
    switch (num()) {
      case UblasType::DENSE:
        *vect.Dense += *vect.Dense;
        break;
      case UblasType::SPARSE:
        *vect.Sparse += *vect.Sparse;
        break;
      default:
        THROW_EXCEPTION("invalid type.");
        break;
    }
    return *this;
  }

  auto vInNum = vIn.num();
  {
    switch (num()) {
      case UblasType::DENSE:
        switch (vInNum) {
          case UblasType::DENSE:
            noalias(*vect.Dense) += *vIn.dense();
            break;
          case UblasType::SPARSE:
            noalias(*vect.Dense) += *vIn.sparse();
            break;
          default:
            THROW_EXCEPTION("invalid type");
            break;
        }
        break;
      case UblasType::SPARSE:
        if (vInNum == UblasType::SPARSE)
          noalias(*vect.Sparse) += *vIn.sparse();
        else
          THROW_EXCEPTION("can not add a dense to a sparse.");
        break;
      default:
        THROW_EXCEPTION("invalid type.");
        break;
    }
  }
  return *this;
}
siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator+=(
    const BlockVector& vIn) {
  unsigned int pos = 0;
  for (auto it : vIn) {
    addBlock(pos, *it);
    pos += it->size();
  }
  return *this;
}

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator-=(
    const SiconosVector& vIn) {
  if (&vIn == this) {
    this->zero();
    return *this;
  }

  auto vInNum = vIn.num();
  {
    switch (num()) {
      case UblasType::DENSE:
        switch (vInNum) {
          case UblasType::DENSE:
            noalias(*vect.Dense) -= *vIn.dense();
            break;
          case UblasType::SPARSE:
            noalias(*vect.Dense) -= *vIn.sparse();
            break;
          default:
            THROW_EXCEPTION("invalid type.");
            break;
        }
        break;
      case UblasType::SPARSE:
        if (vInNum == UblasType::SPARSE)
          noalias(*vect.Sparse) -= *vIn.sparse();
        else
          THROW_EXCEPTION("can not sub a dense to a sparse.");
        break;
      default:
        THROW_EXCEPTION("invalid type.");
        break;
    }
  }
  return *this;
}

siconos::algebra::SiconosVector& siconos::algebra::SiconosVector::operator-=(
    const BlockVector& vIn) {
  unsigned int pos = 0;
  for (auto it = vIn.begin(); it != vIn.end(); ++it) {
    subBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}

//===============
// Comparison
//===============

bool siconos::algebra::operator==(const SiconosVector& m, const SiconosVector& x) {
  DEBUG_PRINTF("norm = %12.8e \n", (m - x).normInf());
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n",
               std::numeric_limits<double>::epsilon());
  DEBUG_EXPR(std::cout << std::boolalpha
                       << ((m - x).normInf() <= std::numeric_limits<double>::epsilon())
                       << std::endl;);
  double atol = 1e-14;
  double rtol = std::numeric_limits<double>::epsilon();
  return ((m - x).normInf() <= atol + rtol * x.normInf());
}

//==================
// y = scalar * x
//==================

siconos::algebra::SiconosVector siconos::algebra::operator*(const SiconosVector& m, double d) {
  auto numM = m.num();

  if (numM == UblasType::DENSE) {
    // Copy m into p and call bindings_blas::scal(d,p), p = d*p.
    DenseVect p = *m.dense();
    bindings_blas::scal(d, p);
    return p;
  } else  // if(numM==UblasType::SPARSE)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

siconos::algebra::SiconosVector siconos::algebra::operator*(double d, const SiconosVector& m) {
  auto numM = m.num();

  if (numM == UblasType::DENSE) {
    // Copy m into p and call bindings_blas::scal(d,p), p = d*p.
    DenseVect p = *m.dense();
    bindings_blas::scal(d, p);
    return p;
  } else  // if(numM==UblasType::SPARSE)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

siconos::algebra::SiconosVector siconos::algebra::operator/(const SiconosVector& m, double d) {
  auto numM = m.num();

  if (numM == UblasType::DENSE) {
    DenseVect p = *m.dense();
    bindings_blas::scal((1.0 / d), p);
    return p;
  }

  else  // if(numM==UblasType::SPARSE){
    return (SparseVect)(*m.sparse() / d);
}

//====================
//  Vectors addition
//====================

siconos::algebra::SiconosVector siconos::algebra::operator+(const SiconosVector& x,
                                                            const SiconosVector& y) {
  if (x.size() != y.size()) THROW_EXCEPTION("inconsistent sizes");

  auto numX = x.num();
  auto numY = y.num();

  if (numX == numY)  // x, y SiconosVector of the same type
  {
    if (numX == UblasType::DENSE) {
      return (DenseVect)(*x.dense() + *y.dense());
    } else
      return (SparseVect)(*x.sparse() + *y.sparse());
  }

  else  // x, y SiconosVector with y and x of different types
  {
    if (numX == UblasType::DENSE)
      return (DenseVect)(*x.dense() + *y.sparse());
    else
      return (DenseVect)(*x.sparse() + *y.dense());
  }
}

void siconos::algebra::add(const SiconosVector& x, const SiconosVector& y, SiconosVector& z) {
  // Computes z = x + y in an "optimized" way (in comparison with operator +)

  if (x.size() != y.size() || x.size() != z.size()) THROW_EXCEPTION("inconsistent sizes");

  auto numX = x.num();
  auto numY = y.num();
  auto numZ = z.num();

  if (&z == &x)  // x, and z are the same object.
  {
    z += y;
  } else if (&z == &y)  // y and z are the same object, different from x
  {
    z += x;
  } else  // No common memory between x,y and z
  {
    if (numZ != UblasType::BLOCK)  // z is a SiconosVector
    {
      if (numX == numY && numX != UblasType::BLOCK)  // x, y SiconosVector of the same type
      {
        if (numX == UblasType::DENSE) {
          if (numZ != UblasType::DENSE)
            THROW_EXCEPTION("Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() + *y.dense();
        } else {
          if (numZ == UblasType::DENSE)
            noalias(*z.dense()) = *x.sparse() + *y.sparse();
          else
            noalias(*z.sparse()) = *x.sparse() + *y.sparse();
        }
      } else if (numX != UblasType::BLOCK &&
                 numY != UblasType::BLOCK)  // x and y of different types => z must be dense.
      {
        if (numZ != UblasType::DENSE) THROW_EXCEPTION("z can not be sparse.");
        if (numX == UblasType::DENSE)
          noalias(*z.dense()) = *x.dense() + *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() + *y.dense();
      }
    }
  }
}

//======================
//  Vectors subtraction
//======================

siconos::algebra::SiconosVector siconos::algebra::operator-(const SiconosVector& x,
                                                            const SiconosVector& y) {
  if (x.size() != y.size()) THROW_EXCEPTION("inconsistent sizes");

  auto numX = x.num();
  auto numY = y.num();

  if (numX == numY)  // x, y SiconosVector of the same type
  {
    if (numX == UblasType::DENSE) {
      return (DenseVect)(*x.dense() - *y.dense());
    } else
      return (SparseVect)(*x.sparse() - *y.sparse());
  } else  // x, y SiconosVector with y and x of different types
  {
    if (numX == UblasType::DENSE)
      return (DenseVect)(*x.dense() - *y.sparse());
    else
      return (DenseVect)(*x.sparse() - *y.dense());
  }
}

void siconos::algebra::sub(const SiconosVector& x, const SiconosVector& y, SiconosVector& z) {
  // Computes z = x - y in an "optimized" way (in comparison with operator +)

  if (x.size() != y.size() || x.size() != z.size()) THROW_EXCEPTION("inconsistent sizes");

  auto numX = x.num();
  auto numY = y.num();
  auto numZ = z.num();

  if (&z == &x)  // x and z are the same object.
  {
    z -= y;
  } else if (&z == &y)  // y and z are the same object
  {
    {
      if (numX == UblasType::DENSE) {
        if (numZ != UblasType::DENSE)
          THROW_EXCEPTION("Subtraction of two dense vectors into a sparse.");
        *z.dense() = *x.dense() - *y.dense();
      } else {
        if (numZ == UblasType::DENSE)
          *z.dense() = *x.sparse() - *y.dense();
        else
          *z.sparse() = *x.sparse() - *y.sparse();
      }
    }
  } else  // No common memory between x or y and z
  {
    if (numZ != UblasType::BLOCK)  // z is a SiconosVector
    {
      if (numX == numY && numX != UblasType::BLOCK)  // x, y SiconosVector of the same type
      {
        if (numX == UblasType::DENSE) {
          if (numZ != UblasType::DENSE)
            THROW_EXCEPTION("Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() - *y.dense();
        } else {
          if (numZ == UblasType::DENSE)
            noalias(*z.dense()) = *x.sparse() - *y.sparse();
          else
            noalias(*z.sparse()) = *x.sparse() - *y.sparse();
        }
      } else if (numX != UblasType::BLOCK &&
                 numY != UblasType::BLOCK)  // x and y of different types => z must be dense.
      {
        if (numZ != UblasType::DENSE) THROW_EXCEPTION("z can not be sparse.");
        if (numX == UblasType::DENSE)
          noalias(*z.dense()) = *x.dense() - *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() - *y.dense();
      }
    }
  }
}

void siconos::algebra::axpby(double a, const SiconosVector& x, double b, SiconosVector& y) {
  // Computes y = ax + by

  if (x.size() != y.size()) THROW_EXCEPTION("inconsistent sizes");

  auto numX = x.num();
  auto numY = y.num();

  if (numX == numY)  // x and y of the same type
  {
    if (numX == UblasType::DENSE)  // all dense
    {
      bindings_blas::scal(b, *y.dense());
      bindings_blas::axpy(a, *x.dense(), *y.dense());
    } else  // all sparse
    {
      *y.sparse() *= b;
      if (&y != &x)
        noalias(*y.sparse()) += a * *x.sparse();
      else
        *y.sparse() += a * *x.sparse();
    }
  }

  else  // x and y of different types
  {
    y *= b;
    {
      if (numX == UblasType::DENSE)
        *y.sparse() += a * *x.dense();
      else
        *y.dense() += a * *x.sparse();
    }
  }
}

void siconos::algebra::axpy(double a, const SiconosVector& x, SiconosVector& y) {
  // Computes y = ax + y

  if (x.size() != y.size()) THROW_EXCEPTION("nconsistent sizes");

  auto numX = x.num();
  auto numY = y.num();

  if (numX == numY)  // x and y of the same type
  {
    if (numX == UblasType::DENSE)  // all dense
      bindings_blas::axpy(a, *x.dense(), *y.dense());

    else  // all sparse
    {
      if (&y != &x)
        noalias(*y.sparse()) += a * *x.sparse();
      else
        *y.sparse() += a * *x.sparse();
    }
  }

  else  // x and y of different types
  {
    {
      if (numX == UblasType::DENSE)
        *y.sparse() += a * *x.dense();
      else
        *y.dense() += a * *x.sparse();
    }
  }
}

double siconos::algebra::inner_prod(const SiconosVector& x, const SiconosVector& m) {
  if (x.size() != m.size()) THROW_EXCEPTION("inconsistent sizes");

  auto numM = m.num();
  auto numX = x.num();

  if (numX == numM) {
    if (numM == UblasType::DENSE)
      return bindings_blas::dot(*x.dense(), *m.dense());
    else
      return ublas::inner_prod(*x.sparse(), *m.sparse());
  } else if (numM == UblasType::DENSE)
    return ublas::inner_prod(*x.sparse(), *m.dense());
  else
    return ublas::inner_prod(*x.dense(), *m.sparse());
}

// outer_prod(v,w) = trans(v)*w
siconos::algebra::SimpleMatrix siconos::algebra::outer_prod(const SiconosVector& x,
                                                            const SiconosVector& m) {
  auto numM = m.num();
  auto numX = x.num();

  if (numM == UblasType::DENSE) {
    if (numX == UblasType::DENSE)
      return (DenseMat)(ublas::outer_prod(*x.dense(), *m.dense()));

    else  // if(numX == UblasType::SPARSE)
      return (DenseMat)(ublas::outer_prod(*x.sparse(), *m.dense()));
  } else  // if(numM == UblasType::SPARSE)
  {
    if (numX == UblasType::DENSE)
      return (DenseMat)(ublas::outer_prod(*x.dense(), *m.sparse()));

    else  // if(numX == UblasType::SPARSE)
      return (DenseMat)(ublas::outer_prod(*x.sparse(), *m.sparse()));
  }
}

void siconos::algebra::scal(double a, const SiconosVector& x, SiconosVector& y, bool init) {
  // To compute y = a *x (init = true) or y += a*x (init = false)

  if (&x == &y) {
    if (init)
      y *= a;
    else {
      y *= (1.0 + a);
    }
  } else {
    unsigned int sizeX = x.size();
    unsigned int sizeY = y.size();

    if (sizeX != sizeY) THROW_EXCEPTION("sizes are not consistent.");

    auto numY = y.num();
    auto numX = x.num();
    if (numX == numY) {
      if (numX == UblasType::DENSE)  // ie if both are Dense
      {
        if (init)
          noalias(*y.dense()) = a * *x.dense();
        else
          noalias(*y.dense()) += a * *x.dense();
      } else  // if both are sparse
      {
        if (init)
          noalias(*y.sparse()) = a * *x.sparse();
        else
          noalias(*y.sparse()) += a * *x.sparse();
      }
    } else {
      if (numY == UblasType::BLOCK || numX == UblasType::BLOCK)  // if y or x is block
      {
        if (init) {
          y = x;
          y *= a;
        } else {
          SiconosVector tmp(x);
          tmp *= a;
          y += tmp;
        }
      } else {
        if (numY == UblasType::DENSE)  // if y is dense
        {
          if (init)
            noalias(*y.dense()) = a * *x.sparse();
          else
            noalias(*y.dense()) += a * *x.sparse();
        } else
          THROW_EXCEPTION("Operation not allowed on non-dense vector.");
      }
    }
  }
}

void siconos::algebra::subscal(double a, const SiconosVector& x, SiconosVector& y,
                               const std::vector<std::size_t>& coord, bool init) {
  // To compute sub_y = a *sub_x (init = true) or sub_y += a*sub_x (init = false)
  // Coord  = [r0x r1x r0y r1y];
  // subX is the sub-vector of x, for row numbers between r0x and r1x-1.
  // The same for y with riy.

  // Check dimensions
  auto dimX = coord[1] - coord[0];
  auto dimY = coord[3] - coord[2];
  if (dimY != dimX) THROW_EXCEPTION("inconsistent sizes between (sub)x and (sub)y.");
  if (dimY > y.size() || dimX > x.size()) THROW_EXCEPTION("input index too large.");

  auto numY = y.num();
  auto numX = x.num();

  if (&x == &y)  // if x and y are the same object
  {
    if (numX == UblasType::DENSE)  // Dense
    {
      ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2]) {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      } else {
        ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    } else  // if (numX == UblasType::SPARSE) // Sparse
    {
      ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2]) {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      } else {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
  } else {
    if (numX == numY) {
      if (numX == UblasType::DENSE)  // ie if both are Dense
      {
        ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      } else  // if both are sparse
      {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
    } else  // x and y of different types ...
    {
      if (numY == UblasType::DENSE)  // y dense, x sparse
      {
        ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      } else  // y sparse, x dense => fails
        THROW_EXCEPTION("Operation not allowed (y sparse, x dense).");
    }
  }
}
void siconos::algebra::cross_product(const SiconosVector& V1, const SiconosVector& V2,
                                     SiconosVector& VOUT) {
  if (V1.size() != 3 || V2.size() != 3 || VOUT.size() != 3)
    THROW_EXCEPTION("allowed only with dim 3.");

  double aux = V1.getValue(1) * V2.getValue(2) - V1.getValue(2) * V2.getValue(1);
  VOUT.setValue(0, aux);

  aux = V1.getValue(2) * V2.getValue(0) - V1.getValue(0) * V2.getValue(2);
  VOUT.setValue(1, aux);

  aux = V1.getValue(0) * V2.getValue(1) - V1.getValue(1) * V2.getValue(0);
  VOUT.setValue(2, aux);
}

//

void siconos::algebra::abs_wise(const SiconosVector& V, SiconosVector& Vabs) {
  for (unsigned int it = 0; it < V.size(); ++it) {
    Vabs.setValue(it, std::abs(V.getValue(it)));
  };
}

//

void siconos::algebra::getMax(const SiconosVector& V, double& maxvalue, unsigned int& idmax) {
  maxvalue = V.getValue(0);
  idmax = 0;
  for (unsigned int it = 1; it < V.size(); ++it) {
    if (V.getValue(it) > maxvalue) {
      maxvalue = V.getValue(it);
      idmax = it;
    };
  };
}

//

void siconos::algebra::getMin(const SiconosVector& V, double& minvalue, unsigned int& idmin) {
  minvalue = V.getValue(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it) {
    if (V.getValue(it) < minvalue) {
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}

// struct exp_op
// {
//   double operator()(double d) const
//   {
//     return std::exp(d);
//   }
// };

void siconos::algebra::setBlock(const SiconosVector& vIn, std::shared_ptr<SiconosVector> vOut,
                                unsigned int sizeB, unsigned int startIn,
                                unsigned int startOut) {
  unsigned int endOut = startOut + sizeB;
  auto numIn = vIn.num();
  auto numOut = vOut->num();
  assert(vOut->size() >= endOut && "The output vector is too small");
  if (numIn == numOut) {
    if (numIn == UblasType::DENSE)  // vIn / vOut are Dense
      noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) =
          ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
    else  // if(numIn == UblasType::SPARSE)// vIn / vOut are Sparse
      noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) =
          ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
  } else  // vIn and vout of different types ...
  {
    if (numIn == UblasType::DENSE)  // vIn Dense
      noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) =
          ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
    else  // if(numIn == UblasType::SPARSE)// vIn Sparse
      noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) =
          ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
  }
}

unsigned int siconos::algebra::SiconosVector::size(void) const {
  if (!_dense) {
    return (vect.Sparse->size());
  } else {
    return (vect.Dense->size());
  }
}

siconos::algebra::SiconosVector& siconos::algebra::operator*=(SiconosVector& v,
                                                              const double& s) {
  if (v._dense)
    *v.dense() *= s;
  else
    *v.sparse() *= s;
  return v;
}

siconos::algebra::SiconosVector& siconos::algebra::operator/=(SiconosVector& v,
                                                              const double& s) {
  if (v._dense)
    *v.dense() /= s;
  else
    *v.sparse() /= s;
  return v;
}

siconos::algebra::SiconosVector::iterator siconos::algebra::SiconosVector::begin() {
  return siconos::algebra::SiconosVector::iterator(*this, 0);
}

siconos::algebra::SiconosVector::const_iterator siconos::algebra::SiconosVector::begin()
    const {
  return siconos::algebra::SiconosVector::const_iterator(*this, 0);
}

siconos::algebra::SiconosVector::iterator siconos::algebra::SiconosVector::end() {
  return siconos::algebra::SiconosVector::iterator(*this, size());
}

siconos::algebra::SiconosVector::const_iterator siconos::algebra::SiconosVector::end() const {
  return siconos::algebra::SiconosVector::const_iterator(*this, size());
}

siconos::algebra::SiconosVector::operator std::vector<double>() {
  std::vector<double> v(size());
  std::copy(begin(), end(), v.begin());
  return v;
}
