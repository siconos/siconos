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
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>

#include <cstddef>
#include <iostream>
#include <vector>

#include "SiconosAlgebraTools.hpp"  // for isComparableTo
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for isComparableTo ...
#include "siconos_debug.h"
// //#define DEBUG_STDOUT
// //#define DEBUG_MESSAGES
//


// =================================================
//                CONSTRUCTORS
// =================================================
siconos::algebra::BlockVector::BlockVector() {
  _tabIndex = std::make_shared<std::vector<std::size_t>>();
}

siconos::algebra::BlockVector::BlockVector(const BlockVector& v) {
  auto nbBlocks = v.numberOfBlocks();
  _tabIndex = std::make_shared<std::vector<std::size_t>>();
  _tabIndex->reserve(nbBlocks);
  _vect.reserve(nbBlocks);
  for (auto& it : v) {
    _vect.push_back(std::make_shared<SiconosVector>(*it));
    _sizeV += (*it).size();
    _tabIndex->push_back(_sizeV);
  }
}

siconos::algebra::BlockVector::BlockVector(std::shared_ptr<SiconosVector> v1,
                                           std::shared_ptr<SiconosVector> v2) {
  // Insert the two vectors in the container
  // NO COPY !!
  if (!v1 && !v2) THROW_EXCEPTION("both vectors are nullptr.");

  _tabIndex = std::make_shared<std::vector<std::size_t>>();

  _tabIndex->reserve(2);
  _vect.reserve(2);

  if (v1) {
    _vect.push_back(v1);
    _sizeV = v1->size();
    _tabIndex->push_back(_sizeV);
  } else
  // If first parameter is a nullptr pointer, then set this(1) to a SiconosVector of the same
  // size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    _sizeV = v2->size();

    _vect.push_back(std::make_shared<SiconosVector>(_sizeV));
    _tabIndex->push_back(_sizeV);
  }
  if (v2) {
    _vect.push_back(v2);
    _sizeV += v2->size();
    _tabIndex->push_back(_sizeV);
  } else  // If second parameter is a nullptr pointer, then set this(2) to a SiconosVector of
          // the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.

    _vect.push_back(std::make_shared<SiconosVector>(v1->size()));
    _sizeV += v1->size();
    _tabIndex->push_back(_sizeV);
  }
}

siconos::algebra::BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int dim) {
  _tabIndex = std::make_shared<std::vector<std::size_t>>();
  _tabIndex->reserve(numberOfBlocks);
  _vect.reserve(numberOfBlocks);
  for (unsigned int i = 0; i < numberOfBlocks; ++i) {
    _vect.push_back(std::make_shared<SiconosVector>(dim));
    _tabIndex->push_back(dim * (i + 1));
  }
  _sizeV = dim * numberOfBlocks;
}

siconos::algebra::BlockVector::BlockVector(unsigned int numberOfBlocks) {
  _tabIndex = std::make_shared<std::vector<std::size_t>>();
  _tabIndex->resize(numberOfBlocks);
  _vect.resize(numberOfBlocks);
}

// ===========================
//      private method
// ===========================

void siconos::algebra::BlockVector::_update() {
  _sizeV = 0;
  _tabIndex = std::make_shared<std::vector<std::size_t>>();

  for (auto it : _vect) {
    if (it) {
      _sizeV += it->size();
    }
    _tabIndex->push_back(_sizeV);
  }
}

// ===========================
//       fill vector
// ===========================

// bool siconos::algebra::BlockVector::isDense() const {
//   return std::find_if(_vect.begin(), _vect.end(), TestDense()) != _vect.end();
// }

void siconos::algebra::BlockVector::zero() {
  for (auto& it : _vect) it->zero();
}

void siconos::algebra::BlockVector::fill(double value) {
  for (auto& it : _vect) {
    if (it) {
      it->fill(value);
    }
  }
}

//=====================
// screen display
//=====================

void siconos::algebra::BlockVector::display() const {
  std::cout << "=======> Block Vector Display (" << _tabIndex->size()
            << " block(s)): " << std::endl;
  for (auto& it : _vect) {
    DEBUG_EXPR(std::cout << "(*it)" << (*it) << std::endl;);
    if (it)
      it->display();
    else
      std::cout << "(*it)-> nullptr" << std::endl;
  }
}

//=====================
// convert to an ostream
//=====================

std::ostream& siconos::algebra::operator<<(std::ostream& os, const BlockVector& bv) {
  os << "[" << bv._vect.size() << "](";
  for (auto& it : bv._vect) {
    if (it)
      os << *it;
    else
      os << "(nil)";
    os << " |";  // separator
  }
  os << ")";
  return os;
}

//=============================
// Elements access (get or set)
//=============================

double siconos::algebra::BlockVector::getValue(unsigned int pos) const {
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size()) blockNum++;

  unsigned int relativePos = pos;

  if (blockNum != 0) relativePos -= (*_tabIndex)[blockNum - 1];

  return (*_vect[blockNum])(relativePos);
}

void siconos::algebra::BlockVector::setValue(unsigned int pos, double value) {
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size()) blockNum++;

  unsigned int relativePos = pos;

  if (blockNum != 0) relativePos -= (*_tabIndex)[blockNum - 1];

  (*_vect[blockNum])(relativePos) = value;
}

double& siconos::algebra::BlockVector::operator()(unsigned int pos) {
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size()) blockNum++;

  unsigned int relativePos = pos;

  if (blockNum != 0) relativePos -= (*_tabIndex)[blockNum - 1];

  return (*_vect[blockNum])(relativePos);
}

double siconos::algebra::BlockVector::operator()(unsigned int pos) const {
  return getValue(pos);
}

//============================================
// Access (get or set) to blocks of elements
//============================================

void siconos::algebra::BlockVector::setVector(unsigned int pos, const SiconosVector& v) {
  assert(pos < _vect.size() && "insertion out of vector size");
  if (!_vect[pos]) THROW_EXCEPTION("this[pos] == nullptr pointer.");
  *_vect[pos] = v;
}

void siconos::algebra::BlockVector::setVectorPtr(unsigned int pos,
                                                 std::shared_ptr<SiconosVector> v) {
  assert(pos < _vect.size() && "insertion out of vector size");
  _vect[pos] = v;
  _update();
}

void siconos::algebra::BlockVector::setAllVect(VectorOfVectors& v) {
  _vect = v;
  _update();
}

unsigned int siconos::algebra::BlockVector::getNumVectorAtPos(unsigned int pos) const {
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size() - 1) blockNum++;
  return blockNum;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator=(
    const BlockVector& vIn) {
  if (&vIn == this)
    return *this;
  else {
    if (siconos::algebra::isComparableTo(*this,
                                         vIn))  // if vIn and this are "block-consistent"
    {
      auto it2 = vIn.begin();
      for (auto& it1 : _vect) {
        (*it1) = (**it2);
        it2++;
      }
    } else {
      for (unsigned int i = 0; i < _sizeV; ++i) (*this)(i) = vIn(i);
    }
    return *this;
  }
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator=(const double* data) {
  unsigned indxPos = 0;

  for (auto& vect : _vect) {
    double* ptr = const_cast<double*>(&data[indxPos]);  // const_cat needed by Eigne::Map
    *vect = Eigen::Map<SiconosVector>(ptr, vect->size());
    indxPos += vect->size();
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator-=(
    const BlockVector& vIn) {
  if (siconos::algebra::isComparableTo(*this, vIn))  // if vIn and this are "block-consistent"
  {
    unsigned int i = 0;
    for (auto& it1 : _vect) {
      *it1 -= *(vIn.vector(i++));
    }
  } else  // use of a temporary SimpleVector... bad way, to be improved. But this case happens
          // rarely ...
  {
    for (unsigned int i = 0; i < _sizeV; ++i) (*this)(i) -= vIn(i);
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator-=(
    const SiconosVector& vIn) {
  unsigned int dim = vIn.size();  // size of the block to be added.
  if (dim > _sizeV) THROW_EXCEPTION("invalid ranges");

  unsigned int currentSize;
  unsigned int index = 0;
  for (auto& v : _vect) {
    currentSize = v->size();
    *v -= vIn.segment(index, currentSize);
    index += currentSize;
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator+=(
    const BlockVector& vIn) {
  if (siconos::algebra::isComparableTo(*this, vIn))  // if vIn and this are "block-consistent"
  {
    unsigned int i = 0;

    for (auto& it1 : _vect) {
      *it1 += *(vIn.vector(i++));
    }
  } else  // use of a temporary SimpleVector... bad way, to be improved. But this case happens
          // rarely ...
  {
    for (unsigned int i = 0; i < _sizeV; ++i) (*this)(i) += vIn(i);
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator+=(
    const SiconosVector& vIn) {
  // Add a part of vIn (starting from index) to the current vector.
  // vIn must be a SimpleVector.

  // At the end of the present function, index is equal to index + the dim. of the added
  // sub-vector.

  unsigned int dim = vIn.size();  // size of the block to be added.
  if (dim > _sizeV) THROW_EXCEPTION("invalid ranges");
  unsigned int currentSize;
  unsigned int index = 0;

  for (auto& it : _vect) {
    currentSize = it->size();
    *it += vIn.segment(index, currentSize);
    index += currentSize;
  }
  return *this;
}

void siconos::algebra::BlockVector::insertPtr(std::shared_ptr<SiconosVector> v) {
  if (!v) THROW_EXCEPTION("v is a nullptr vector.");

  _sizeV += v->size();
  _vect.push_back(v);
  _tabIndex->push_back(_sizeV);
}

void siconos::algebra::BlockVector::setBlock(const SiconosVector& vIn, unsigned int sizeB,
                                             unsigned int startIn, unsigned int startOut) {
  // Check dim ...
  unsigned int endOut = startOut + sizeB;

  assert(startIn < vIn.size());
  assert(startOut < size());
  assert((startIn + sizeB) <= vIn.size());
  assert(endOut <= size());

  // We look for the block of vOut that include index startOut
  unsigned int blockOutStart = 0;
  while (startOut >= (*_tabIndex)[blockOutStart] && blockOutStart < _tabIndex->size())
    blockOutStart++;
  // Relative position in the block blockOutStart.
  unsigned int posOut = startOut;
  if (blockOutStart != 0) posOut -= (*_tabIndex)[blockOutStart - 1];

  // We look for the block of vOut that include index endOut
  unsigned int blockOutEnd = blockOutStart;
  while (endOut > (*_tabIndex)[blockOutEnd] && blockOutEnd < _tabIndex->size()) blockOutEnd++;

  // => the block to be set runs from block number blockOutStart to block number blockOutEnd.

  if (blockOutEnd == blockOutStart)  //
  {
    _vect[blockOutStart]->segment(posOut, sizeB) = vIn.segment(startIn, sizeB);
  } else  // More that one block of vOut are concerned
  {
    // The current considered block ...
    auto currentBlock = _vect[blockOutStart];

    // Size of the subBlock of vOut to be set.
    size_t subSizeB = currentBlock->size() - posOut;
    unsigned int posIn = startIn;

    // Set first sub-block (currentBlock) values, between index posOut and posOut+subSizeB,
    // with vIn values from posIn to posIn+subSizeB.
    currentBlock->segment(posOut, sizeB) = vIn.segment(posIn, sizeB);

    // Other blocks, except number blockOutEnd.
    unsigned int currentBlockNum = blockOutStart + 1;
    while (currentBlockNum != blockOutEnd) {
      posIn += subSizeB;
      currentBlock = _vect[currentBlockNum];
      subSizeB = currentBlock->size();
      currentBlock->segment(0, subSizeB) = vIn.segment(posIn, subSizeB);
      currentBlockNum++;
    }
    // set last subBlock ...
    currentBlock = _vect[blockOutEnd];

    posIn += subSizeB;

    // Size of the considered sub-block
    subSizeB = endOut - (*_tabIndex)[blockOutEnd - 1];

    currentBlock->segment(0, subSizeB) = vIn.segment(posIn, subSizeB);
  }
}

double siconos::algebra::BlockVector::norm2() const {
  double d = 0;
  for (auto& it : _vect) {
    assert(it);
    d += pow(it->norm2(), 2);
  }
  return sqrt(d);
}

double siconos::algebra::BlockVector::normInf() const {
  double d = 0;
  for (auto& it : _vect) {
    assert(it);
    d = fmax(it->lpNorm<Eigen::Infinity>(), d);
  }
  return d;
}

std::shared_ptr<siconos::algebra::SiconosVector>
siconos::algebra::BlockVector::toSiconosVector() const {
  {
    if (_tabIndex->size() > 1) {
      size_t total_size = 0;
      for(auto t : *_tabIndex) {
        total_size += t;
      }
      auto vectOut = std::make_shared<SiconosVector>(total_size);
      for(auto v : _vect) {
        *vectOut << *v;
      }
      return vectOut;
    } else {
      // No copy, just a ref.
      return _vect[0];
    }
  }
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator=(
    const SiconosVector& vIn) {
  setBlock(vIn, _sizeV, 0, 0);
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator*=(double s) {
  for (auto it = begin(); it != end(); ++it) {
    (**it) *= s;
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator/=(double s) {
  for (auto it = begin(); it != end(); ++it) {
    (**it) /= s;
  }
  return *this;
}

bool siconos::algebra::isComparableTo(const BlockVector& v1, const BlockVector& v2) {
  // return:
  //  - true if both are block but with blocks which are facing each other of the same size.
  //  - false in other cases
  //
  auto& I1 = *v1.tabIndex();
  auto& I2 = *v2.tabIndex();

  return (I1 == I2);
}
