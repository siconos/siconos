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

#include "BlockVector.hpp"

#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>

#include <cstddef>
#include <iostream>
#include <vector>

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "siconos_debug.h"
// //#define DEBUG_STDOUT
// //#define DEBUG_MESSAGES
//

siconos::algebra::BlockVector::BlockVector(std::size_t numberOfBlocks,
                                           Eigen::Index blockSize) {
  blockStartIndices_.reserve(numberOfBlocks + 1);
  blocks_.reserve(numberOfBlocks);
  totalSize_ = 0;
  for (std::size_t i = 0; i < numberOfBlocks; ++i) {
    blocks_.emplace_back(std::make_shared<SiconosVector>(blockSize));
    totalSize_ += blockSize;
    blockStartIndices_.push_back(totalSize_);
  }
}

siconos::algebra::BlockVector::BlockVector(std::size_t numberOfBlocks) {
  blockStartIndices_.resize(numberOfBlocks + 1);
  blocks_.resize(numberOfBlocks);
  blockStartIndices_[0] = 0;
}

void siconos::algebra::BlockVector::_update() {
  totalSize_ = 0;
  blockStartIndices_.clear();
  blockStartIndices_.reserve(blocks_.size() + 1);
  blockStartIndices_.push_back(0);
  for (const auto& it : blocks_) {
    if (it) {
      totalSize_ += it->size();
    }
    blockStartIndices_.push_back(totalSize_);
  }
}

// ===========================
//       fill vector
// ===========================

void siconos::algebra::BlockVector::setZero() {
  for (auto& it : blocks_) it->setZero();
}

void siconos::algebra::BlockVector::setConstant(double value) {
  for (auto& it : blocks_) {
    if (it) {
      it->setConstant(value);
    }
  }
}

void siconos::algebra::BlockVector::fill(std::size_t size_of_data, const double* data) {
  assert(data && "BlockVector::fill: input data pointer is null");
  assert(static_cast<Eigen::Index>(size_of_data) == totalSize_ &&
         "total size of data does not match expected BlockVector size");
  // Warning: we do not check data and assumes it is properly allocated
  std::size_t offset = 0;
  for (std::size_t nb = 0; nb < blocks_.size(); ++nb) {
    auto& block = *(blocks_[nb]);
    Eigen::Index size = block.size();
    ConstMapVectorType src{data + offset, size};
    block = src;  // copy data viewed by the map into the current block
    offset += size;
  }
}

void siconos::algebra::BlockVector::fill(const SiconosVector& input) {
  assert(input.size() == totalSize_ &&
         "total size of input vector does not match expected BlockVector size");
  // Warning: we do not check data and assumes it is properly allocated
  std::size_t offset = 0;
  for (std::size_t nb = 0; nb < blocks_.size(); ++nb) {
    auto& block = *(blocks_[nb]);
    auto blockSize = block.size();
    block = input.segment(offset, blockSize);  // copy
    offset += blockSize;
  }
}

double siconos::algebra::BlockVector::operator()(Eigen::Index globalIndex) const {
  // Find the block ...
  for (std::size_t i = 0; i < blockStartIndices_.size() - 1; ++i) {
    if (globalIndex >= blockStartIndices_[i] && globalIndex < blockStartIndices_[i + 1]) {
      // Then the index in the block
      std::size_t localIndex = globalIndex - blockStartIndices_[i];
      return (*blocks_[i])(localIndex);
    }
  }
  throw std::out_of_range("Index out of range in BlockVector");
}

void siconos::algebra::BlockVector::setVectorPtr(std::size_t pos,
                                                 std::shared_ptr<SiconosVector> v) {
  assert(pos < blocks_.size() && "insertion out of vector size");
  blocks_[pos] = v;
  _update();
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator=(
    const BlockVector& other) {
  if (this == &other) return *this;
  assert(totalSize_ == other.size());
  assert(blocks_.size() == other.numberOfBlocks());
  assert(blockStartIndices_ == other.blockStartIndices_);
  // totalSize_ = other.totalSize_;
  //  blockStartIndices_ = other.blockStartIndices_;

  // Deep copy if blocks
  blocks_.clear();
  blocks_.reserve(other.blocks_.size());
  for (const auto& block : other.blocks_) {
    if (block)
      blocks_.emplace_back(std::make_shared<SiconosVector>(*block));
    else
      blocks_.emplace_back(nullptr);
  }

  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator+=(
    const BlockVector& input) {
  assert(totalSize_ == input.size());
  assert(blocks_.size() == input.numberOfBlocks());

  std::size_t pos = 0;
  for (auto& block : blocks_) {
    *block += *(input.vector(pos));
    pos++;
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator+=(
    const SiconosVector& input) {
  assert(totalSize_ == input.size());

  Eigen::Index currentSize;
  Eigen::Index offset = 0;

  for (auto& it : blocks_) {
    currentSize = it->size();
    *it += input.segment(offset, currentSize);
    offset += currentSize;
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator-=(
    const SiconosVector& input) {
  assert(totalSize_ == input.size());

  Eigen::Index currentSize;
  Eigen::Index offset = 0;

  for (auto& it : blocks_) {
    currentSize = it->size();
    *it -= input.segment(offset, currentSize);
    offset += currentSize;
  }
  return *this;
}

siconos::algebra::BlockVector& siconos::algebra::BlockVector::operator*=(double s) {
  for (auto& block : blocks_) *block *= s;
  return *this;
}

void siconos::algebra::BlockVector::insertPtr(std::shared_ptr<SiconosVector> v) {
  if (!v) THROW_EXCEPTION("v is a nullptr vector.");

  totalSize_ += v->size();
  blocks_.push_back(v);
  blockStartIndices_.push_back(totalSize_);
}

double siconos::algebra::BlockVector::norm() const {
  double d = 0;
  for (auto& v : blocks_) {

    d += v->squaredNorm();
  }
  return std::sqrt(d);
}

siconos::algebra::SiconosVector siconos::algebra::BlockVector::toSiconosVector() const {
  SiconosVector result(totalSize_);
  size_t currentIndex = 0;
  for (const auto& block : blocks_) {
    if (block) {
      auto vecSize = block->size();
      result.segment(currentIndex, vecSize) = *block;
      currentIndex += vecSize;
    }
  }
  return result;
}

// Free functions

double siconos::algebra::normInf(const siconos::algebra::BlockVector& vect) {
  double d = 0;
  for (auto& v : vect) {
    d = std::max(v->lpNorm<Eigen::Infinity>(), d);
  }
  return d;
}

//=====================
// screen display
//=====================

void siconos::algebra::print(const siconos::algebra::BlockVector& vect) {
  std::cout << "=======> Block Vector Display (" << vect.numberOfBlocks()
            << " block(s)): " << std::endl;
  for (auto& it : vect) {
    DEBUG_EXPR(std::cout << "(*it)" << (*it) << std::endl;);
    if (it) {
      siconos::algebra::print(*it);
      std::cout << std::endl;
    } else
      std::cout << "(*it)-> nullptr" << std::endl;
  }
}