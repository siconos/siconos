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

/*! \file BlockVector.hpp
  \brief Object to handle vectors of vectors
*/

#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include <memory>
#include <vector>

#include "SiconosSerialization.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

/**
   "Block" vector : container (list) of SiconosVector

   A block vector is a stl vector that handles pointers to SiconosVector.

   Insertion of nullptr std::shared_ptr<SiconosVector> is not allowed.

*/
class BlockVector {
 private:
  ACCEPT_SERIALIZATION(BlockVector);

  /** Size (ie total number of scalar elements, not number of blocks) */
  Eigen::Index totalSize_ = 0;

  /** A container of pointers on SiconosVector. */
  siconos::algebra::blocks::SharedVector blocks_;

  /** blockStartIndices_[i] = blockStartIndices_[i-1] + ni, ni being the size of block[i]. */
  std::vector<Eigen::Index> blockStartIndices_ = {0};

  /* recompute the totalSize_ and blockStartIndices_ */
  void _update();

  // Rule of five
  BlockVector(const BlockVector&) = delete;
  BlockVector(BlockVector&&) noexcept = delete;
  BlockVector& operator=(BlockVector&&) noexcept = delete;

 public:
  /** default contructor
   */
  BlockVector() = default;

  /** contructor with a BlockVector of n (numberOfBlocks) blocks
   *  of the same size (dim) filled with a new vector
   *
   *  \param numberOfBlocks number of blocks
   *  \param blockSize unique dimension of all vectors in the BlockVector
   */
  BlockVector(std::size_t numberOfBlocks, Eigen::Index blockSize);

  /** contructor with a BlockVector of n (numberOfBlocks) blocks that point on nullptr
   *
   *  \param numberOfBlocks number of blocks
   */
  explicit BlockVector(std::size_t numberOfBlocks);

  /** destructor
   */
  ~BlockVector() noexcept = default;

  /** \return the size of the vector (sum of the sizes of all its blocks) */
  auto size() const { return totalSize_; };

  /** \return an iterator pointing to the first block in the container. */
  inline siconos::algebra::blocks::SharedVector::iterator begin() { return blocks_.begin(); };

  /**  \return an iterator referring to the past-the-end element in the container. */
  inline siconos::algebra::blocks::SharedVector::iterator end() { return blocks_.end(); };

  /** \return an iterator pointing to the first block in the container. */
  inline siconos::algebra::blocks::SharedVector::const_iterator begin() const {
    return blocks_.begin();
  };

  /**  \return an iterator referring to the past-the-end element in the container. */
  inline siconos::algebra::blocks::SharedVector::const_iterator end() const {
    return blocks_.end();
  };

  /** \return the complete stl container */
  inline siconos::algebra::blocks::SharedVector getAllVect() const { return blocks_; }

  /** \return the number of SiconosVectors in the container */
  inline auto numberOfBlocks() const { return blocks_.size(); };

  //   /** \return true if all SiconosVector in the container are dense **/
  //   bool isDense() const;

  /** sets all the values of the vector to 0.0 */
  void setZero();

  /** set all values of the vector component to value.
   *
   *  \param a double
   */
  void setConstant(double a);

  /** set the content of the BlockVector from an array
   *  \param size_of_data size of input array (must fit with total size of the BlockVector)
   *  \param data input array
   */
  void fill(std::size_t size_of_data, const double* data);

  /** set the content of the BlockVector from a SiconosVector (copy)
   *  \param vec input vector
   */
  void fill(const SiconosVector& vec);

  /** \return value at a given global position
   *
   *  \param globalIndex index of the required component
   */
  double operator()(Eigen::Index globalIndex) const;
  double& operator()(Eigen::Index globalIndex);

  /** get a block (SiconosVector) of the vector
   *
   *  \param pos index of the required block
   *  \return the expected block
   */
  inline auto vector(std::size_t pos) const { return blocks_[pos]; };

  /** set a block with a given vector (pointer link!)
   *
   *  \param pos index of the block to set
   *  \param v source vector to be inserted at position i
   */
  void setVectorPtr(std::size_t pos, std::shared_ptr<SiconosVector> v);

  /**
     Assignment operator
     Full copy -> each block is a new, no shared memory

      \param input the vector to be copied
      \return  BlockVector&
   */
  BlockVector& operator=(const BlockVector& input);

  /**
      Add in place operator

      \param input rhs of the operator
      \return BlockVector&
  */
  BlockVector& operator+=(const BlockVector&);

  /**
      Add a SiconosVector in place

      \param input rhs of the operator
      \return BlockVector&
  */
  BlockVector& operator+=(const SiconosVector& input);

  /**
      Subtract a SiconosVector in place

      \param input rhs of the operator
      \return BlockVector&
  */
  BlockVector& operator-=(const SiconosVector& input);

  /**
      multiply by a scalar, result in place

      \param s the scalar factor
      \return BlockVector&
  */
  BlockVector& operator*=(double s);

  /** Insert a new block (no allocation and nor copy)
   *
   *  \param v the vector to be inserted
   */
  void insertPtr(std::shared_ptr<SiconosVector> v);

  /** \return the Euclidian norm of the vector */
  double norm() const;

  /**
     Tranform a BlockVector into a SiconosVector.

     Required for plugins, that need contiguous memory for their parameters.

     \return a vector (the result depends on the number of blocks in input.
     1 block : link to first component of the container, more : copy of all components into a
     SiconosVector)
  */
  SiconosVector toSiconosVector() const;
};

// Free functions

/** \return the infinite norm of a block vector
 * \param v the input vector
 */
double normInf(const BlockVector& v);

/** display data on standard output
 *  \param v the input vector
 */
void print(const BlockVector& v);

}  // namespace siconos::algebra

#endif
