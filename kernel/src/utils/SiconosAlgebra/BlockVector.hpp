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

/*! \file BlockVector.hpp
  \brief Object to handle vectors of vectors
*/

#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include "SiconosAlgebraTypeDef.hpp"

class SiconosVector;

/** "Block" vector : container (list) of SiconosVector
 *
 * A block vector is a stl vector that handles pointers to SiconosVector.
 *
 * Insertion of nullptr SP::SiconosVector is not allowed.
 *
 */
class BlockVector
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(BlockVector);

  /** Size (ie total number of scalar elements, not number of blocks) */
  unsigned int _sizeV = 0;

  /** A container of pointers on SiconosVector. */
  VectorOfVectors _vect;

  /** tabindex[i] = tabindex[i-1] + ni, ni being the size of block[i]. */
  SP::Index _tabIndex;

  /* recompute the _sizeV and _tabIndex */
  void _update();

public:

  /** default contructor
   */
  BlockVector();

  /** copy contructor
   *  \param v BlockVector&
   */
  BlockVector(const BlockVector& v);

  /** contructor with a 2 SiconosVectors
   *  \param v1 first vector
   *  \param v2 second vector
   */
  BlockVector(SP::SiconosVector v1, SP::SiconosVector v2);

  /** contructor with a BlockVector of n (numberOfBlocks) blocks
   * of the same size (dim) filled with a new vector
   *  \param numberOfBlocks number of blocks
   *  \param dim dimension of the vector
   */
  BlockVector(unsigned int numberOfBlocks, unsigned int dim);

  /** contructor with a BlockVector of n (numberOfBlocks) blocks that point on nullptr
   *  \param numberOfBlocks number of blocks
   */
  BlockVector(unsigned int numberOfBlocks);

  /** destructor
   */
  ~BlockVector(){};

  /** Set a subblock of the current vector with the content (copy) of a SiconosVector
      \param input the vector to be copied
      \param size_block size of the block to be copied
      \param start_in starting position in input vector of the block to be copied
      \param start_out starting position in current vector of the block to be filled in.
   */
  void setBlock(const SiconosVector& input, unsigned int size_block, unsigned int start_in, unsigned int start_out);

  /** \return the size of the vector (sum of the sizes of all its blocks) */
  unsigned int size() const
  {
    return _sizeV;
  };


  /** \return an iterator pointing to the first block in the container. */
  inline VectorOfVectors::iterator begin()
  {
    return _vect.begin();
  };

  /**  \return an iterator referring to the past-the-end element in the container. */
  inline VectorOfVectors::iterator end()
  {
    return _vect.end();
  };

  /** \return an iterator pointing to the first block in the container. */
  inline VectorOfVectors::const_iterator begin() const
  {
    return _vect.begin();
  };

  /**  \return an iterator referring to the past-the-end element in the container. */
  inline VectorOfVectors::const_iterator end() const
  {
    return _vect.end();
  } ;

  /** \return the complete stl container */
  inline VectorOfVectors getAllVect() const
  {
    return _vect;
  }

  /** \return the number of SiconosVectors in the container */
  inline unsigned int numberOfBlocks() const
  {
    return _tabIndex->size();
  };

  /** \return true if all SiconosVector in the container are dense **/
  bool isDense() const;

  /** sets all the values of the vector to 0.0 */
  void zero();

  /** set all values of the vector component to value.
   * \param a double
   */
  void fill(double a);

  /** display data on standard output */
  void display(void) const;

  /** put data of the vector into a std::string
   * \return std::string
   */
  std::string toString() const;

  /** Get a component of the vector
   *  \param i index of the required component
   *  \return the component value
   */
  double getValue(unsigned int i) const;

  /** set a component of the vector
   *  \param i index of the required component
   *  \param value of the component
   */
  void setValue(unsigned int i, double value);

  /** get a component of the vector
   *  \param i index of the required component
   *  \return value of the component
   */
  double& operator()(unsigned int i) ;

  /** get a component of the vector
   *  \param i index of the required component
   *  \return value of the component
   */
  double operator()(unsigned int i) const;

  /** get a block (SiconosVector) of the vector
   * \param pos index of the required block
   * \return the expected block
   */
  inline SP::SiconosVector vector(unsigned int pos)
  {
    return _vect[pos];
  };

  /** gets a block (SiconosVector) of the vector
   * \param pos index of the required block
   * \return the expected block
   */
  inline SPC::SiconosVector vector(unsigned int pos) const
  {
    return _vect[pos];
  };

  /** get a block (SiconosVector) of the vector
   * \param pos index of the required block
   * \return the expected block
   */
  SP::SiconosVector operator [](unsigned int pos) ;

  /** get a block (SiconosVector) of the vector
   * \param pos index of the required block
   * \return the expected block
   */
  SPC::SiconosVector operator [](unsigned int pos) const;

  /** set a block with a given vector (copy!)
   * \param pos index of the block to set
   * \param v source vector to be copied at position i
   */
  void setVector(unsigned int pos, const SiconosVector& v);

  /** set a block with a given vector (pointer link!)
   * \param pos index of the block to set
   * \param v source vector to be inserted at position i
   */
  void setVectorPtr(unsigned int pos, SP::SiconosVector v);

  /** Fill the container with a list of SiconosVector.
      Warning: pointer links, no copy
      \param v the vectors to be inserted
   */
  void setAllVect(VectorOfVectors& v);

  /** \return a pointer to the index tab
   */
  inline const SP::Index tabIndex() const
  {
    return _tabIndex;
  }

  /** get the number of the vector that handles element at position "pos"
   *  \param pos unsigned int, position of the element
   *  \return unsigned int number of the searched vector
   */
  unsigned int getNumVectorAtPos(unsigned int pos) const;

  /** Assignment operator
      \param vIn the vector to be copied
      \return  BlockVector&
  */
  BlockVector& operator =(const BlockVector& vIn);

  /** Assignment operator
      \param data data to put in the BlockVector
      \return  BlockVector&
  */
  BlockVector& operator = (const double* data);

  /** Assignment operator
   *  \param vIn the vector to be copied
   * \return  BlockVector&
   */
  BlockVector& operator =(const SiconosVector& vIn);

  /** Subtract in place operator
      \param vIn rhs of the operator
      \return BlockVector&
  */
  BlockVector& operator -=(const BlockVector& vIn);

  /** Add in place operator
      \param vIn rhs of the operator
      \return BlockVector&
  */
  BlockVector& operator +=(const BlockVector&);

  /** Add in place operator
      \param vIn rhs of the operator
      \return BlockVector&
  */
  BlockVector& operator += (const SiconosVector& vIn);

  /** Subtract in place operator
      \param vIn rhs of the operator
      \return BlockVector&
  */
  BlockVector& operator -= (const SiconosVector& vIn);

  /** multiply by a scalar, result in place
      \param s the scalar factor
      \return BlockVector&
  */
  BlockVector& operator *= (double s);

  /** divide by a scalar, result in place
      \param s the scalar factor
      \return BlockVector&
  */
  BlockVector& operator /= (double s);

  // /** Insert a new block (allocation and copy)
  // *  \param v the vector to be inserted
  // */
  // void insert(const SiconosVector& v) ;

  /** Insert a new block (no allocation and nor copy)
   *  \param v the vector to be inserted
   */
  void insertPtr(SP::SiconosVector v);

  /** \return the Euclidian norm of the vector */
  double norm2() const;

  /** \return the infinite norm of the vector */
  double normInf() const;

  /** Tranform a BlockVector into a SiconosVector.
      Required for plugins, that need contiguous memory for their parameters.
      \return a vector (the result depends on the number of blocks in input.
      1 block : link to first component of the container, more : copy of all components into a SiconosVector)
  */
  SP::SiconosVector prepareVectorForPlugin() const;

  /** \defgroup BlockVectorFriends

      List of friend functions of the BlockVector class

      @{
  */

  /** offstream operator
   * \param os An output stream
   * \param bv a BlockVector
   * \return The same output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const BlockVector& bv);

   /** End of Friend functions group @} */

  ACCEPT_NONVIRTUAL_VISITORS();

};

#endif
