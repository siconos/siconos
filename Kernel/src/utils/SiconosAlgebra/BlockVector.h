/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

/*! \file BlockVector.h

*/

#ifndef MYBLOCKVECTOR_H
#define MYBLOCKVECTOR_H

#include "SiconosVector.h"

class SimpleVector;


/** Object to handle block-vectors (ie list of SiconosVector*)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *
 * A block vector is a stl vector that handles pointers to SiconosVector.
 *
 * Insertion of NULL SiconosVector* is not allowed.
 *
 */
class BlockVector : public SiconosVector
{
private:
  // A container of pointers on SiconosVector.
  BlocksVect vect;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  std::deque<bool> isBlockAllocatedIn;

  /** tabindex[i] = tabindex[i-1] + ni, ni being the size of svref[i]. */
  Index tabIndex;

public:

  /** default contructor
   */
  BlockVector();

  /** contructor from data by function read call
   *  \param a string
   *  \param a bool
   */
  BlockVector(const std::string&, const bool);

  /** contructor with a SiconosVector (copy)
   *  \param SiconosVector& v
   */
  BlockVector(const SiconosVector&);

  /** copy contructor
   *  \param BlockVector& v
   */
  BlockVector(const BlockVector&);

  /** contructor with a 2 SiconosVectors
   *  \param SiconosVector* v1
   *  \param SiconosVector* v2
   */
  BlockVector(SiconosVector*, SiconosVector*);

  /** constructor with the number of Blocks and their dimension (ie all Blocks have the same dim AND are SimpleVectors)
   *  \param unsigned int : number of Blocks
   *  \param unsigned int : dim of each Block
   */
  BlockVector(unsigned int, unsigned int);

  /** constructor for block of blocks (all subblocks have the same number of blocks)
   *  \param unsigned int n1 : number of main blocks
   *  \param unsigned int n2 : number of blocks in each main block
   *  \param vector<SiconosVector*> the blocks (n2 first vectors are in first block, n2 next in second block and so on)
   */
  BlockVector(unsigned int, unsigned int, std::vector<SiconosVector*>);

  /** destructor
   */
  ~BlockVector();

  /** iterator equal to vect.begin
      \return a BlocksVectIterator
  */
  inline BlockVectIterator begin()
  {
    return vect.begin();
  };

  /** iterator equal to vect.end
      \return a BlocksVectIterator
  */
  inline BlockVectIterator end()
  {
    return vect.end();
  };

  /** const iterator equal to vect.begin
      \param a BlocksVectIterator
  */
  inline ConstBlockVectIterator begin() const
  {
    return vect.begin();
  };

  /** const iterator equal to vect.end
      \param a BlocksVectIterator
  */
  inline ConstBlockVectIterator end() const
  {
    return vect.end();
  } ;

  /** get vect, ie all the vectors of the object
   * \return a BlocksVect
   */
  inline BlocksVect getAllVect() const
  {
    return vect;
  }

  /** get the attribute num of current vector
   * \return an unsigned int.
   */
  const unsigned int getNum() const;

  /** get the attribute if it's type is DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect
   */
  const DenseVect getDense(unsigned int = 0) const;

  /** get the attribute if it's type is SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect
   */
  const SparseVect getSparse(unsigned int = 0) const;

  /** get a pointer on DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect*
   */
  DenseVect* getDensePtr(unsigned int = 0) const;

  /** get a pointer on SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect*
   */
  SparseVect* getSparsePtr(unsigned int = 0) const;

  /** return the array of double values of the vector
   *  \exception SiconosVectorException
   *  \param unsigned int: vector position (only for block vector)
   *  \return double* : the pointer on the array
   */
  double* getArray(unsigned int = 0) const;

  /** get block starting at pos and of size block.size()
   *  \param an int, position of the first element of the required block
   *  \param a SiconosVector, in-out parameter.
   */
  void getBlock(unsigned int, SiconosVector&) const;

  /** sets all the values of the vector to 0.0
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   */
  void zero();

  /** get the number of SimpleVector-Blocks
   *  \return unsigned int
   */
  inline const unsigned int getNumberOfBlocks() const
  {
    return tabIndex.size();
  };

  /** resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   *  \exception SiconosVectorException
   */
  void resize(unsigned int, bool = true);

  /** compute the infinite norm of the vector
   *  \return a double
   */
  const double normInf() const;

  /** return the Euclidian norm of the vector
   *  \return a double
   */
  const double norm2() const ;

  /** display data on standard output
   */
  void display(void) const;

  /** return i-eme SiconosVector of vect
   * \param unsigned int: block number
   * \return a SimpleVector
   */
  SimpleVector getVector(unsigned int) const;

  /** return i-eme SiconosVector of vect
   * \param unsigned int: block number
   * \return a pointer to a SiconosVector
   */
  SiconosVector* getVectorPtr(unsigned int);

  /** set i-eme SiconosVector of vect (copy)
   * \param unsigned int: block number
   * \param a SiconosVector
   */
  void setVector(unsigned int, const SiconosVector&);

  /** set i-eme SiconosVector of vect (pointer link)
   * \param unsigned int: block number
   * \param a pointer to a SiconosVector
   */
  void setVectorPtr(unsigned int, SiconosVector*);

  /** set all values of the vector component to value.
   * \param a double
   */
  void fill(double);

  /** get the index tab
   * \return a standard vector of int
   */
  inline Index getTabIndex() const
  {
    return tabIndex;
  }

  /** put data of the vector into a string
   */
  const std::string toString() const;

  /** return the element vector[i]
   *  \param an unsigned int i
   *  \return a double
   */
  const double getValue(unsigned int);

  /** set the element vector[i]
   *  \param an unsigned int i
   *  \param the value
   */
  void setValue(unsigned int, double);

  /** get the element at position i (warning: absolute position.)
   *  \param an unsigned integer i
   *  \return a reference to a double
   */
  double& operator()(unsigned int) ;

  /** get the element at position i (warning: absolute position.)
   *  \param an unsigned integer i
   *  \return a double
   */
  const double operator()(unsigned int) const;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SiconosVector*
   */
  SiconosVector* operator [](unsigned int) ;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SiconosVector*
   */
  const SiconosVector* operator [](unsigned int) const;

  /** operator =
   *  \param SiconosVector : the vector to be copied
   */
  BlockVector& operator =(const SiconosVector&);

  /** operator =
   *  \param SiconosVector : the vector to be copied
   */
  BlockVector& operator =(const BlockVector&);

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  BlockVector& operator = (const DenseVect&);

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  BlockVector& operator = (const SparseVect&);

  /** operator +=
   *  \param SiconosVector : a vector to add
   */
  BlockVector& operator +=(const SiconosVector&);

  /** operator -=
   *  \param SiconosVector : a vector to subtract
   */
  BlockVector& operator -=(const SiconosVector&);

  /** operator /=
   *  \param double, a scalar
   */
  BlockVector& operator /=(double);

  /** operator /=
   *  \param int, a scalar
   */
  BlockVector& operator /=(int);

  /** multiply the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   */
  BlockVector& operator *=(double);

  /** multiply the current vector with an int
   *  \param an int
   *  \exception SiconosVectorException
   */
  BlockVector& operator *=(int);

  /** add a subvector in this vector: allocation and copy
   *  \param SiconosVector& v : the vector to add
   */
  void add(const  SiconosVector&) ;

  /** add a pointer to a subvector in this vector: no reallocation nor copy.
   *  \param a pointer to SiconosVector*
   */
  void addPtr(SiconosVector*) ;
};

#endif
