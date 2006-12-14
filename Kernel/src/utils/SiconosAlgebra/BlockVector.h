/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
 *  \version 2.0.1.
 *
 */
class BlockVector : public SiconosVector
{
private:
  // A container of pointers on SiconosVector (that are to be SimpleVector : no Block of Block allowed
  BlocksVect vect;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  std::deque<bool> isBlockAllocatedIn;

  // tabindex[i] = tabindex[i-1] + ni, ni being the size of svref[i].
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

  /** constructor with the number of Blocks and their dimension (ie all Blocks have the same dim)
  *  \param unsigned int : number of Blocks
  *  \param unsigned int : dim of each Block
  */
  BlockVector(unsigned int, unsigned int);

  /** destructor
  */
  ~BlockVector();

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
  unsigned int getNum() const;

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

  /** return false if one of the block is a null pointer
  * \return a bool
  */
  bool check() const;

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

  /** get the vector size, ie the total number of (double)
  *  elements in the vector
  *  \return unsigned int
  */
  inline unsigned int size() const
  {
    return tabIndex[tabIndex.size() - 1];
  };

  /** get the number of SimpleVector-Blocks
  *  \return unsigned int
  */
  inline unsigned int getNumberOfBlocks() const
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
  const double norm() const ;

  /** display data on standard output
  */
  void display(void) const;

  /** return i-eme SiconosVector of vect
  * \return a SimpleVector
  */
  SimpleVector getVector(unsigned int) const;

  /** return i-eme SiconosVector of vect
  * \return a pointer to a SiconosVector
  */
  SiconosVector* getVectorPtr(unsigned int);

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
  std::string toString() const;

  /** return the element vector[i]
   *  \param an unsigned int i
   *  \return a double
   */
  double getValue(unsigned int);

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
  double operator()(unsigned int) const;

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
  *  \exception SiconosVectorException
  */
  void add(const  SiconosVector&) ;

  /** add a pointer to a subvector in this vector: no reallocation nor copy.
  *  \param a pointer to SiconosVector*
  *  \exception SiconosVectorException
  */
  void addPtr(SiconosVector*) ;
};

#endif
