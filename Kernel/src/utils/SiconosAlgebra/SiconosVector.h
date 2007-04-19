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


/*! \file SiconosVector.h
  \brief SiconosVector class

*/

#ifndef __SiconosVector__
#define __SiconosVector__

#include "SiconosAlgebra.h"

/** Union of DenseVect pointer and SparseVect pointer
 */
union Vect
{
  DenseVect *Dense; // num = 1
  SparseVect *Sparse; // num = 4
};

/** Abstract class to provide interface for vectors handling
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \date (creation) 07/21/2006
 *  \version 2.0.1.
 *  vector can be either a SimpleVector or a BlockVector, ie a container of several SimpleVector
 *  See documentation of these derivated classes for more details
 *
 */
class SiconosVector
{
protected:

  /**\var isBlockVector
   * bool to check the type of the current vector; true if block else false.
   */
  bool isBlockVector;

  /** Size (ie total number of scalar elements, not number of blocks) */
  unsigned int sizeV;

  /** default constructor */
  SiconosVector(bool = false);

public:

  /** true if the vector is block else false.
   * \return a bool.
   */
  inline const bool isBlock(void) const
  {
    return isBlockVector;
  }

  /** Destructor.
   */
  virtual ~SiconosVector();

  /** get the attribute num of current vector
   * \return an unsigned int.
   */
  virtual const unsigned int getNum() const = 0;

  /** reserved to BlockVector
      \return a BlocksVectIterator
  */
  virtual BlockVectIterator begin();

  /** reserved to BlockVector
      \return a BlocksVectIterator
  */
  virtual BlockVectIterator end();

  /** reserved to BlockVector
      \return a ConstBlocksVectIterator
  */
  virtual ConstBlockVectIterator begin() const;

  /** reserved to BlockVector
      \return a ConstBlocksVectIterator
  */
  virtual ConstBlockVectIterator end() const;

  /** get the attribute if it's type is DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect
   */
  virtual const DenseVect getDense(unsigned int = 0) const = 0;

  /** get the attribute if it's type is SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect
   */
  virtual const SparseVect getSparse(unsigned int = 0) const = 0;

  /** get a pointer on DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect*
   */
  virtual DenseVect* getDensePtr(unsigned int = 0) const = 0;

  /** get a pointer on SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect*
   */
  virtual SparseVect* getSparsePtr(unsigned int = 0) const = 0;

  /** return the array of double values of the vector
   *  \exception SiconosVectorException
   *  \param unsigned int: vector position (only for block vector)
   *  \return double* : the pointer on the array
   */
  virtual double* getArray(unsigned int = 0) const = 0;

  /** get block starting at pos and of size block.size()
   *  \param an int, position of the first element of the required block
   *  \param a SiconosVector, in-out parameter.
   */
  virtual void getBlock(unsigned int, SiconosVector&) const = 0;

  /** sets all the values of the vector to 0.0
   */
  virtual void zero() = 0;

  /** get the vector size, ie the total number of (double) elements in the vector
   *  \return unsigned int
   */
  inline const unsigned int size() const
  {
    return sizeV;
  };

  /** resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   *  \param: dim of the resized vector
   *  \param: a bool, true if old values are preserved else false. Default = true.
   *  \exception for Block Vector, resizing not allowed.
   */
  virtual void resize(unsigned int, bool = true) = 0;

  /** compute the infinite norm of the vector
   *  \return a double
   */
  virtual const double normInf() const = 0;

  /** return the Euclidian norm of the vector
   *  \return a double
   */
  virtual const double norm2() const = 0 ;

  /** display data on standard output
   */
  virtual void display(void)const = 0;

  /** if this is a block vector return i-eme SimpleVector, else return this.
   * \return a pointer to a SimpleVector
   */
  virtual SiconosVector* getVectorPtr(unsigned int) = 0;

  /** set all values of the vector component to value.
   * \param a double
   */
  virtual void fill(double) = 0;

  /** get the index tab
   * \return a standard vector of int
   */
  virtual Index getTabIndex() const ;

  /** get the number of SimpleVector-Blocks - only usefull for BlockVector.
   *  \return unsigned int
   */
  virtual const unsigned int getNumberOfBlocks() const
  {
    return 1;
  };

  /** put data of the vector into a string
   */
  virtual const std::string toString() const = 0;

  /* Note: in the following functions, index is a general one;
   that means that for a SimpleVector v, v(i) is index i element but
   for a BlockVector w that contains 2 SiconosVector of size 3
   w(4) corresponds to the first element of the second vector. */

  /** return the element vector[i]
   *  \param an unsigned int i
   *  \return a double
   */
  virtual const double getValue(unsigned int) = 0;

  /** set the element vector[i]
   *  \param an unsigned int i
   *  \param the value
   */
  virtual void setValue(unsigned int, double) =  0;

  /** get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  virtual double& operator()(unsigned int) = 0;

  /** get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  virtual const double operator()(unsigned int) const = 0;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SiconosVector*
   */
  virtual SiconosVector* operator [](unsigned int) = 0;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SiconosVector*
   */
  virtual const SiconosVector* operator [](unsigned int) const = 0;

  /** operator =
   *  \param SiconosVector : the vector to be copied
   */
  virtual SiconosVector& operator =(const SiconosVector&) = 0;

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  virtual SiconosVector& operator = (const DenseVect&) = 0;

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  virtual SiconosVector& operator = (const SparseVect&) = 0;

  /** operator +=
   *  \param SiconosVector : a vector to add
   */
  virtual SiconosVector& operator +=(const SiconosVector&) = 0;


  /** operator -=
   *  \param SiconosVector : a vector to subtract
   */
  virtual SiconosVector& operator -=(const SiconosVector&) = 0;

  /** operator /=

  *  \param double, a scalar
  */
  virtual SiconosVector& operator /=(double) = 0;


  /** operator /=
   *  \param int, a scalar
   */
  virtual SiconosVector& operator /=(int) = 0;

  /** multiply the current vector with a double
   *  \param a double
   */
  virtual SiconosVector& operator *=(double) = 0;

  /** multiply the current vector with an int
   *  \param an int
   */
  virtual SiconosVector& operator *=(int) = 0;

  /** add a subvector in this vector: allocation and copy
   *  \param SiconosVector& v : the vector to add
   */
  virtual void add(const  SiconosVector&) ;

  /** add a pointer to a subvector in this vector: no reallocation nor copy.
   *  \param a pointer to SiconosVector*
   */
  virtual void addPtr(SiconosVector*) ;

  /** computes this = x + y with atlas xpy .
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  virtual void xpy(const SiconosVector &, const SiconosVector &);

  /** computes this = ax + by with atlas axpby .
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  virtual void axpby(double, const SiconosVector&, double, const SiconosVector&);

  /** computes this = a*x with atlas scal.
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  virtual void scal(double, const SiconosVector&);

  /** computes this = a*x + y with atlas axpy.
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  virtual void axpy(double, const SiconosVector&, const SiconosVector&);

};

#endif
