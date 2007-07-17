/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
  \brief Interface for vectors handling.

*/

#ifndef __SiconosVector__
#define __SiconosVector__

#include "SiconosAlgebra.h"
#include "SiconosVectorException.h"

/** Union of DenseVect and SparseVect pointers -
    DENSE, num = 1,
    SPARSE, num = 4
 */
union VECTOR_UBLAS_TYPE
{
  DenseVect *Dense; // num = 1
  SparseVect *Sparse; // num = 4
};

/** Abstract class to provide interface for vectors handling
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \date (creation) 07/21/2006
 *  \version 2.1.1.
 *
 *  In Siconos, a "vector" can be either a SimpleVector or a BlockVector, ie a container of several pointers to SiconosVector.
 *
 * Remark: a lot of functions in the interface are "reserved to ..." Block or Simple.
 * See the documentation in these classes for more details.
 *
 * You can find an overview on how to build and use vectors and matrices in \ref GS_SicAlgebra .
 *
 */
class SiconosVector
{

protected:

  /** Size (ie total number of scalar elements, not number of blocks) */
  unsigned int sizeV;

  /** A number to specify the type of the vector: (block or ublas-type)
   * 0 -> BlockVector, 1 -> DenseVect, 4 -> SparseVect
   * Note: 4 for sparse to keep the same num as for matrices.
   */
  unsigned int num;

  /** default constructor */
  SiconosVector() {};

  /** constructor with type-number
      \param unsigned int, type-number of the vector
  */
  SiconosVector(unsigned int);

  /** basic constructor
      \param unsigned int, type-number of the vector
      \param size of the vector */
  SiconosVector(unsigned int, unsigned int);

public:

  /** Destructor.
   */
  virtual ~SiconosVector() {};

  /** true if the vector is block else false.
   * \return a bool.
   */
  inline const bool isBlock() const
  {
    if (num == 0) return true ;
    else return false;
  }

  /** get the vector size, ie the total number of (double) elements in the vector
   *  \return unsigned int
   */
  inline const unsigned int size() const
  {
    return sizeV;
  };

  /** Get the type number of the current vector.
   * \return an unsigned int
   */
  inline const unsigned int getNum() const
  {
    return num;
  };

  /** get the number of SimpleVector-Blocks - only usefull for BlockVector.
   *  \return unsigned int
   */
  inline virtual const unsigned int getNumberOfBlocks() const
  {
    return 1;
  };

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

  /** get the ublas embedded vector if it's type is Dense
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect
   */
  virtual const DenseVect getDense(unsigned int = 0) const = 0;

  /** get the ublas embedded vector if it's type is Sparse
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect
   */
  virtual const SparseVect getSparse(unsigned int = 0) const = 0;

  /** get a pointer to the ublas embedded vector if it's type is Dense
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect*
   */
  virtual DenseVect* getDensePtr(unsigned int = 0) const = 0;

  /** get a pointer to the ublas embedded vector if it's type is Sparse
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect*
   */
  virtual SparseVect* getSparsePtr(unsigned int = 0) const = 0;

  /** return the array of double values of the vector
   *  \param unsigned int: vector position (only for block vector)
   *  \return double* : the pointer on the array
   */
  virtual double* getArray(unsigned int = 0) const = 0;

  /** get block starting at "pos" (first argument) and write it in v (second arg)
   *  \param pos an int, position of the first element of the required block
   *  \param v a SiconosVector*, in-out parameter.
   */
  //  virtual void getBlock(unsigned int, SiconosVector*) const = 0;

  /** sets all the values of the vector to 0.0
   */
  virtual void zero() = 0;

  /** resize the vector with nbcol columns. The existing elements of the vector are preseved when specified.
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

  /** if this is a block vector return SiconosVector* number i (arg), else return this.
   * \param i, unsigned int
   * \return a pointer to a SiconosVector
   */
  virtual SiconosVector* getVectorPtr(unsigned int) = 0;

  /** if this is a block vector return SiconosVector* number i (arg), else return this.
   * \param i, unsigned int
   * \return a pointer to a SiconosVector
   */
  virtual const SiconosVector* getVectorPtr(unsigned int) const = 0;

  /** set SiconosVector number i (copy) with v (second arg) - Useful only for BlockVector (else equivalent to a single copy)
   * \param i, unsigned int, block number (0 for SimpleVector)
   * \param v, a SiconosVector
   */
  virtual void setVector(unsigned int, const SiconosVector&) = 0;

  /** set SiconosVector number i (pointer link) with v (second arg) - Useful only for BlockVector
   * \param i, unsigned int: block number (0 for SimpleVector)
   * \param v, a pointer to a SiconosVector
   */
  virtual void setVectorPtr(unsigned int, SiconosVector*) = 0;

  /** set all values of the vector component to input value.
   * \param a double
   */
  virtual void fill(double) = 0;

  /** reserved to BlockVector - get the index tab
   * \return a standard vector of int
   */
  virtual Index getTabIndex() const ;

  /** reserved to BlockVector - get the index tab
   * \return a pointer to a standard vector of int
   */
  virtual const Index * getTabIndexPtr() const ;

  /** put data of the vector into a string
   */
  virtual const std::string toString() const = 0;

  /* Note: in the following functions, index is a general one;
   that means that for a SimpleVector v, v(i) is index i element but
   for a BlockVector w that contains 2 SiconosVector of size 3
   w(4) corresponds to the first element of the second vector. */

  /** get the element at position i in the vector
   *  \param an unsigned int i
   *  \return a double
   */
  virtual const double getValue(unsigned int) const = 0;

  /** set the element at position i in the vector.
   *  \param an unsigned int i
   *  \param the value
   */
  virtual void setValue(unsigned int, double) =  0;

  /** get the element at position i in the vector
   *  \param an integer i
   *  \return a double
   */
  virtual double& operator()(unsigned int) = 0;

  /** get the element at position i in the vector
   *  \param an integer i
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

  /** assignment
   *  \param SiconosVector : the vector to be copied
   */
  virtual SiconosVector& operator = (const SiconosVector&) = 0;

  /** assignment
   *  \param a DenseVect : the vector to be copied
   */
  virtual SiconosVector& operator = (const DenseVect&) = 0;

  /** assignment
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

  /** multiply the current vector with a scalar
   *  \param template, double, int ...
   */
  template <class T> SiconosVector& operator *= (const T& s)
  {
    if (num == 0)
    {
      BlockVectIterator it;
      for (it = begin(); it != end(); ++it)
        (**it) *= s;
    }
    else if (num == 1)
      //atlas::scal((double)m,*vect.Dense);
      *getDensePtr() *= s;
    else
      *getSparsePtr() *= s;
    return *this;
  }

  /** divide the current vector with a scalar
   *  \param template, double, int ...
   */
  template <class T> SiconosVector& operator /= (const T& s)
  {
    if (num == 0)
    {
      BlockVectIterator it;
      for (it = begin(); it != end(); ++it)
        (**it) /= s;
    }
    else if (num == 1)
      //atlas::scal((double)m,*vect.Dense);
      *getDensePtr() /= s;
    else
      *getSparsePtr() /= s;
    return *this;
  }

  /** reserved to BlockVector - Insert a subvector in this vector: allocation and copy
   *  \param SiconosVector& v : the vector to be inserted
   */
  virtual inline void insert(const  SiconosVector&)
  {
    SiconosVectorException::selfThrow("SiconosVector::insert() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
  };

  /** reserved to BlockVector - Insert a pointer to a subvector in this vector: no reallocation nor copy.
   *  \param a pointer to SiconosVector*
   */
  virtual inline void insertPtr(SiconosVector*)
  {
    SiconosVectorException::selfThrow("SiconosVector::insertPtr() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
  };

  /** Compare two (block) vectors: true if they have the same number of blocks and if
      blocks which are facing each other have the same size;
      always true if one of the two is a SimpleVector.
      \param a SiconosVector*.
      \param a SiconosVector*.
  */
  friend const bool isComparableTo(const SiconosVector *, const SiconosVector *);

  /** Swap x and y contents, using atlas swap.*/
  friend void swap(SiconosVector&, SiconosVector&);

};

#endif
