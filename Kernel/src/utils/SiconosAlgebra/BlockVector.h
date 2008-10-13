/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
  \brief Object to handle vectors of vectors ( ... of vectors)
*/

#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include "SiconosVector.h"

class SimpleVector;

/** Object to handle block-vectors (ie list of SP::SiconosVector)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *
 * A block vector is a stl vector that handles pointers to SiconosVector.
 *
 * Insertion of NULL SP::SiconosVector is not allowed.
 *
 */
class BlockVector : public SiconosVector , public boost::enable_shared_from_this<BlockVector>
{
private:

  /** A container of pointers on SiconosVector. */
  VectorOfVectors vect;

  /** tabindex[i] = tabindex[i-1] + ni, ni being the size of svref[i]. */
  Index * tabIndex;

public:

  /** default contructor
   */
  BlockVector();

  /** contructor from data by function read call
   *  \param a string
   *  \param a bool
   */
  BlockVector(const std::string&, const bool);

  /** copy contructor
   *  \param BlockVector& v
   */
  BlockVector(const BlockVector&);

  /** contructor with a SiconosVector (copy)
   *  \param SiconosVector& v
   */
  BlockVector(const SiconosVector&);

  /** contructor with a 2 SiconosVectors
   *  \param SP::SiconosVector v1
   *  \param SP::SiconosVector v2
   */
  BlockVector(SP::SiconosVector, SP::SiconosVector);

  /** constructor with the number of Blocks and their dimension (ie all Blocks have the same dim AND are SimpleVectors)
   *  \param unsigned int : number of Blocks
   *  \param unsigned int : dim of each Block
   */
  BlockVector(unsigned int, unsigned int);

  /** destructor
   */
  ~BlockVector();

  /** iterator equal to vect.begin
      \return a VectorOfVectors::iterator
  */
  inline VectorOfVectors::iterator begin()
  {
    return vect.begin();
  };

  /** iterator equal to vect.end
      \return a VectorOfVectors::iterator
  */
  inline VectorOfVectors::iterator end()
  {
    return vect.end();
  };

  /** const iterator equal to vect.begin
      \param a VectorOfVectors::iterator
  */
  inline VectorOfVectors::const_iterator begin() const
  {
    return vect.begin();
  };

  /** const iterator equal to vect.end
      \param a VectorOfVectors::iterator
  */
  inline VectorOfVectors::const_iterator end() const
  {
    return vect.end();
  } ;

  /** get vect, ie all the vectors of the object
   * \return a VectorOfVectors
   */
  inline VectorOfVectors getAllVect() const
  {
    return vect;
  }

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

  /** get the number of Blocks
   *  \return unsigned int
   */
  inline const unsigned int getNumberOfBlocks() const
  {
    return tabIndex->size();
  };

  /** return the array of double values of the vector
   *  \exception SiconosVectorException
   *  \param unsigned int: vector position (only for block vector)
   *  \return double* : the pointer on the array
   */
  double* getArray(unsigned int = 0) const;

  /** sets all the values of the vector to 0.0
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   */
  void zero();

  /** set all values of the vector component to value.
   * \param a double
   */
  void fill(double);

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

  /** put data of the vector into a string
   */
  const std::string toString() const;

  /** return the element vector[i]
   *  \param an unsigned int i
   *  \return a double
   */
  const double getValue(unsigned int) const;

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

  /** return i-eme SiconosVector of vect
   * \param unsigned int: block number
   * \return a SimpleVector
   */
  SimpleVector getVector(unsigned int) const;

  /** return i-eme SiconosVector of vect
   * \param unsigned int: block number
   * \return a pointer to a SiconosVector
   */
  inline SP::SiconosVector getVectorPtr(unsigned int pos)
  {
    return vect.at(pos);
  };

  /** return i-eme SiconosVector of vect
   * \param unsigned int: block number
   * \return a pointer to a SiconosVector
   */
  inline SPC::SiconosVector getVectorPtr(unsigned int pos) const
  {
    return vect.at(pos);
  };

  /** set i-eme SiconosVector of vect (copy)
   * \param unsigned int: block number
   * \param a SiconosVector
   */
  void setVector(unsigned int, const SiconosVector&);

  /** set i-eme SiconosVector of vect (pointer link)
   * \param unsigned int: block number
   * \param a pointer to a SiconosVector
   */
  void setVectorPtr(unsigned int, SP::SiconosVector);

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SP::SiconosVector
   */
  SP::SiconosVector operator [](unsigned int) ;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SP::SiconosVector
   */
  SPC::SiconosVector operator [](unsigned int) const;

  /** get the index tab
   * \return a standard vector of int
   */
  inline Index getTabIndex() const
  {
    return *tabIndex;
  }

  /** get a pointer to the index tab
   * \return Index*
   */
  inline const Index* getTabIndexPtr() const
  {
    return tabIndex;
  }

  /** get an iterator that points to the first element of tabIndex
   * \return an Index::iterator
   */
  inline Index::iterator tabIndexBegin()
  {
    return tabIndex->begin();
  }

  /** get an iterator that points to tabIndex.end()
   * \return an Index::iterator
   */
  inline Index::iterator tabIndexEnd()
  {
    return tabIndex->end();
  }

  /** get an iterator that points to the first element of tabIndex
   * \return an Index::iterator
   */
  inline Index::const_iterator tabIndexBegin() const
  {
    return tabIndex->begin();
  }

  /** get an iterator that points to tabIndex.end()
   * \return an Index::iterator
   */
  inline Index::const_iterator tabIndexEnd() const
  {
    return tabIndex->end();
  }

  /** get the number of the vector that handles element at position "pos"
      \param unsigned int, position of the element
      \return unsigned int number of the searched vector
  */
  unsigned int getNumVectorAtPos(unsigned int) const;

  /** add a part of the input vector (starting from pos. i) to the current vector
   *  \param an unsigned int i (in-out)
   *  \param a SiconosVector (in-out)
   */
  void addSimple(unsigned int&, const SiconosVector&);

  /** subtract a part of the input vector (starting from pos. i) to the current vector
   *  \param an unsigned int i (in-out)
   *  \param a SiconosVector (in-out)
   */
  void subSimple(unsigned int&, const SiconosVector&);

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

  /** Insert a subvector in this vector: allocation and copy
   *  \param SiconosVector& v : the vector to be inserted
   */
  void insert(const SiconosVector&) ;

  /** Insert a pointer to a subvector in this vector: no reallocation nor copy.
   *  \param a pointer to SP::SiconosVector
   */
  void insertPtr(SP::SiconosVector) ;

};

#endif
