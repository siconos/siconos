/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

/** \class MySiconosVector
 *  \brief This is an abstract class to provide interface for vector handling
 *  vector can be either a MySimpleVector or a MyBlockVector, ie a container of several MySimpleVector
 *  See documentation of these derivated classes for more details
 *  \author SICONOS Development Team - copyright INRIA
 *  \date (creation) 07/21/2006
 *  \version 1.3.0.
 *
 */

#ifndef __MySiconosVector__
#define __MySiconosVector__

#include "SiconosAlgebra.h"
#include "SiconosVectorException.h"

/**\brief MyVect is an union of DenseVect pointer and SparseVect pointer
 */
union MyVect
{
  DenseVect *Dense; // num = 1
  SparseVect *Sparse; // num = 4
};

class MySiconosVector
{
protected:

  /**\var isBlockVector
   * bool to check the type of the current vector; true if block else false.
   */
  bool isBlockVector;

  /**\fn MySiconosVector (bool = false)
   * \brief default constructor */
  MySiconosVector(bool = false);

public:

  /**\fn bool isBlock ()
   * \brief true if the vector is block else false.
   * \return a bool.
   */
  inline bool isBlock(void) const
  {
    return isBlockVector;
  }

  /**\fn ~MySiconosVector ()
   * \brief Destructor.
   */
  virtual ~MySiconosVector(void) = 0;

  /** \fn unsigned int getNum() const
   *  \brief get the attribute num of current vector
   * \return an unsigned int.
   */
  virtual unsigned int getNum() const = 0;

  /** \fn DenseVect getDense(unsigned int = 0)
   *  \brief get the attribute if it's type is DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect
   */
  virtual const DenseVect getDense(unsigned int = 0) const = 0;

  /** \fn SparseVect getSparse(unsigned int = 0)
   *  \brief get the attribute if it's type is SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect
   */
  virtual const SparseVect getSparse(unsigned int = 0) const = 0;

  /** \fn DenseVect* getDensePtr(unsigned int = 0)
   *  \brief get a pointer on DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect*
   */
  virtual DenseVect* getDensePtr(unsigned int = 0) const = 0;

  /** \fn SparseVect* getSparsePtr(unsigned int = 0)
   *  \brief get a pointer on SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect*
   */
  virtual SparseVect* getSparsePtr(unsigned int = 0) const = 0;

  /** \fn void zero();
   *  \brief sets all the values of the vector to 0.0
   */
  virtual void zero() = 0;

  /** \fn unsigned int size() const
   *  \brief get the vector size, ie the total number of (double)
   *  elements in the vector
   *  \return unsigned int
   */
  virtual unsigned int size() const = 0;

  /** \fn  void resize (unsigned int nbcol, bool val = true)const
   *  \brief resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   *  \param: dim of the resized vector
   *  \param: a bool, true if old values are preserved else false. Default = true.
   *  \exception for Block Vector, resizing not allowed.
   */
  virtual void resize(unsigned int, bool = true) = 0;

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the vector
   *  \return a double
   */
  virtual const double normInf() const = 0;

  /** \fn double norm()
   *  \brief return the Euclidian norm of the vector
   *  \return a double
   */
  virtual const double norm() const = 0 ;

  /** \fn void display();
   *  \brief display data on standard output
   */
  virtual void display(void)const = 0;

  /** \fn MySiconosVector* getVectorPtr(unsigned int);
   *  \brief if this is a block vector return i-eme MySimpleVector, else return this.
   * \return a pointer to a SimpleVector
   */
  virtual MySiconosVector* getVectorPtr(unsigned int) = 0;

  /** \fn void fill(double value);
   *  \brief set all values of the vector component to value.
   * \param a double
   */
  virtual void fill(double) = 0;

  /** \fn unsigned int getNumberOfBlocks() const
   *  \brief get the number of SimpleVector-Blocks - only usefull for BlockVector.
   *  \return unsigned int
   */
  virtual unsigned int getNumberOfBlocks() const
  {
    return 1;
  };

  // Note: in the following functions, index is a general one;
  // that means that for a SimpleVector v, v(i) is index i element but
  // for a BlockVector w that contains 2 SiconosVector of size 3
  // w(4) corresponds to the first element of the second vector.

  /** \fn double& operator ()(unsigned int i)
   *  \brief get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  virtual double& operator()(unsigned int) = 0;

  /** \fn double operator ()(unsigned int i)const
   *  \brief get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  virtual double operator()(unsigned int) const = 0;

  /** \fn MySiconosVector* operator[] (unsigned int i)
   *  \brief get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a MySiconosVector*
   */
  virtual MySiconosVector* operator [](unsigned int) = 0;

  /** \fn MySiconosVector* operator[] (unsigned int i)
   *  \brief get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a MySiconosVector*
   */
  virtual const MySiconosVector* operator [](unsigned int) const = 0;

  /** \fn operator = (const MySiconosVector&)
   *  \param MySiconosVector : the vector to be copied
   */
  virtual MySiconosVector& operator =(const MySiconosVector&) = 0;

  /** \fn operator += (const MySiconosVector&)
   *  \param MySiconosVector : a vector to add
   */
  virtual MySiconosVector& operator +=(const MySiconosVector&) = 0;

  /** \fn operator -= (const MySiconosVector&)
   *  \param MySiconosVector : a vector to subtract
   */
  virtual MySiconosVector& operator -=(const MySiconosVector&) = 0;

  /** \fn operator /= (double)
   *  \param double, a scalar
   */
  virtual MySiconosVector& operator /=(double) = 0;

  /** \fn operator /= (int)
   *  \param int, a scalar
   */
  virtual MySiconosVector& operator /=(int) = 0;

  /** \fn operator*=(double)
   *  \brief multiply the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   */
  virtual MySiconosVector& operator *=(double) = 0;

  /** \fn operator*=(int)
   *  \brief multiply the current vector with an int
   *  \param an int
   *  \exception SiconosVectorException
   */
  virtual MySiconosVector& operator *=(int) = 0;
};

#endif
