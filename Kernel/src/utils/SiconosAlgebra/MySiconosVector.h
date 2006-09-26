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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "SiconosVectorException.h"


using namespace boost::numeric::ublas;

/** \brief DenseVect is a define of boost::ublas::numeric::vector<double, std::vector<double> >
 */
#define DenseVect  vector<double,std::vector<double> >
/** \brief SparseVect is a define of boost::ublas::numeric::mapped<double>
 */
#define SparseVect mapped_vector<double>
const double tolerance = 1e-10; // value used to compare matrices. Matrices A and B are equal when (A-B).normInf()<tolerance.

/** \brief TYP is an enumerated type of DENSE, TRIANGULAR, SYMMETRIC, SPARSE, BANDED. TYP is used to describe the type of matrix or vector we want to construct.
 */
enum TYP {DENSE = 1, TRIANGULAR, SYMMETRIC, SPARSE, BANDED};

/**\brief MyVect is an union of DenseVect pointer and SparseVect pointer
 */
union MyVect
{
  DenseVect *Dense;
  SparseVect *Sparse;
};

class MySiconosVector
{
protected:

  /**\var isBlockVector
   * bool to check the type of the current vector; true if block else false.
   */
  bool isBlockVector;

  /**\fn setIsBlock (bool)
   * \brief set the value of isBlockVector.
   * \param a bool to set isBlockVector (optional, default = false)
   */
  void setIsBlock(bool val)
  {
    isBlockVector = val;
  }

public:

  /**\fn bool isBlock ()
   * \brief true if the vector is block else false.
   * \return a bool.
   */
  inline  bool isBlock(void) const
  {
    return isBlockVector;
  }

  /**\fn ~MySiconosVector ()
   * \brief Destructor.
   */
  virtual ~MySiconosVector(void) = 0;

  /** \fn int getNum() const
   *  \brief get the attribute num of current vector
   * \return an int.
   */
  virtual int getNum(void)const = 0;

  /** \fn DenseVect getDense()
   *  \brief get the attribute if it's type is DenseVect
   *  \return a DenseVect
   */
  virtual const DenseVect getDense(void)const = 0;

  /** \fn SparseVect getSparse()
   *  \brief get the attribute if it's type is SparseVect
   *  \return a SparseVect
   */
  virtual const SparseVect getSparse(void)const = 0;

  /** \fn DenseVect* getDensePtr()
   *  \brief get a pointer on DenseVect
   *  \return a DenseVect*
   */
  virtual const DenseVect* getDensePtr(void)const = 0;

  /** \fn SparseVect* getSparsePtr()
   *  \brief get a pointer on SparseVect
   *  \return a SparseVect*
   */
  virtual const SparseVect* getSparsePtr(void)const = 0;

  /** \fn void setNum(int n)
   *  \brief set the attribute num of current vector with n
   */
  virtual void setNum(int) = 0;

  /** \fn void zero();
   *  \brief sets all the values of the vector to 0.0
   */
  virtual void zero(void) = 0;

  /** \fn unsigned int size() const
   *  \brief get the vector size, ie the total number of (double)
   *  elements in the vector
   *  \return int
   */
  virtual int size(void)const = 0;

  /** \fn  void resize (int nbcol, bool val = true)const
   *  \brief resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   *  \exception SiconosVectorException
   */
  virtual void resize(int, bool = true) = 0;

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the vector
   *  \return a double
   */
  virtual const double normInf(void)const = 0;

  /** \fn void display();
   *  \brief display data on standard output
   */
  virtual void display(void)const = 0;

  // Note: in the following functions, index is a general one;
  // that means that for a SimpleVector v, v(i) is index i element but
  // for a BlockVector w that contains 2 SiconosVector of size 3
  // w(4) corresponds to the first element of the second vector.

  /** \fn double& operator ()(int i)
   *  \brief get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  virtual double& operator()(int) = 0;

  /** \fn double operator ()(int i)const
   *  \brief get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  virtual double operator()(int)const = 0;

  /** \fn operator * (const MySiconosVector &v)
   *  \brief multiply the current vector with the vector v
   *  \param a MySiconosVector
   *  \return a MySiconosVector, the result of the multiplication
   */
  virtual double operator * (const MySiconosVector&) = 0;

  /** \fn operator = (const MySiconosVector&)
   *  \param MySiconosVector : the vector to be copied
   */
  virtual const MySiconosVector& operator =(const MySiconosVector&) = 0;

  /** \fn operator += (const MySiconosVector&)
   *  \param MySiconosVector : a vector to add
   */
  virtual const MySiconosVector& operator +=(const MySiconosVector&) = 0;

  /** \fn operator -= (const MySiconosVector&)
   *  \param MySiconosVector : a vector to subtract
   */
  virtual const MySiconosVector& operator -=(const MySiconosVector&) = 0;

  /** \fn operator /= (double)
   *  \param double, a scalar
   */
  virtual const MySiconosVector& operator /=(double) = 0;

  /** \fn operator /= (int)
   *  \param int, a scalar
   */
  virtual const MySiconosVector& operator /=(int) = 0;

  /** \fn operator*=(double)
   *  \brief multiply the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   */
  virtual const MySiconosVector& operator *=(double) = 0;

  /** \fn operator*=(int)
   *  \brief multiply the current vector with an int
   *  \param an int
   *  \exception SiconosVectorException
   */
  virtual const MySiconosVector& operator *=(int) = 0;
};

#endif
