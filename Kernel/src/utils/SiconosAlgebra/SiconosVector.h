/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
#ifndef __SiconosVector__
#define __SiconosVector__

#include "SiconosVectorException.h"
#include "SiconosConst.h"
#include "check.h"
#include <iostream>
#include <lapack++.h>
#include <fstream>
#include <vector>
#include <string>
#include<deque>

const char N_DOUBLE_PRECISION[] = "%1.52e "; // double mantisse precision /!\ MACHINE DEPENDE
const unsigned int M_MAXSIZEFORDISPLAY = 10;
const std::string DEFAULT_FORMAT = "ascii";

/** \class SiconosVector
 *  \brief This is an abstract class to provide interface for vector handling
 *  vector can be either a SimpleVector or a CompositeVector, ie a container of several SimpleVector
 *  See documentation of these derivated classes for more details
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.2.
 *
 */

class SiconosMatrix;

class SiconosVector
{
protected:

  // true if composite
  bool composite;

public:

  /** \fn SiconosVector()
   *  \brief contructor
   *  \return SiconosVector
   */
  SiconosVector();

  /** \fn SiconosVector (const SiconosVector & );
   *  \brief copy contructor
   */
  SiconosVector::SiconosVector(const SiconosVector &);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  virtual ~SiconosVector();

  /***********************************************************************************************/

  /** \fn bool isComposite()
   *  \brief test whether the present vector is composite or not
   * \return a bool
   */
  inline bool isComposite() const
  {
    return composite;
  }

  // !!! WARNING : all the following functions are to be implemented in derivated classes !!!

  /** \fn std::vector<int> getTabIndex() const
   *  \brief get the index tab (usefull only for composite, should not be used for simple) => avoid downcast
   * \return a standard vector of int
   */
  virtual std::vector<int> getTabIndex() const = 0;

  /** \fn std::string toString();
   *  \brief put data of the vector into a string
   */
  virtual std::string toString() const = 0;

  /** \fn void display();
   *  \brief display data on standard output
   */
  virtual void display() const = 0 ;

  // Note: in the following functions, index is a general one;
  // that means that for a SimpleVector v, v(i) is index i element but
  // for a CompositeVector w that contains 2 SiconosVector of size 3
  // w(4) corresponds to the first element of the second vector.

  /** \fn operator (int index)
   *  \brief get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  virtual double& operator()(const int unsigned&) const = 0;

  /** \fn void setValue(const int unsigned index, const double d)
   *  \brief set the value of one element of the vector
   *  \param double d : the new value
   *  \param int index : the position of the element which is set
   */
  virtual void setValue(const int unsigned&, const double&) = 0 ;

  /** \fn double getValue(const int unsigned i)
   *  \brief get the value of index i element of the vector
   *  \param int index : the position of the element
   */
  virtual const double getValue(const int unsigned& index) const = 0;

  /** \fn void setValues(const vector<double> v, const int& = 0)
   *  \brief set the values of the vector to a new set of value
   *  \param vector<double> v
   *  \param optional, only for composite, to set values of vector number i
   */
  virtual void setValues(const std::vector<double>& v, const unsigned int& = 0) = 0;

  /** \fn const LaVectorDouble getValues(const int& i) const
   *  \brief get the values saved in vector (depends on vector type)
   *  \param optional, only for composite, to get values of vector number i
   *  \return a LaVectorDouble
   */
  virtual const LaVectorDouble getValues(const int unsigned& = 0) const = 0;

  /** \fn unsigned int size() const
   *  \brief get the vector size, ie the total number of (double)
   *  elements in the vector
   * \param (optional). =0 -> number of element in the vector
   *                    =1 -> number of subvectors if composite
   *  \return int
   */
  virtual unsigned int size(const unsigned int& = 0) const = 0 ;

  /** \fn bool read(std::string fileName, std::string mode = ASCII)
   *  \brief write the vector in a file
   *  \param std::string fileName : the file to read
   *  \param std::string mode : ASCII or BINARY
   *  \exception SiconosVectorException
   *  \return true if no error
   */
  virtual bool read(const std::string& , const std::string& = DEFAULT_FORMAT) = 0;

  /** \fn bool write(std::string fileName, std::string mode = ASCII)
   *  \brief write the vector in a file
   *  \param std::string fileName : the file to read
   *  \param std::string mode : ASCII (default mode) or BINARY
   *  \exception SiconosVectorException
   *  \return true if no error
   */
  virtual bool write(const std::string& , const std::string& = DEFAULT_FORMAT) const  = 0 ;

  /** \fn double* getArray()
   *  \brief return the array of double values of the vector
   *  \exception SiconosVectorException
   *  \return double* : the pointer on the array
   */
  virtual double* getArray() const = 0;

  /** \fn bool add(const SiconosVector& v)
   *  \brief add a sub Vector in this vector - Usefull only for composite
   *  \param SiconosVector& v : the vector to add
   *  \exception SiconosVectorException
   */
  //virtual void add(const SiconosVector &) = 0 ;

  /** \fn void zero();
   *  \brief set the values to 0.0
   */
  virtual void zero() = 0;

  /** \fn bool addPtr(SiconosVector* v)
   *  \brief add a pointer to sub Vector in this vector - Usefull only for composite
   *  \param SiconosVector*
   *  \exception SiconosVectorException
   */
  //virtual void addPtr(SiconosVector *) = 0 ;

  // OPERATORS ---------------------------------------------------------------

  // internal
  /** \fn operator+=(const SiconosVector &)
   *  \brief add a vector to the current one
   *  \param a SiconosVector
   *  \exception SiconosVectorException
   *  \return the result of the addition in the current vector
   */
  virtual SiconosVector &operator+=(const SiconosVector &) = 0;

  /** \fn operator-=(const SiconosVector &)
   *  \brief subtract a vector to the current one
   *  \param a SiconosVector
   *  \exception SiconosVectorException
   *  \return the result of the subtraction in the current vector
   */
  virtual SiconosVector &operator-=(const SiconosVector &) = 0;

  /** \fn operator*=(const double)
   *  \brief multiply the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   *  \return the result of the multiplication in the current vector
   */
  virtual SiconosVector &operator*=(const double&) = 0;

  /** \fn operator/=(const double)
   *  \brief divide the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   *  \return the result of the division in the current vector
   */
  virtual SiconosVector &operator/=(const double&) = 0;

  // affectation operation
  /** \fn operator=(const SiconosVector& v)
   *  \brief affectation of the operator.
   *  \param SiconosVector&
   *  \exception SiconosVectorException
   *  \return SiconosVector & : current vector, equlas to the parameter
   */
  virtual SiconosVector& operator = (const SiconosVector& v) = 0;

  // Logical operators
  /** \fn bool operator == (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  virtual bool operator == (const SiconosVector& v) const = 0 ;

  /** \fn bool operator != (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  virtual bool operator != (const SiconosVector& v) const = 0 ;

};
#endif

