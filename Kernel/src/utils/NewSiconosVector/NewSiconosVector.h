//$Id: NewSiconosVector.h,v 1.13 2005/02/15 15:15:33 charlety Exp $

#ifndef __NewSiconosVector__
#define __NewSiconosVector__

#include <iostream>
#include <lapack++.h>
#include <fstream>
#include <stdio.h>
#include <vector>
#include "SiconosVectorException.h"
#include "check.h"

using namespace std;

const string N_ASCII = "ascii";
const string N_BINARY = "binary";
const char N_DOUBLE_PRECISION[] = "%1.52e "; // double mantisse precision /!\ MACHINE DEPENDE

const short M_MAXSIZEFORDISPLAY = 10;


class SiconosVector
{
protected:

  bool composite;

public:

  /** \fn SiconosVector()
   *  \brief contructor
   *  \return SiconosVector
   */
  SiconosVector();

  /** \fn SiconosVector(const vector<double> v)
   *  \brief contructor with a vector
   *  \param vector<double> v
   *  \return SiconosVector
   */
  SiconosVector(const vector<double> v);


  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  virtual ~SiconosVector();

  /***********************************************************************************************/

  inline bool isComposite() const
  {
    return this->composite;
  }

  /** \fn void zero();
   *  \brief set the values to 0.0
   */
  virtual void zero();


  /** \fn string toString();
  *  \brief put datas of the vector in a string
  */
  string toString();

  /** \fn void display();
   *  \brief display data on standard output
   */
  virtual void display() const = 0 ;

  /** \fn operator (int index)
   *  \brief set the element vector[i]
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return the element vector[i]
   */
  virtual double& operator()(const int unsigned index) = 0;

  /** \fn operator (int index)
   *  \brief get the element vector[i]
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return the element vector[i]
   */
  virtual double operator()(const int unsigned index) const = 0 ;

  /** \fn void setValue(const int unsigned index, const double d)
   *  \brief set the value of one element of the vector
   *  \param double d : the new value
   *  \param int index : the position of the element which is set
   */
  inline void setValue(const int unsigned index, const double d)
  {
    (*this)(index) = d;
  }

  /** \fn double getValue(const int unsigned index)
   *  \brief get the value of one element of the vector
   *  \param int index : the position of the element
   */
  inline double getValue(const int unsigned index) const
  {
    return (*this)(index);
  }

  /** \fn void setValues(const vector<double> v)
   *  \brief set the values of the vector to a new set of value
   *  \param vector<double> v
   */
  virtual void setValues(const vector<double> v) = 0;

  /** \fn int size() const
   *  \brief get the vector size
   *  \exception to be defined
   *  \return int : the vector size
   */
  virtual int size() const = 0 ;

  /** \fn bool read(string fileName, string mode = ASCII)
   *  \brief write the vector in a file
   *  \param string fileName : the file to read
   *  \param string mode : ASCII or BINARY
   *  \exception SiconosVectorException
   *  \return true if no error
   */
  virtual bool read(string fileName, string mode = N_ASCII) = 0;

  /** \fn bool write(string fileName, string mode = ASCII)
   *  \brief write the vector in a file
   *  \param string fileName : the file to read
   *  \param string mode : ASCII (default mode) or BINARY
   *  \exception SiconosVectorException
   *  \return true if no error
   */
  virtual bool write(string fileName, string mode = N_ASCII) const = 0 ;

  /** \fn double* getArray()
   *  \brief return the array of double values of the vector
   *  \exception SiconosVectorException
   *  \return double* : the pointer on the array
   */
  virtual double* getArray() = 0;

  // OPERATORS ---------------------------------------------------------------

  // internes
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
  virtual SiconosVector &operator*=(const double) = 0;

  /** \fn operator/=(const double)
   *  \brief divide the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   *  \return the result of the division in the current vector
   */
  virtual SiconosVector &operator/=(const double) = 0;

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

//$Log: NewSiconosVector.h,v $
//Revision 1.13  2005/02/15 15:15:33  charlety
//
//_ modified some very slow functions to increase performance
//
//Revision 1.12  2005/02/11 13:30:39  charlety
//_ added or modified some doxygen comments in SiconosMatrix.
//_ the function "getArray" is vector is now pure virtual, and implemented in CompositeVector (it returns NULL if the vector is composite).
//
//Revision 1.11  2005/02/02 15:54:51  jbarbier
//- sample RollingBalls added
//
//- function getArray() added to SimpleVector to return the pointer on the array of double values
//
//Revision 1.10  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.9  2004/09/22 11:16:29  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.8  2004/09/14 13:24:53  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.7  2004/09/09 14:32:46  charlety
//
//_ New tests for operators of multiplication between vectors and matrices.
//
//Revision 1.6  2004/08/19 15:21:27  charlety
//
//_ SimpleVector and CompositeVector in progress.
//_ for the operators, we prefer now using directly functions of Blas1++ instead
//  of these of Blas++.h
//
//Revision 1.5  2004/08/13 10:36:11  charlety
//
//_ tests for simpleVector in progress
//
//Revision 1.4  2004/08/12 11:15:30  charlety
//
//_ SimpleVector developed at 90% (the functions of computations with matrices remain to do.)
//
//Revision 1.3  2004/08/11 14:16:07  charlety
//
//_ NewSiconosVector in progress...(NewSiconosVector is an abstract class and
//  SimpleVector inherits of NewSiconosVector).
//
//Revision 1.2  2004/07/30 14:21:54  charlety
//
//_ new functions and tests for the new SiconosVector
//
