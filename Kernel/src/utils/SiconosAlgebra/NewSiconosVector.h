
#ifndef __NewSiconosVector__
#define __NewSiconosVector__

#include <iostream>
#include <lapack++.h>
#include <fstream>
#include <stdio.h>
#include <vector>
#include "SiconosVectorException.h"
#include "check.h"

//using namespace std;

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


  /** \fn double norm()
   *  \brief return the Euclidian norm of the vector
   *  \return a double
   */
  //virtual double norm() = 0;


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

