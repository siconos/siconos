//$Id: CompositeVector.h,v 1.9 2005/02/11 13:30:40 charlety Exp $

#ifndef COMPOSITEVECTOR_H
#define COMPOSITEVECTOR_H

#include "NewSiconosVector.h"
#include "SimpleVector.h"

class CompositeVector : public SiconosVector
{
public:

  CompositeVector();

  CompositeVector(const string file, const bool ascii);

  /** \fn CompositeVector(const vector<double> v)
  *  \brief contructor with a vector
  *  \param vector<double> v
  *  \return CompositeVector
  */
  CompositeVector(const vector<double> v);

  /** \fn CompositeVector(const SiconosVector& v)
   *  \brief contructor with a SiconosVector
   *  \param SiconosVector& v
   *  \exception SiconosVectorException
   *  \return CompositeVector
   */
  CompositeVector(const SiconosVector& v);
  CompositeVector(const CompositeVector& v);

  /** \fn CompositeVector(const int size)
   *  \brief contructor with a size given in parameter. All the elements are initialized to 0.0
   *  \param int size
   *  \exception SiconosVectorException if size < 0
   *  \return CompositeVector
   */
  CompositeVector(const int size);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  ~CompositeVector();

  /********************************************************************************************/

  inline vector<SiconosVector*> getSvref()
  {
    return this->getSvref();
  };

  /** \fn void display();
     *  \brief display data on standard output
     */
  void display() const  ;

  /** \fn operator (int index)
   *  \brief set the element vector[i]
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return the element vector[i]
   */
  double& operator()(const int unsigned index) ;

  /** \fn operator (int index)
   *  \brief get the element vector[i]
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return the element vector[i]
   */
  double operator()(const int unsigned index) const  ;


  /** \fn bool add(const SiconosVector& v)
   *  \brief add a sub Vector in this vector
   *  \param SiconosVector& v : the vector to add
   *  \exception SiconosVectorException
   *  \return bool : true if no error
   */
  void add(const SiconosVector& v) ;

  /** \fn void setValues(const vector<double> v)
   *  \brief set the values of the vector to a new set of value
   *  \param vector<double> v
   */
  void setValues(const vector<double> v) ;

  /** \fn int size() const
   *  \brief get the vector size
   *  \exception to be defined
   *  \return int : the vector size
   */
  int size() const  ;

  /** \fn bool read(string fileName, string mode = ASCII)
   *  \brief write the vector in a file
   *  \param string fileName : the file to read
   *  \param string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool read(string fileName, string mode = N_ASCII) ;

  /** \fn bool write(string fileName, string mode = ASCII)
   *  \brief write the vector in a file
   *  \param string fileName : the file to read
   *  \param string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool write(string fileName, string mode = N_ASCII) const  ;

  /** \fn double* getArray()
  *  \brief return the array of double values of the vector
  *  \exception SiconosVectorException
  *  \return double* : the pointer on the array
  */
  double* getArray();


  // generic internal operators
  CompositeVector &operator+=(const SiconosVector &) ;
  CompositeVector &operator-=(const SiconosVector &) ;
  CompositeVector &operator = (const SiconosVector& v) ;

  // Logical operators
  /** \fn bool operator == (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  bool operator == (const SiconosVector& v) const  ;

  /** \fn bool operator != (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  bool operator != (const SiconosVector& v) const  ;

  // generic internal operator for mixed operations
  CompositeVector addition(const SiconosVector&) const;
  CompositeVector subtraction(const SiconosVector&) const;

  // specific internal operators
  CompositeVector &operator+=(const CompositeVector &) ;
  CompositeVector &operator-=(const CompositeVector &) ;
  CompositeVector &operator*=(const double) ;
  CompositeVector &operator/=(const double) ;
  CompositeVector &operator = (const CompositeVector& v) ;

  // generic external operators
  //  friend CompositeVector operator + (const CompositeVector& v1, /*const*/ SiconosVector& v2);
  //  friend CompositeVector operator - (const CompositeVector& v1, const SiconosVector& v2);


  // specific external operators
  friend CompositeVector operator * (const CompositeVector& v, const double d) ;
  friend CompositeVector operator * (const double d, const CompositeVector& v);
  friend CompositeVector operator / (const CompositeVector&  v, const double d);
  friend CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2);
  friend CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2);
  friend SimpleVector operator * (/*const*/ SiconosMatrix &m, /*const*/ CompositeVector &v);

  friend SimpleVector matTransVecMult(SiconosMatrix &m, SiconosVector &v);

  //
private:
  vector<SiconosVector*> svref;
  vector<int> tabindex;

};



#endif // COMPOSITEVECTOR_H
//$Log: CompositeVector.h,v $
//Revision 1.9  2005/02/11 13:30:40  charlety
//_ added or modified some doxygen comments in SiconosMatrix.
//_ the function "getArray" is vector is now pure virtual, and implemented in CompositeVector (it returns NULL if the vector is composite).
//
//Revision 1.8  2004/09/09 14:32:45  charlety
//
//_ New tests for operators of multiplication between vectors and matrices.
//
//Revision 1.7  2004/08/24 11:29:20  charlety
//
//_ methods replacing generic operators for mixed operations done.
//
//Revision 1.6  2004/08/24 07:35:07  charlety
//
//_ CompositeVector finished at 95%.
//
//TODO :
//_ solve the problem of operators (gcc cannot choice between theoperators of simpleVector and compositeVector in mixed operations).
//_ the behavior of specific operator= of CompositeVector is not ok. Should copy values only.
//
//Revision 1.5  2004/08/20 12:32:25  charlety
//
//_ operators for CompositeVector, read / write functions in file.
//
//Revision 1.4  2004/08/19 15:21:27  charlety
//
//_ SimpleVector and CompositeVector in progress.
//_ for the operators, we prefer now using directly functions of Blas1++ instead
//  of these of Blas++.h
//
//Revision 1.3  2004/08/18 14:53:21  charlety
//
//_ use of Lapack routines for operations on SimpleVector (in progress)
//