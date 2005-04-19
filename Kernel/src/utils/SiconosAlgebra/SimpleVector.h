#ifndef __SimpleVector__
#define __SimpleVector__

#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include <string>

//using namespace std;

class SiconosMatrix;

class SimpleVector : public SiconosVector
{
private:
  LaVectorDouble lavd;

public:

  /** \fn SimpleVector()
   *  \brief contructor
   *  \return SimpleVector
   */
  SimpleVector();

  /** \fn SimpleVector (string file, bool ascii)
   *  \brief contructor with an input file
   *  \param a string which contain the file path
   *  \param a boolean to indicate if the file is in ascii
   *  \return SimpleVector
   */
  SimpleVector(const string file, const bool ascii);

  /** \fn SimpleVector(const vector<double> v)
   *  \brief contructor with a vector
   *  \param vector<double> v
   *  \return SimpleVector
   */
  SimpleVector(const vector<double> v);

  /** \fn SimpleVector(const SiconosVector& v)
   *  \brief contructor with a SiconosVector
   *  \param SiconosVector& v
   *  \exception SiconosVectorException
   *  \return SimpleVector
   */
  SimpleVector(const SiconosVector& v);
  SimpleVector(const SimpleVector& v);

  /** \fn SimpleVector(const int size)
   *  \brief contructor with a size given in parameter. All the elements are initialized to 0.0
   *  \param int size
   *  \exception SiconosVectorException if size < 0
   *  \return SimpleVector
   */
  SimpleVector(const int size);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  ~SimpleVector();

  /********************************************************************************************/


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

  /** \fn void zero();
   *  \brief set the values to 0.0
   */
  virtual void zero();

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

  // internal generic operators
  SimpleVector &operator+=(const SiconosVector &) ;
  SimpleVector &operator-=(const SiconosVector &) ;
  SimpleVector &operator = (const SiconosVector& v) ;

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
  SimpleVector addition(const SiconosVector&) const;
  SimpleVector subtraction(const SiconosVector&) const;

  // internal specific operators
  SimpleVector &operator+=(const SimpleVector &) ;
  SimpleVector &operator-=(const SimpleVector &) ;
  SimpleVector &operator*=(const double) ;
  SimpleVector &operator/=(const double) ;
  SimpleVector &operator = (const SimpleVector& v) ;






  // generic external operators
  friend SimpleVector operator + (const SiconosVector& v1, const SiconosVector& v2);
  friend SimpleVector operator - (const SiconosVector& v1, const SiconosVector& v2);

  // specific external operators
  friend SimpleVector operator * (const SimpleVector& v, const double d) ;
  friend SimpleVector operator * (const double d, const SimpleVector& v);
  friend SimpleVector operator / (const SimpleVector&  v, const double d);
  friend SimpleVector operator + (const SimpleVector& v1, const SimpleVector& v2);
  friend SimpleVector operator - (const SimpleVector& v1, const SimpleVector& v2);
  friend SimpleVector operator * (const SiconosMatrix &m, const SimpleVector &v);

  friend SimpleVector matTransVecMult(SiconosMatrix &m, SimpleVector &v);

};
#endif

