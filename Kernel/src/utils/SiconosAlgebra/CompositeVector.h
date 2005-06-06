#ifndef COMPOSITEVECTOR_H
#define COMPOSITEVECTOR_H

#include "NewSiconosVector.h"
#include "SimpleVector.h"

class SimpleVector;

class CompositeVector : public SiconosVector
{
public:

  CompositeVector();

  CompositeVector(const std::string&, const bool&);

  /** \fn CompositeVector(const vector<double> v)
  *  \brief contructor with a vector
  *  \param vector<double> v
  *  \return CompositeVector
  */
  CompositeVector(const std::vector<double>&);

  /** \fn CompositeVector(const SiconosVector& v)
   *  \brief contructor with a SiconosVector
   *  \param SiconosVector& v
   *  \exception SiconosVectorException
   *  \return CompositeVector
   */
  CompositeVector(const SiconosVector&);
  CompositeVector(const CompositeVector&);

  /** \fn CompositeVector(const int size)
   *  \brief contructor with a size given in parameter. All the elements are initialized to 0.0
   *  \param int size
   *  \exception SiconosVectorException if size < 0
   *  \return CompositeVector
   */
  CompositeVector(const int&);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  ~CompositeVector();

  /********************************************************************************************/

  inline std::vector<SiconosVector*> getSvref()
  {
    return getSvref();
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
  void add(const SiconosVector &v) ;
  void add(SiconosVector *v) ;

  /** \fn void setValues(const vector<double> v)
   *  \brief set the values of the vector to a new set of value
   *  \param vector<double> v
   */
  void setValues(const std::vector<double>& v) ;

  /** \fn int size() const
   *  \brief get the vector size
   *  \exception to be defined
   *  \return int : the vector size
   */
  int size() const  ;

  /** \fn bool read(const std::string& fileName,const std::string& mode = ASCII)
   *  \brief write the vector in a file
   *  \param std::string fileName : the file to read
   *  \param std::string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool read(const std::string &, const std::string& = N_ASCII) ;

  /** \fn bool write(const string& fileName,const string& mode = ASCII)
   *  \brief write the vector in a file
   *  \param string fileName : the file to read
   *  \param string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool write(const std::string& , const std::string& = N_ASCII) const  ;

  /** \fn double* getArray()
  *  \brief return the array of double values of the vector
  *  \exception SiconosVectorException
  *  \return double* : the pointer on the array
  */
  double* getArray() const;


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
  friend SimpleVector operator * (const SiconosMatrix &m, const CompositeVector &v);

  friend SimpleVector matTransVecMult(SiconosMatrix &m, SiconosVector &v);

  //
private:
  std::vector<SiconosVector*> svref;
  std::vector<int> tabindex;

};



#endif // COMPOSITEVECTOR_H
