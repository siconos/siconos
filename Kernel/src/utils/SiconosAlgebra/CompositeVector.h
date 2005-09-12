#ifndef COMPOSITEVECTOR_H
#define COMPOSITEVECTOR_H

#include "NewSiconosVector.h"
#include "SimpleVector.h"

class SimpleVector;

/** \class CompositeVector
 *  \brief This class describes Composite Vectors, containers of several Siconos Vectors (that should be SimpleVectors)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *
 */

class CompositeVector : public SiconosVector
{
public:

  // CONSTRUCTORS
  /** \fn CompositeVector()
   *  \brief default contructor
   */
  CompositeVector();

  /** \fn CompositeVector(const std::string&, const bool&)
   *  \brief contructor from data by function read call
   *  \param a string
   *  \param a bool
   */
  CompositeVector(const std::string&, const bool&);

  /** \fn CompositeVector(const SiconosVector& v)
   *  \brief contructor with a SiconosVector
   *  \param SiconosVector& v
   *  \exception SiconosVectorException
   */
  CompositeVector(const SiconosVector&);

  /** \fn CompositeVector(const SiconosVector& v)
   *  \brief contructor with a SimpleVector
   *  \param SimpleVector& v
   */
  CompositeVector(const SimpleVector& v);

  /** \fn CompositeVector(const SiconosVector& v)
   *  \brief copy contructor
   *  \param SiconosVector& v
   */
  CompositeVector(const CompositeVector&);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  ~CompositeVector();

  // GETTERS/SETTERS
  /** \fn std::vector<SiconosVector*> getSvref() const
   *  \brief get svref
   * \return a standard vector ofm SiconosVector
   */
  inline std::vector<SiconosVector*> getSvref() const
  {
    return svref;
  }

  /** \fn std::vector<int> getTabIndex() const
   *  \brief get the index tab
   * \return a standard vector of int
   */
  inline std::vector<int> getTabIndex() const
  {
    return tabindex;
  }

  /** \fn std::string toString();
   * \brief put datas of the vector in a std::string
   * useless for composite ?
   */
  inline std::string toString() const
  {
    return "composite!!";
  }

  /** \fn void display();
   *  \brief display data on standard output
   */
  void display() const  ;

  /** \fn operator (int index)
   *  \brief get the element vector[i]
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return the element vector[i]
   */
  double& operator()(const unsigned int&) const ;

  /** \fn void setValue(const int unsigned index, const double d)
   *  \brief set the value of one element of the vector
   *  \param double d : the new value
   *  \param int index : the position of the element which is set
   */
  inline void setValue(const int unsigned& i, const double& d)
  {
    (*this)(i) = d;
  }

  /** \fn double getValue(const int unsigned i)
   *  \brief get the value of index i element of the vector
   *  \param int index : the position of the element
   */
  inline const double getValue(const int unsigned& index) const
  {
    return (*this)(index);
  }

  /** \fn void setValues(const vector<double> v)
  *  \brief set the values of the vector to a new set of value
  *  \param vector<double> v
  *  \param optional, the index of required vector in svref
  */
  void setValues(const std::vector<double>& v, const int unsigned& = 0) ;

  /** \fn const LaVectorDouble getValues() const
   *  \brief get the values saved in vector (depends on vector type)
   *  \param optional, the index of required vector in svref
   *  \return a LaVectorDouble
   */
  const LaVectorDouble getValues(const int unsigned& = 0) const ;

  /** \fn unsigned int size() const
   *  \brief get the vector size
   *  \exception to be defined
   *  \return int : the vector size
   */
  unsigned int size(const unsigned int& = 0) const  ;

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

  /** \fn bool add(const SiconosVector& v)
   *  \brief add a sub Vector in this vector
   *  \param SiconosVector& v : the vector to add
   *  \exception SiconosVectorException
   */
  void add(const SiconosVector &v) ;

  /** \fn bool addPtr(SiconosVector* v)
   *  \brief add a sub Vector in this vector
   *  \param SiconosVector*
   *  \exception SiconosVectorException
   */
  void addPtr(SiconosVector *v) ;

  // generic internal operators
  CompositeVector &operator+=(const SiconosVector &) ;
  CompositeVector &operator-=(const SiconosVector &) ;
  CompositeVector &operator = (const SiconosVector& v) ;
  // specific internal operators
  CompositeVector &operator+=(const CompositeVector &) ;
  CompositeVector &operator-=(const CompositeVector &) ;
  CompositeVector &operator*=(const double&) ;
  CompositeVector &operator/=(const double&) ;
  CompositeVector &operator = (const CompositeVector& v) ;

  // Logical operators
  /** \fn bool operator == (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  bool operator == (const SiconosVector& v) const  ;
  bool operator == (const CompositeVector& v) const  ;

  /** \fn bool operator != (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  bool operator != (const SiconosVector& v) const  ;
  bool operator != (const CompositeVector& v) const  ;

  // generic internal operator for mixed operations
  CompositeVector addition(const SiconosVector&) const;
  CompositeVector subtraction(const SiconosVector&) const;


  // generic external operators
  //  friend CompositeVector operator + (const CompositeVector& v1, /*const*/ SiconosVector& v2);
  //  friend CompositeVector operator - (const CompositeVector& v1, const SiconosVector& v2);

  // specific external operators
  friend CompositeVector operator * (const CompositeVector&, const double&) ;
  friend CompositeVector operator * (const double&, const CompositeVector&);
  friend CompositeVector operator / (const CompositeVector&, const double&);
  friend CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2);
  friend CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2);
  friend SimpleVector operator * (const SiconosMatrix &m, const CompositeVector &v);

  friend SimpleVector matTransVecMult(SiconosMatrix &, CompositeVector &);


  //
private:
  // A container of pointers on SiconosVector (that are to be SimpleVector : no Composite of Composite allowed
  std::vector<SiconosVector*> svref;
  //
  std::vector<int> tabindex;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  std::vector<bool> isSvrefAllocatedIn;
};

#endif
