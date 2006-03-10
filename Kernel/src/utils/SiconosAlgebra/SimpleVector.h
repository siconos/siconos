/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#ifndef __SimpleVector__
#define __SimpleVector__

#include "SiconosVector.h"
#include "SiconosMatrix.h"
#include <string>

class SiconosMatrix;

class SimpleVector : public SiconosVector
{
private:
  LaVectorDouble lavd;

  /** \fn SimpleVector()
   *  \brief default contructor
   *  \return SimpleVector
   */
  SimpleVector();

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** \fn SimpleVector (std::string file, bool ascii)
   *  \brief contructor from an input file
   *  \param a std::string which contains the file path
   *  \param a boolean to indicate if the file is an ascii one
   *  \return SimpleVector
   */
  SimpleVector(const std::string& , const bool&);

  /** \fn SimpleVector(const std::vector<double> v)
   *  \brief contructor with a vector
   *  \param vector<double> v
   *  \return SimpleVector
   */
  SimpleVector(const std::vector<double>&);

  /** \fn SimpleVector(const SiconosVector& v)
   *  \brief copy constructor
   *  \param SimpleVector&
   */
  SimpleVector(const SimpleVector& v);

  /** \fn SimpleVector(const SiconosVector& v)
   *  \brief contructor with a SiconosVector
   *  \param SiconosVector& v
   *  \exception SiconosVectorException
   *  \return SimpleVector
   */
  SimpleVector(const SiconosVector& v);

  /** \fn SimpleVector(const int size)
   *  \brief contructor with a size given in parameter. All the elements are initialized to 0.0
   *  \param int size
   *  \exception SiconosVectorException if size < 0
   *  \return SimpleVector
   */
  SimpleVector(const unsigned int& size);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  ~SimpleVector();

  // --- OTHER FUNCTIONS ---

  /** \fn SimpleVector* getVectorPtr(const unsigned int&) const;
   *  \brief return the current object. This function is really usefull only for composite.
   * \return a pointer to a SimpleVector
   */
  inline SimpleVector* getVectorPtr(const unsigned int&)
  {
    return this;
  };

  /** \fn std::vector<int> getTabIndex() const
   *  \brief get the index tab (usefull only for composite, should not be used for simple) => avoid downcast
   * \return a standard vector of int
   */
  std::vector<unsigned int> getTabIndex() const;

  /** \fn void zero();
   *  \brief set the values to 0.0
   */
  void zero();

  /** \fn std::string toString();
   *  \brief put datas of the vector in a std::string
   */
  std::string toString() const;

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
  double& operator()(const unsigned int &)  const;

  /** \fn void setValue(const int&, const double&)
   *  \brief set the value of one element of the vector
   *  \param double d : the new value
   *  \param int index : the position of the element to be set
   */
  void setValue(const unsigned int&, const double&);

  /** \fn const double getValue(const int&) const
   *  \brief get the value of one element of the vector
   *  \param int index : the position of the element
   */
  const double getValue(const unsigned int&) const;




  /** \fn void setValues(const vector<double> v,  const int& = 0)
   *  \brief set the values of the vector to a new set of value
   *  \param vector<double> v
   *  \param: int, not used for Simple, only for composite
   */
  void setValues(const std::vector<double>& v, const unsigned int& = 0) ;

  /** \fn const LaVectorDouble getValues() const
   *  \brief get lavd vector
   *  \return a LaVectorDouble
   *  \param: int, not used for Simple, only for composite
   */
  inline const LaVectorDouble getValues(const unsigned int& i = 0) const
  {
    return lavd;
  }

  /** \fn void getBlock(const vector<unsigned int>& index,, SimpleVector& block)
   *  \brief get block corresponding to indexes given in index
   *  \param a vector<unsigned int> for indexes and a SimpleVector (in-out paramater)
   */
  void getBlock(const std::vector<unsigned int>& , SimpleVector&) const;

  /** \fn unsigned int size() const
   *  \brief get the vector size
   *  \exception to be defined
   *  \return int : the vector size
   */
  unsigned int size(const unsigned int& = 0) const  ;

  /** \fn bool read(std::string fileName, std::string mode = ASCII)
   *  \brief read the vector in a file
   *  \param std::string fileName : the file to read
   *  \param std::string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool read(const std::string&, const std::string& = DEFAULT_FORMAT) ;

  /** \fn bool write(std::string fileName, std::string mode = ASCII)
   *  \brief write the vector in a file
   *  \param std::string fileName : the file to read
   *  \param std::string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool write(const std::string& , const std::string& = DEFAULT_FORMAT) const  ;

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
  inline void add(const SiconosVector &v)
  {
    SiconosVectorException::selfThrow(" SimpleVector::add function should not be used for SimpleVector");
  }

  /** \fn bool addPtr(SiconosVector* v)
   *  \brief add a sub Vector in this vector
   *  \param SiconosVector*
   *  \exception SiconosVectorException
   */
  inline void addPtr(SiconosVector *v)
  {
    SiconosVectorException::selfThrow(" SimpleVector::add function should not be used for SimpleVector");
  }

  /** \fn double norm()
   *  \brief return the Euclidian norm of the vector
   *  \return a double
   */
  double norm() const ;

  // internal generic operators
  SimpleVector &operator+=(const SiconosVector &) ;
  SimpleVector &operator-=(const SiconosVector &) ;
  SimpleVector &operator = (const SiconosVector& v) ;
  // internal specific operators
  SimpleVector &operator+=(const SimpleVector &) ;
  SimpleVector &operator-=(const SimpleVector &) ;
  SimpleVector &operator*=(const double&) ;
  SimpleVector &operator/=(const double&) ;
  SimpleVector &operator = (const SimpleVector& v) ;

  /** \fn bool operator == (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  bool operator == (const SiconosVector& v) const  ;
  bool operator == (const SimpleVector& v) const  ;

  /** \fn bool operator != (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */
  bool operator != (const SiconosVector& v) const  ;
  bool operator != (const SimpleVector& v) const  ;


  // generic internal operator for mixed operations
  SimpleVector addition(const SiconosVector&) const;
  SimpleVector subtraction(const SiconosVector&) const;

  // generic external operators
  friend SimpleVector operator + (const SiconosVector& v1, const SiconosVector& v2);
  friend SimpleVector operator - (const SiconosVector& v1, const SiconosVector& v2);

  // specific external operators
  friend SimpleVector operator * (const SimpleVector& v, const double& d) ;
  friend SimpleVector operator * (const double& d, const SimpleVector& v);
  friend SimpleVector operator / (const SimpleVector&  v, const double& d);
  //friend SimpleVector operator + (const SimpleVector& v1, const SimpleVector& v2);
  //friend SimpleVector operator - (const SimpleVector& v1, const SimpleVector& v2);
  friend SimpleVector operator * (const SiconosMatrix &m, const SimpleVector &v);
  friend SimpleVector operator * (const SiconosMatrix &m, const SiconosVector &v);
  friend SimpleVector operator * (const SimpleVector &v, const SiconosMatrix &m);
  friend SimpleVector operator * (const SiconosVector &v, const SiconosMatrix &m);

  friend SimpleVector matTransVecMult(SiconosMatrix &, SimpleVector &);
  friend SimpleVector matTransVecMult(SiconosMatrix &, SiconosVector &);

};
#endif

