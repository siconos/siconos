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
#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include "SiconosVector.h"
#include "SimpleVector.h"

class SimpleVector;

/** \class BlockVector
 *  \brief This class describes Block Vectors, containers of several Siconos Vectors (that should be SimpleVectors)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *
 */

class BlockVector : public SiconosVector
{
private:
  // A container of pointers on SiconosVector (that are to be SimpleVector : no Block of Block allowed
  std::vector<SimpleVector*> svref;
  // tabindex[i] = tabindex[i-1] + ni, ni being the size of svref[i].
  std::vector<unsigned int> tabindex;

  /** vector of double that contains the double values pointer by SimpleVectors in svref - Necessary to have contiguous values in memory (for Fortran routines ...) */
  std::vector<double> workVector;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  std::deque<bool> isSvrefAllocatedIn;

public:

  // CONSTRUCTORS
  /** \fn BlockVector()
   *  \brief default contructor
   */
  BlockVector();

  /** \fn BlockVector(const std::string, const bool)
   *  \brief contructor from data by function read call
   *  \param a string
   *  \param a bool
   */
  BlockVector(const std::string, const bool);

  /** \fn BlockVector(const SimpleVector& v)
   *  \brief contructor with a SimpleVector
   *  \param SimpleVector& v
   */
  BlockVector(const SimpleVector&);

  /** \fn BlockVector(vector<SimpleVector*> v)
   *  \brief contructor with a list of SimpleVector*
   *  \param a vector<SimpleVector*>
   */
  BlockVector(std::vector<SimpleVector*>);

  /** \fn BlockVector(SimpleVector* v1, SimpleVector* v2)
   *  \brief contructor with a 2 SimpleVectors
   *  \param SimpleVector* v1
   *  \param SimpleVector* v2
   */
  BlockVector(SimpleVector*, SimpleVector*);

  /** \fn BlockVector(const BlockVector& v)
   *  \brief copy contructor
   *  \param BlockVector& v
   */
  BlockVector(const BlockVector&);

  /** \fn BlockVector(unsigned int i, unsigned int j)
   *  \brief constructor with the number of blocks and their dimension (ie all blocks have the same dim)
   *  \param unsigned int : number of blocks
   *  \param unsigned int : dim of each block
   */
  BlockVector(unsigned int, unsigned int);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  ~BlockVector();

  // GETTERS/SETTERS
  /** \fn std::vector<SimpleVector*> getSvref() const
   *  \brief get svref
   * \return a standard vector of SimpleVector
   */
  inline std::vector<SimpleVector*> getSvref() const
  {
    return svref;
  }

  /** \fn SimpleVector* getVectorPtr(const unsigned int) const;
   *  \brief return i-eme SimpleVector of svref
   * \return a pointer to a SimpleVector
   */
  SimpleVector* getVectorPtr(const unsigned int);

  /** \fn std::vector<unsigned int> getTabIndex() const
   *  \brief get the index tab
   * \return a standard vector of int
   */
  inline std::vector<unsigned int> getTabIndex() const
  {
    return tabindex;
  }

  /** \fn std::string toString();
   * \brief put datas of the vector in a std::string
   * useless for block ?
   */
  inline std::string toString() const
  {
    return "block!!";
  }

  /** \fn void display();
   *  \brief display data on standard output
   */
  void display() const  ;

  /** \fn operator (unsigned int index)
   *  \brief get the element vector[i]
   *  \param an unsigned integer i
   *  \exception SiconosVectorException
   *  \return the element vector[i]
   */
  double& operator()(const unsigned int) const ;

  /** \fn void setValue(const unsigned int index, const double d)
   *  \brief set the value of one element of the vector
   *  \param double d : the new value
   *  \param int index : the position of the element which is set
   */
  inline void setValue(const unsigned int i, const double d)
  {
    (*this)(i) = d;
  }

  /** \fn double getValue(const unsigned int i)
   *  \brief get the value of index i element of the vector
   *  \param unsigned int index : the position of the element
   */
  inline const double getValue(const unsigned int index) const
  {
    return (*this)(index);
  }

  /** \fn void setValues(const vector<double> v)
  *  \brief set the values of the vector to a new set of value
  *  \param vector<double> v
  *  \param optional, the index of required vector in svref
  */
  void setValues(const std::vector<double>& v, const unsigned int = 0) ;

  /** \fn const LaVectorDouble getValues() const
   *  \brief get the values saved in vector (depends on vector type)
   *  \param optional, the index of required vector in svref
   *  \return a LaVectorDouble
   */
  const LaVectorDouble getValues(const unsigned int = 0) const ;

  /** \fn unsigned int size(unsigned int  index) const
   *  \brief if index = 0, get the vector total-size (ie number of elements in vector)
   *         if index = 1, get the number of element in svref.
   *  \exception to be defined
   *  \return int : the vector size
   */
  unsigned int size(const unsigned int = 0) const  ;

  /** \fn bool read(const std::string fileName,const std::string mode = ASCII)
   *  \brief write the vector in a file
   *  \param std::string fileName : the file to read
   *  \param std::string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool read(const std::string , const std::string = DEFAULT_FORMAT) ;

  /** \fn bool write(const string fileName,const string mode = ASCII)
   *  \brief write the vector in a file
   *  \param string fileName : the file to read
   *  \param string mode : ASCII or BINARY
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool write(const std::string , const std::string = DEFAULT_FORMAT) const  ;

  /** \fn double* getArray()
   *  \brief return the array of double values of the vector
   *  \exception SiconosVectorException
   *  \return double* : the pointer on the array
   */
  double* getArray() const;

  /** \fn bool add(const SimpleVector& v)
   *  \brief add a sub Vector in this vector
   *  \param SimpleVector& v : the vector to add
   *  \exception SiconosVectorException
   */
  void add(const  SimpleVector&) ;

  /** \fn bool addPtr(SimpleVector* v)
   *  \brief add a sub Vector in this vector
   *  \param SiconosVector*
   *  \exception SiconosVectorException
   */
  void addPtr(SimpleVector*) ;

  /** \fn void zero();
   *  \brief set the values to 0.0
   */
  void zero();

  // generic internal operators
  BlockVector &operator+=(const SiconosVector &) ;
  BlockVector &operator-=(const SiconosVector &) ;
  BlockVector &operator = (const SiconosVector& v) ;
  // specific internal operators
  BlockVector &operator+=(const BlockVector &) ;
  BlockVector &operator-=(const BlockVector &) ;
  BlockVector &operator*=(const double) ;
  BlockVector &operator/=(const double) ;
  BlockVector &operator = (const BlockVector& v) ;

  /** \fn bool operator != (const SiconosVector& v) const;
   *  \brief compares two vectors (sizes and values).
   *  \return bool
   */

  // generic internal operator for mixed operations
  BlockVector addition(const SiconosVector&) const;
  BlockVector subtraction(const SiconosVector&) const;


  // generic external operators
  //  friend BlockVector operator + (const BlockVector& v1, /*const*/ SiconosVector& v2);
  //  friend BlockVector operator - (const BlockVector& v1, const SiconosVector& v2);

  // specific external operators
  friend BlockVector operator * (const BlockVector&, const double) ;
  friend BlockVector operator * (const double, const BlockVector&);
  friend BlockVector operator / (const BlockVector&, const double);
  friend BlockVector operator + (const BlockVector& v1, const BlockVector& v2);
  friend BlockVector operator - (const BlockVector& v1, const BlockVector& v2);
  friend SimpleVector operator * (const SiconosMatrix &m, const BlockVector &v);

  friend SimpleVector matTransVecMult(SiconosMatrix &, BlockVector &);


  //
};

#endif
