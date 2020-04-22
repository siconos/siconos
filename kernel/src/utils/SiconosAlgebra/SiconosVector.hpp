/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file SiconosVector.hpp
 */

#ifndef __SiconosVector__
#define __SiconosVector__

#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosVectorFriends.hpp"
#include "SiconosVectorException.hpp"

struct SiconosVectorIterator;
struct SiconosVectorConstIterator;

/** Union to gather all types of ublas vectors used in Siconos */
union VECTOR_UBLAS_TYPE
{
  DenseVect *Dense; // num = 1
  SparseVect *Sparse; // num = 4
};


/** Vectors of double. (Interface to various types of Boost-Ublas vectors).

    Two possible types: Siconos::DENSE (default) and Siconos:SPARSE.

    \rst
    Reference :ref:`siconos_algebra` in Siconos users' guide.
    \endrst

*/
class SiconosVector : public std::enable_shared_from_this<SiconosVector>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosVector);

  bool _dense = true;

  /**
   * Union of pointers to the ublas vector type (dense or sparse)
   */
  VECTOR_UBLAS_TYPE vect;

public:

  /***************************** CONSTRUCTORS ****************************/

  /** Creates a zero-size vector. */
  SiconosVector();

  /** creates a vector, all components set to zero.
      \param row the size of the vector
      \param type the type of vector (dense or sparse)
  */
  SiconosVector(unsigned row, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** creates a vector and initializes its content with a single value
      \param row size of the new vector
      \param val value to initialize its content
      \param type type of vector (dense or sparse)
  */
  SiconosVector(unsigned row, double val, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** creates a dense vector from a copy of a stl vector.
      \param vec vector to be copied
      \param type of the vector (dense or sparse)
  */
  SiconosVector(const std::vector<double>& vec, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** copy constructor
      \param v source vector to be copied
  */
  SiconosVector(const SiconosVector& v);

  /** creates a dense vector, with a copy.
   *  \param v source vector (ublas dense)
   */
  SiconosVector(const DenseVect& v);

  /** creates a sparse vector, with a copy.
   *  \param v source vector (ublas sparse)
   */
  SiconosVector(const SparseVect& v);

  /** creates a vector from data in a file
   *  \param filename file name (possibly with path)
   *  \param is_ascii file format (true if ascii, false if binary)
   */
  SiconosVector(const std::string& filename, bool is_ascii);

  /** constructor from the concatenation of two vectors
   * \param v1 the first vector
   * \param v2 the second vector
   */
  SiconosVector(const SiconosVector& v1, const SiconosVector& v2);

  /** constructor from a BlockVector.
   * explicit to forbid implicit conversion/conversion constructor.
   * \param input source vector
   */
  explicit SiconosVector(const BlockVector& input);//, bool = false);

  /** destructor
   */
  ~SiconosVector();

  /** Copy a the content of a BlockVector into a SiconosVector

      The aim of this function is to be able to handle contiguous memory
      (which is not guaranteed in a BlockVector).

      This could be have been done with a constructor but
      doing so leads to implicit copy-construction in operators call.
   */
  void block2contiguous(const BlockVector & vIn);

  /** get the vector size, ie the total number of (double) elements in the vector
   *  \return unsigned int
   */
  unsigned int size() const;

  /** Get the type number of the current vector.
   * \return an unsigned int
   */
  unsigned int num() const
  {
    if (_dense) return 1;
    else return 4;
  }

  /** get a pointer to the ublas embedded vector if it's type is Dense
   *  \return a DenseVect*
   */
  inline DenseVect* dense() const
  {
    return vect.Dense;
  };

  /** get a pointer to the ublas embedded vector if it's type is Sparse
   *  \return a SparseVect*
   */
  SparseVect* sparse() const;

  /** \return the array of double values of the vector
   */
  double* getArray() const;

  /** sets all the values of the vector to 0.0 */
  void zero();

  /** Resize the vector. The existing elements may be preseved if specified.
   * \param size new size of the vector
   * \param preserve true if the content of the vector must be preserved.
   */
  void resize(unsigned int size, bool preserve= true);

  /** \return the infinite norm of the vector */
  double normInf()const;

  /** \return the Euclidian norm of the vector */
  double norm2() const ;

  /** \return the sum of all elements of the vector */
  double vector_sum() const;

  /** display vector content */
  void display(void) const;

  /** set all values of the vector to input value.
   * \param a input value
   */
  void fill(double a);

  /** \return the content of the vector as a string */
  std::string toString() const;

  /** for iterator interface */
  typedef SiconosVectorIterator iterator;

  /** for iterator interface */
  typedef SiconosVectorConstIterator const_iterator;

  /** \return an iterator pointing to the first element in the vector. */
  iterator begin();

  /** \return an iterator pointing to the first element in the vector. */
  const_iterator begin() const;

  /**  \return an iterator referring to the past-the-end element in the vector container. */
  iterator end();

  /**  \return an iterator referring to the past-the-end element in the vector container. */
  const_iterator end() const;

  /** cast a SiconosVector into a std::vector<double> (performs copy) */
  operator std::vector<double>();

  //************************** VECTORS HANDLING AND OPERATORS *******************************

  /** Get a component of the vector
   *  \param i index of the required component
   *  \return the component value
   */
  double getValue(unsigned int i) const ;

  /** set a component of the vector
   *  \param i index of the required component
   *  \param value of the component
   */
  void setValue(unsigned int i, double value);

  /** get a component of the vector
   *  \param i index of the required component
   *  \return value of the component
   */
  double& operator()(unsigned int i);

  /** get a component of the vector
   *  \param i index of the required component
   *  \return value of the component
   */
  double operator()(unsigned int i) const;

  /** set a sub-block of the current vector
      \param i the beginning of the destination range
      \param v vector to be copied
  */
  void setBlock(unsigned int i, const SiconosVector& v);

  /** copy a part of the vector into another 
      \param vOut destination vector
      \param sizeB number of the elements to copy
      \param startIn the beginning of the range of elements to copy from
      \param startOut the beginning of the destination range
  */
  void toBlock(SiconosVector& vOut, unsigned int sizeB,
               unsigned int startIn, unsigned int startOut) const;

  /** add the input vector to a sub-block of the current vector
      \param i the beginning of the destination range
      \param v the source vector to be added
   */
  void addBlock(unsigned int i, const SiconosVector& v);

  /** subtract the input vector to a sub-block of the current vector
      \param i the beginning of the destination range
      \param v the source vector to be added
   */
  void subBlock(unsigned int i, const SiconosVector& v);

  /** copy the vector into an array
   * \param data the memory where to copy the data
   * \return the number of element written (size of the vector)
   */
  unsigned copyData(double* data) const;

  /** operator =
   * \param v the vector to be copied
   * \return  SiconosVector&
   */
  SiconosVector& operator = (const SiconosVector& v);

  /** operator =
   * \param b the vector to be copied
   * \return  SiconosVector&
   */
  SiconosVector& operator = (const BlockVector& b);

  /** operator =
   * \param v the vector to be copied
   * \return  SiconosVector&
   */
  SiconosVector& operator = (const DenseVect& v);

  /** operator =
   *  \param sp the vector to be copied
   * \return  SiconosVector&
   */
  SiconosVector& operator = (const SparseVect& sp);

  /** operator =
   *  \param d data to put the in vector
   * \return  SiconosVector&
   */
  SiconosVector& operator = (const double* d);

  /** operator +=
   * \param v the vector to add
   * \return  SiconosVector&
   */
  SiconosVector& operator +=(const SiconosVector& v);

  /** operator +=
   * \param v the vector to add
   * \return  SiconosVector&
   */
  SiconosVector& operator +=(const BlockVector& v);

  /** operator -=
   * \param  v the vector to subtract
   * \return  SiconosVector&
   */
  SiconosVector& operator -=(const SiconosVector& v);
  /** operator -=
   * \param  v the vector to subtract
   * \return  SiconosVector&
   */
  SiconosVector& operator -=(const BlockVector& v);

  /** \defgroup SiconosVectorFriends

      List of friend functions of the SiconosVector class

      @{
  */

  /** send data of the vector to an ostream
   * \param os An output stream
   * \param sv a SiconosVector
   * \return The same output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const SiconosVector& sv);

  friend SiconosVector& operator *= (SiconosVector& v, const double& s);

  friend SiconosVector& operator /= (SiconosVector& v, const double& s);

  friend bool operator ==(const SiconosVector&, const SiconosVector&);

  friend SiconosVector operator * (double, const SiconosVector&);

  friend SiconosVector operator * (const SiconosVector&, double);

  friend SiconosVector operator / (const SiconosVector&, double);

  friend SiconosVector operator + (const SiconosVector&, const SiconosVector&);

  friend void add(const SiconosVector&, const SiconosVector&, SiconosVector&);

  friend SiconosVector operator - (const SiconosVector&, const SiconosVector&);

  friend void sub(const SiconosVector&, const SiconosVector&, SiconosVector&);

  friend void axpby(double, const SiconosVector&, double, SiconosVector&);

  friend void axpy(double, const SiconosVector&, SiconosVector&);

  friend double inner_prod(const SiconosVector&, const SiconosVector&);

  friend SimpleMatrix outer_prod(const SiconosVector&, const SiconosVector&);

  friend void scal(double, const SiconosVector&, SiconosVector&, bool);

  friend void subscal(double, const SiconosVector&, SiconosVector&, const Index&, bool);

  friend void cross_product(const SiconosVector&, const SiconosVector&, SiconosVector&);

  friend void abs_wise(const SiconosVector&, SiconosVector&);

  friend void getMax(const SiconosVector&, double &, unsigned int &);

  friend void  getMin(const SiconosVector&, double &, unsigned int &);

  friend struct IsDense;

  friend struct IsSparse;

  friend struct IsBlock;

  friend class TestDense;

  /** End of Friend functions group @} */

  //  temporary workaround, the visitor has to be removed or rework -- xhub
  ACCEPT_NONVIRTUAL_VISITORS();

};

/* functor/predicate used to test vectors type containers such as BlockVector */
class TestDense
{
public:
  bool operator()(SP::SiconosVector input) const
  {
    return input->_dense;
  }
};


#endif
