/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/** Union of DenseVect and SparseVect pointers -
    Siconos::DENSE, num = 1,
    SPARSE, num = 4
 */
union VECTOR_UBLAS_TYPE
{
  DenseVect *Dense; // num = 1
  SparseVect *Sparse; // num = 4
};

/** Vectors of double. (Interface to various types of Boost-Ublas vectors).
 *
 * \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 07/21/2006
 *
 * Used to handle vectors of double.
 *
 * Two possible types: Siconos::DENSE (default) and SPARSE.
 *
 * You can find an overview on how to build and use vectors and matrices in siconos users' guide .
 *
 */
class SiconosVector : public std11::enable_shared_from_this<SiconosVector>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosVector);


  bool _dense;

  /**
   * Union of pointers to the ublas vector type (dense or sparse)
   */
  VECTOR_UBLAS_TYPE vect;

public:
  /** Default constructor
   */
  SiconosVector();

  /***************************** CONSTRUCTORS ****************************/

  /** constructor with the type and the dimension of the Boost vector
   *  \param row the size of the vector
   *  \param type the type of vector
   */
  SiconosVector(unsigned row, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** constructor with the dimension of the Boost vector, a default value and the type.
   *  \param row the size of the vector
   *  \param val value for initiliazing the vector
   *  \param type the type of vector
   */
  SiconosVector(unsigned row, double val, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** constructor with  a std::vector of the values and the type of the boost vector.
   *  \param vec a std::vector<double>
   *  \param type  an Siconos::UBLAS_TYPE
   */
  SiconosVector(const std::vector<double>& vec, Siconos::UBLAS_TYPE type = Siconos::DENSE);

  /** copy constructor
   *  \param v SiconosVector
   */
  SiconosVector(const SiconosVector& v);

  SiconosVector(const BlockVector & vIn);
  /** constructor with a DenseVect vector (see SiconosAlgebra.hpp for details)
   *  \param v a DenseVect
   */
  SiconosVector(const DenseVect& v);

  /** constructor with a SparseVect vector (see SiconosAlgebra.hpp for details)
   *  \param v a SparseVect
   */
  SiconosVector(const SparseVect& v);

  /** constructor with an input file
   *  \param filename a std::string which contain the file path
   *  \param is_ascii a boolean to indicate if the file is in ascii
   */
  SiconosVector(const std::string& filename, bool is_ascii);

  /** constructor for the concatenation of two vectors
   * \param v1 the first vector
   * \param v2 the second vector
   */
  SiconosVector(const SiconosVector& v1, const SiconosVector& v2);

  /** destructor
   */
  ~SiconosVector();

  /** get the vector size, ie the total number of (double) elements in the vector
   *  \return unsigned int
   */
  unsigned int size() const;

  /** true if the vector is block else false.
   * \return a bool.
   */
  bool isBlock() const
  {
    return false;
  };


  /** Get the type number of the current vector.
   * \return an unsigned int
   */
  unsigned int num() const
  {
    if (_dense) return 1;
    else return 4;
  }

  /** get the ublas embedded vector if it's type is Dense
   *  \param pos unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a DenseVect
   */
  const DenseVect getDense(unsigned int pos = 0) const;

  /** get the ublas embedded vector if it's type is Sparse
   *  \param pos unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a SparseVect
   */
  const SparseVect getSparse(unsigned int pos = 0) const;

  /** get a pointer to the ublas embedded vector if it's type is Dense
   *  \param pos unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a DenseVect*
   */
  inline DenseVect* dense(unsigned int pos = 0) const
  {
    return vect.Dense;
  };

  /** get a pointer to the ublas embedded vector if it's type is Sparse
   *  \param pos unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a SparseVect*
   */
  SparseVect* sparse(unsigned int pos = 0) const;

  /** return the array of double values of the vector
   *  \return a pointer to the array
   */
  double* getArray() const;

  /** sets all the values of the vector to 0.0
   */
  void zero();

  /** resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   * \param nbcol
   * \param preserve
   * \exception SiconosVectorException
   */
  void resize(unsigned int nbcol, bool preserve= true);

  /** compute the infinite norm of the vector
   *  \return a double
   */
  double normInf()const;

  /** return the Euclidian norm of the vector
   *  \return a double
   */
  double norm2() const ;

  /** return the sum of all elements of the vector
   * \return a double
   */
  double vector_sum() const;

  /** display data on standard output
   */
  void display(void) const;

  /** return the current object. This function is really useful only for block vector
   * \return a pointer to a SiconosVector
   * \param pos
   */
  inline SP::SiconosVector vector(unsigned int pos)
  {
    return shared_from_this();
  };

  /** return the current object. This function is really useful only for block vector
   * \param pos
   * \return a pointer to a SiconosVector
   */
  inline SPC::SiconosVector vector(unsigned int pos) const
  {
    return shared_from_this();
  };

  /** set SiconosVector number i (copy) with v (second arg) - Useful only for BlockVector (else equivalent to a single copy)
   * \param num unsigned int: block number (0 for SiconosVector)
   * \param v a SiconosVector
   */
  void setVector(unsigned int num, const SiconosVector& v);

  /** set SiconosVector number i (pointer link) with v (second arg) - Useful only for BlockVector
   * \param num unsigned int: block number (0 for SiconosVector)
   * \param v a pointer to a SiconosVector
   */
  inline void setVectorPtr(unsigned int num, SP::SiconosVector v)
  {
    SiconosVectorException::selfThrow("SiconosVector::setVectorPtr(num,v), not allowed for SiconosVector.");
  };

  /** set all values of the vector component to input value.
   * \param a double
   */
  void fill(double a);

  /** put data of the vector into a std::string
   * \return std::string
   */
  const std::string toString() const;

  //************************** VECTORS HANDLING AND OPERATORS *******************************

  /** get the element at position i in the vector
   *  \param i an unsigned int
   *  \return a double
   */
  double getValue(unsigned int i) const ;

  /** set the element at position i in the vector.
   *  \param i an unsigned int
   *  \param value
   */
  void setValue(unsigned int i, double value);

  /** get the element at position i in the vector
   *  \param i an integer
   *  \return a double
   */
  double& operator()(unsigned int i);

  /** get the element at position i in the vector
   *  \param i an integer
   *  \return a double
   */
  double operator()(unsigned int i)const;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param i an unsigned integer
   *  \return a SP::SiconosVector
   */
  inline SP::SiconosVector operator [](unsigned int i)
  {
    return shared_from_this();
  };

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param i an unsigned integer
   *  \return a SP::SiconosVector
   */
  inline SPC::SiconosVector operator [](unsigned int i) const
  {
    return shared_from_this();
  };

  /** set the elements starting from position i with input vector
   *  \param i an unsigned int
   *  \param v a ref to a SiconosVector
   */
  void setBlock(unsigned int i, const SiconosVector& v);

  /** copy a part of the vector into another (usefull for converting a SiconosVector
   * into a BlockVector
   * \param vOut
   * \param sizeB
   * \param startIn
   * \param startOut
   */
  void toBlock(SiconosVector& vOut, unsigned int sizeB,
               unsigned int startIn, unsigned int startOut) const;

  /** add the input vector to the elements starting from position i.
   *  \param i an unsigned int
   *  \param v a SiconosVector
   */
  void addBlock(unsigned int i, const SiconosVector& v);

  /** subtract the input vector to the elements starting from position i.
   *  \param i an unsigned int i
   *  \param v a SiconosVector
   */
  void subBlock(unsigned int i, const SiconosVector& v);

  /** copy the data to double* (dense only)
   * \param data the memory where to copy the data
   * \return the number of element written (size of the vector)
   */
  unsigned copyData(double* data) const;

  /** operator =
   * \param v SiconosVector, the vector to be copied
   * \return  SiconosVector&
   */
  SiconosVector& operator = (const SiconosVector& v);

  /** operator =
   * \param b BlockVector : the vector to be copied
   * \return  SiconosVector&
   */
  SiconosVector& operator = (const BlockVector& b);

  /** operator =
   * \param v a DenseVect : the vector to be copied
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
   * \param v SiconosVector : a vector to add
   * \return  SiconosVector&
   */
  SiconosVector& operator +=(const SiconosVector& v);
  SiconosVector& operator +=(const BlockVector& v);

  /** operator -=
   * \param  v SiconosVector : a vector to subtract
   * \return  SiconosVector&
   */
  SiconosVector& operator -=(const SiconosVector& v);
  SiconosVector& operator -=(const BlockVector& v);

  /** component-wise exponential of a vector
      \param SiconosVector input, such that result (this) = exp(input)
  */
  void exp(SiconosVector& input);

  /** component-wise exponential of a vector, in-place.
      this = exp(this)
  */
  void exp_in_place();

  friend SiconosVector& operator *= (SiconosVector& v, const double& s);

  friend SiconosVector& operator /= (SiconosVector& v, const double& s);

  /** Copy a part of a vector into a sublock of another vector
   *   \param vIn SiconosVector input vector (read only)
   *   \param vOut SiconosVector vector to be overwritten
   *   \param sizeB int size of the block to be copied
   *   \param startIn int starting index of the block to be copied (in input vector)
   *   \param startOut int starting index of the block to be overwritten (in output vector)
   */
  friend void setBlock(const SiconosVector& vIn, SP::SiconosVector vOut,
                       unsigned int sizeB, unsigned int startIn, unsigned int startOut);

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

  /*
  friend SiconosVector abs_wise(const SiconosVector&);

  friend void getMax(const SiconosVector&, double &,unsigned int &);

  friend void  getMin(const SiconosVector&, double &, unsigned int &);
  */

  friend struct IsDense;

  friend struct IsSparse;

  friend struct IsBlock;

  //  temporary workaround, the visitor has to be removed or rework -- xhub
  ACCEPT_NONVIRTUAL_VISITORS();

};



#endif
