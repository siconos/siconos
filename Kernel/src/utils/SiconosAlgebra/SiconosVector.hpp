/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
 * You can find an overview on how to build and use vectors and matrices in \ref GS_SicAlgebra .
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
   *  \param an unsigned int, dimension
   *  \param an_UBLAS_TYPE
   */
  SiconosVector(unsigned int , Siconos::UBLAS_TYPE = Siconos::DENSE);

  /** constructor with the dimension of the Boost vector, a default value and the type.
   *  \param an unsigned int, dimension
   *  \param double a, so that *this = [a a a ...]
   *  \param an Siconos::UBLAS_TYPE (default = dense)
   */
  SiconosVector(unsigned int , double, Siconos::UBLAS_TYPE = Siconos::DENSE);

  /** constructor with  a std::vector of the values and the type of the boost vector.
   *  \param a std::vector<double>
   *  \param an Siconos::UBLAS_TYPE
   */
  SiconosVector(const std::vector<double>&, Siconos::UBLAS_TYPE = Siconos::DENSE);

  /** copy constructor
   *  \param SiconosVector
   */
  SiconosVector(const SiconosVector&);

  SiconosVector(const BlockVector & vIn);
  /** constructor with a DenseVect vector (see SiconosAlgebra.hpp for details)
   *  \param a DenseVect
   */
  SiconosVector(const DenseVect&);

  /** constructor with a SparseVect vector (see SiconosAlgebra.hpp for details)
   *  \param a SparseVect
   */
  SiconosVector(const SparseVect&);

  /** constructor with an input file
   *  \param a std::string which contain the file path
   *  \param a boolean to indicate if the file is in ascii
   */
  SiconosVector(const std::string&, bool = true);

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
  unsigned int getNum() const
  {
    if (_dense) return 1;
    else return 4;
  }

  /** get the ublas embedded vector if it's type is Dense
   *  \param unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a DenseVect
   */
  const DenseVect getDense(unsigned int = 0) const;

  /** get the ublas embedded vector if it's type is Sparse
   *  \param unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a SparseVect
   */
  const SparseVect getSparse(unsigned int = 0) const;

  /** get a pointer to the ublas embedded vector if it's type is Dense
   *  \param unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a DenseVect*
   */
  inline DenseVect* dense(unsigned int = 0) const
  {
    return vect.Dense;
  };

  /** get a pointer to the ublas embedded vector if it's type is Sparse
   *  \param unsigned int: position of the required vector (useless for SiconosVector, default = 0)
   *  \return a SparseVect*
   */
  SparseVect* sparse(unsigned int = 0) const;

  /** return the array of double values of the vector
   *  \param unsigned int: vector position (only for block vector)
   *  \return double* : the pointer on the array
   */
  double* getArray(unsigned int = 0) const;

  /** get block starting at "pos" (first argument) and of size block.size() and write it in v (second arg)
   *  \param pos an int, position of the first element of the required block
   *  \param v a SiconosVector *, in-out parameter.
   */
  //  void getBlock(unsigned int, SP::SiconosVector) const;

  /** sets all the values of the vector to 0.0
   */
  void zero();

  /** resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   *  \exception SiconosVectorException
   */
  void resize(unsigned int, bool = true);

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
  double sum() const;

  /** display data on standard output
   */
  void display(void) const;

  /** return the current object. This function is really usefull only for block vector
   * \return a pointer to a SiconosVector
   */
  inline SP::SiconosVector vector(unsigned int)
  {
    return shared_from_this();
  };

  /** return the current object. This function is really usefull only for block vector
   * \return a pointer to a SiconosVector
   */
  inline SPC::SiconosVector vector(unsigned int) const
  {
    return shared_from_this();
  };

  /** set SiconosVector number i (copy) with v (second arg) - Useful only for BlockVector (else equivalent to a single copy)
   * \param unsigned int: block number (0 for SiconosVector)
   * \param a SiconosVector
   */
  void setVector(unsigned int, const SiconosVector&);

  /** set SiconosVector number i (pointer link) with v (second arg) - Useful only for BlockVector
   * \param unsigned int: block number (0 for SiconosVector)
   * \param a pointer to a SiconosVector
   */
  inline void setVectorPtr(unsigned int, SP::SiconosVector)
  {
    SiconosVectorException::selfThrow("SiconosVector::setVectorPtr(num,v), not allowed for SiconosVector.");
  };

  /** set all values of the vector component to input value.
   * \param a double
   */
  void fill(double);

  /** put data of the vector into a string
   */
  const std::string toString() const;

  //************************** VECTORS HANDLING AND OPERATORS *******************************

  /** get the element at position i in the vector
   *  \param an unsigned int i
   *  \return a double
   */
  double getValue(unsigned int) const ;

  /** set the element at position i in the vector.
   *  \param an unsigned int i
   *  \param the value
   */
  void setValue(unsigned int, double);

  /** get the element at position i in the vector
   *  \param an integer i
   *  \return a double
   */
  double& operator()(unsigned int);

  /** get the element at position i in the vector
   *  \param an integer i
   *  \return a double
   */
  double operator()(unsigned int)const;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SP::SiconosVector
   */
  inline SP::SiconosVector operator [](unsigned int)
  {
    return shared_from_this();
  };

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SP::SiconosVector
   */
  inline SPC::SiconosVector operator [](unsigned int) const
  {
    return shared_from_this();
  };

  /** set the elements starting from position i with input vector
   *  \param an unsigned int i
   *  \param a ref to a SiconosVector
   */
  void setBlock(unsigned int, const SiconosVector&);

  /** copy a part of the vector into another (usefull for converting a SiconosVector
   * into a BlockVector
   */
  void toBlock(SiconosVector& vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut) const;

  /** add the input vector to the elements starting from position i.
   *  \param an unsigned int i
   *  \param a SiconosVector
   */
  void addBlock(unsigned int, const SiconosVector&);

  /** subtract the input vector to the elements starting from position i.
   *  \param an unsigned int i
   *  \param a SiconosVector
   */
  void subBlock(unsigned int, const SiconosVector&);

  /** operator =
   *  \param SiconosVector : the vector to be copied
   */
  SiconosVector& operator = (const SiconosVector&);

  /** operator =
   *  \param BlockVector : the vector to be copied
   */
  SiconosVector& operator = (const BlockVector&);

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  SiconosVector& operator = (const DenseVect&);

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  SiconosVector& operator = (const SparseVect&);

  /** operator +=
   *  \param SiconosVector : a vector to add
   */
  SiconosVector& operator +=(const SiconosVector&);
  SiconosVector& operator +=(const BlockVector&);

  /** operator -=
   *  \param SiconosVector : a vector to subtract
   */
  SiconosVector& operator -=(const SiconosVector&);
  SiconosVector& operator -=(const BlockVector&);


  friend SiconosVector& operator *= (SiconosVector& v, const double& s);

  friend SiconosVector& operator /= (SiconosVector& v, const double& s);

  /** Copy a part of a vector into a sublock of another vector
      \param SiconosVector input vector (read only)
      \param SiconosVector vector to be overwritten
      \param int size of the block to be copied
      \param int starting index of the block to be copied (in input vector)
      \param int starting index of the block to be overwritten (in output vector)
  **/
  friend void setBlock(const SiconosVector&, SP::SiconosVector, unsigned int, unsigned int, unsigned int);

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

  friend struct VectorNum;

  friend struct IsDense;

  friend struct IsSparse;

  friend struct IsBlock;

  //  temporary workaround, the visitor has to be removed or rework -- xhub
  ACCEPT_NONVIRTUAL_VISITORS();

};

struct VectorNum : public Question<unsigned int>
{

  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    if (v._dense) answer = 1;
    else answer = 4;
  }

  void visit(std11::shared_ptr<SiconosVector> v)
  {
    if (v->_dense) answer = 1;
    else answer = 4;
  }

  void visit(const BlockVector& v)
  {
    answer = 0;
  }

};

struct IsDense : public Question<bool>
{
  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    answer = v._dense;
  }

  void visit(const BlockVector& v)
  {
    answer = false;
  }
};

struct IsSparse : public Question<bool>
{

  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    answer = !v._dense;
  }

  void visit(const BlockVector& v)
  {
    answer = false;
  }
};

struct IsBlock : public Question<bool>
{
  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    answer = false;
  }

  void visit(const BlockVector& v)
  {
    answer = true;
  }
};


#endif
