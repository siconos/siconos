/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

/*! \file SimpleVector.hpp
 */

#ifndef __SimpleVector__
#define __SimpleVector__

#include "SiconosVector.hpp"
#include "SimpleVectorFriends.hpp"

/** Vectors of double. (Interface to various types of Boost-Ublas vectors).
 *
 * \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 07/21/2006
 *
 * Used to handle vectors of double.
 *
 * Two possible types: DENSE (default) and SPARSE.
 *
 * You can find an overview on how to build and use vectors and matrices in \ref GS_SicAlgebra .
 *
 */
class SimpleVector: public SiconosVector , public boost::enable_shared_from_this<SimpleVector>
{
protected:

  bool _dense;

  /**
   * Union of pointers to the ublas vector type (dense or sparse)
   */
  VECTOR_UBLAS_TYPE vect;

public:
  /** Default constructor
   */
  SimpleVector();

  /***************************** CONSTRUCTORS ****************************/

  /** constructor with the type and the dimension of the Boost vector
   *  \param an unsigned int, dimension
   *  \param an_UBLAS_TYPE
   */
  SimpleVector(unsigned int , UBLAS_TYPE = DENSE);

  /** constructor with the dimension of the Boost vector, a default value and the type.
   *  \param an unsigned int, dimension
   *  \param double a, so that *this = [a a a ...]
   *  \param an UBLAS_TYPE (default = dense)
   */
  SimpleVector(unsigned int , double, UBLAS_TYPE = DENSE);

  /** constructor with  a std::vector of the values and the type of the boost vector.
   *  \param a std::vector<double>
   *  \param an UBLAS_TYPE
   */
  SimpleVector(const std::vector<double>&, UBLAS_TYPE = DENSE);

  /** copy constructor
   *  \param SimpleVector
   */
  SimpleVector(const SimpleVector&);

  /** copy constructor
   *  \param SiconosVector
   */
  SimpleVector(const SiconosVector&);

  /** constructor with a DenseVect vector (see SiconosAlgebra.hpp for details)
   *  \param a DenseVect
   */
  SimpleVector(const DenseVect&);

  /** constructor with a SparseVect vector (see SiconosAlgebra.hpp for details)
   *  \param a SparseVect
   */
  SimpleVector(const SparseVect&);

  /** constructor with an input file
   *  \param a std::string which contain the file path
   *  \param a boolean to indicate if the file is in ascii
   */
  SimpleVector(const std::string&, bool = true);

  /** destructor
   */
  ~SimpleVector();

  /** get the vector size, ie the total number of (double) elements in the vector
   *  \return unsigned int
   */
  unsigned int size() const
  {
    if (!_dense)
    {
      return (vect.Sparse->size());
    }
    else
    {
      return (vect.Dense->size());
    }
  }

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
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect
   */
  const DenseVect getDense(unsigned int = 0) const;

  /** get the ublas embedded vector if it's type is Sparse
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect
   */
  const SparseVect getSparse(unsigned int = 0) const;

  /** get a pointer to the ublas embedded vector if it's type is Dense
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect*
   */
  inline DenseVect* dense(unsigned int = 0) const
  {
    return vect.Dense;
  };

  /** get a pointer to the ublas embedded vector if it's type is Sparse
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
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
   * \param unsigned int: block number (0 for SimpleVector)
   * \param a SiconosVector
   */
  void setVector(unsigned int, const SiconosVector&);

  /** set SiconosVector number i (pointer link) with v (second arg) - Useful only for BlockVector
   * \param unsigned int: block number (0 for SimpleVector)
   * \param a pointer to a SiconosVector
   */
  inline void setVectorPtr(unsigned int, SP::SiconosVector)
  {
    SiconosVectorException::selfThrow("SimpleVector::setVectorPtr(num,v), not allowed for SimpleVector.");
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
  SimpleVector& operator = (const SiconosVector&);

  /** operator =
   *  \param SimpleVector : the vector to be copied
   */
  SimpleVector& operator = (const SimpleVector&);

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  SimpleVector& operator = (const DenseVect&);

  /** operator =
   *  \param a DenseVect : the vector to be copied
   */
  SimpleVector& operator = (const SparseVect&);

  /** operator +=
   *  \param SiconosVector : a vector to add
   */
  SimpleVector& operator +=(const SiconosVector&);

  /** operator -=
   *  \param SiconosVector : a vector to subtract
   */
  SimpleVector& operator -=(const SiconosVector&);


  SiconosVector& operator *= (double s)
  {
    if (_dense)
      //atlas::scal((double)m,*vect.Dense);
      *dense() *= s;
    else
      *sparse() *= s;
    return *this;
  }


  SiconosVector& operator /= (double s)
  {
    if (_dense)
      //atlas::scal((double)m,*vect.Dense);
      *dense() /= s;
    else
      *sparse() /= s;
    return *this;
  }

  friend void setBlock(const SiconosVector&, SP::SiconosVector, unsigned int, unsigned int, unsigned int);

  friend bool operator ==(const SiconosVector&, const SiconosVector&);

  friend SimpleVector operator * (double, const SiconosVector&);

  friend SimpleVector operator * (const SiconosVector&, double);

  friend SimpleVector operator / (const SimpleVector&, double);

  friend SimpleVector operator + (const SiconosVector&, const SiconosVector&);

  friend void add(const SiconosVector&, const SiconosVector&, SiconosVector&);

  friend SimpleVector operator - (const SiconosVector&, const SiconosVector&);

  friend void sub(const SiconosVector&, const SiconosVector&, SiconosVector&);

  friend void axpby(double, const SiconosVector&, double, SiconosVector&);

  friend void axpy(double, const SiconosVector&, SiconosVector&);

  friend double inner_prod(const SiconosVector&, const SiconosVector&);

  friend SimpleMatrix outer_prod(const SiconosVector&, const SiconosVector&);

  friend void scal(double, const SiconosVector&, SiconosVector&, bool);

  friend void subscal(double, const SiconosVector&, SiconosVector&, const Index&, bool);

  friend void cross_product(const SiconosVector&, const SiconosVector&, SiconosVector&);

  friend class VectorNum;

  friend class IsDense;

  friend class IsSparse;

  friend class IsBlock;

  ACCEPT_VISITORS();

};


struct VectorNum : public SiconosVisitor
{
private :
  unsigned int _answer;

public:

  void visit(const SimpleVector& v)
  {
    if (v._dense) _answer = 1;
    else _answer = 4;
  }

  void visit(boost::shared_ptr<SimpleVector> v)
  {
    if (v->_dense) _answer = 1;
    else _answer = 4;
  }

  void visit(const BlockVector& v)
  {
    _answer = 0;
  }

  void visit(boost::shared_ptr<BlockVector> v)
  {
    _answer = 0;
  }

};

struct IsDense : public SiconosVisitor
{

  bool _answer;

  void visit(const SimpleVector& v)
  {
    _answer = v._dense;
  }

  void visit(boost::shared_ptr<SimpleVector> v)
  {
    _answer = v->_dense;
  }

  void visit(const BlockVector& v)
  {
    _answer = false;
  }

  void visit(boost::shared_ptr<BlockVector> v)
  {
    _answer = false;
  }

};

struct IsSparse : public SiconosVisitor
{

  bool _answer;

  void visit(const SimpleVector& v)
  {
    _answer = !v._dense;
  }

  void visit(boost::shared_ptr<SimpleVector> v)
  {
    _answer = !v->_dense;
  }

  void visit(const BlockVector& v)
  {
    _answer = false;
  }

  void visit(boost::shared_ptr<BlockVector> v)
  {
    _answer = false;
  }

};

struct IsBlock : public SiconosVisitor
{

  bool _answer;

  void visit(const SimpleVector& v)
  {
    _answer = false;
  }

  void visit(boost::shared_ptr<SimpleVector> v)
  {
    _answer = false;
  }

  void visit(const BlockVector& v)
  {
    _answer = true;
  }

  void visit(boost::shared_ptr<BlockVector> v)
  {
    _answer = true;
  }
};

DEFINE_SPTR(SimpleVector);

#endif
