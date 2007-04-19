/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

/*! \file SimpleVector.h
 */

#ifndef __SimpleVector__
#define __SimpleVector__

#include "SiconosVector.h"

class SimpleMatrix;

/** Vectors (embedded various types of Boost-Ublas vectors of double).
 *
 * \author SICONOS Development Team - copyright INRIA
 *   \version 2.0.1.
 *   \date (Creation) 07/21/2006
 *
 *
 *
 *
 * SimpleVector is used in the platform to store vectors (mathematical object).
 *
 */
class SimpleVector: public SiconosVector
{
private:

  /** unsigned int num
   * an unsigned int which make a correspondance with Boost Vector: 1 -> DenseVect, 4 -> SparseVect
   * Note: 4 for sparse to keep the same num as for matrices.
   */
  unsigned int num;

  /** Vect vect (See SiconosMatrix.h for more details on Mat type);
   * union of The Boost Matrices : DenseVect, SparseVect are encapsulated.
   */
  Vect vect;

  /** constructor with the type of the Boost vector
   *  \param TYP
   */
  SimpleVector();

public:
  /***************************** CONSTRUCTORS ****************************/

  /** constructor with the type and the dimension of the Boost vector
   *  \param an unsigned int, dimension
   *  \param a TYP
   */
  SimpleVector(unsigned int , TYP = DENSE);

  /** constructor with the type and the dimension of the Boost vector and a default value
   *  \param an unsigned int, dimension
   *  \param double a, so that *this = [a a a ...]
   *  \param a TYP (default = dense)
   */
  SimpleVector(unsigned int , double, TYP = DENSE);

  /** constructor with the type of the boost vector, a std::vector of the values and the dimension of the vector
   *  \param a std::vector<double>
   *  \param an unsigned int
   *  \param a TYP
   */
  SimpleVector(const std::vector<double>, TYP = DENSE);

  /** copy constructor
   *  \param SimpleVector
   */
  SimpleVector(const SimpleVector&);

  /** copy constructor
   *  \param SiconosVector
   */
  SimpleVector(const SiconosVector&);

  /** constructor with a DenseVect vector (see SiconosMatrix.h for details)
   *  \param a DenseVect
   */
  SimpleVector(const DenseVect&);

  /** constructor with a SparseVect vector (see SiconosMatrix.h for details)
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
  ~SimpleVector(void);

  /******************************** METHODS ******************************/

  /** get the attribute num of current vector
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   * \return an unsigned int.
   */
  inline const unsigned int getNum() const
  {
    return num;
  };

  /** get the attribute if it's type is DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect
   */
  const DenseVect getDense(unsigned int = 0) const;

  /** get the attribute if it's type is SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect
   */
  const SparseVect getSparse(unsigned int = 0) const;

  /** get a pointer on DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect*
   */
  DenseVect* getDensePtr(unsigned int = 0) const;

  /** get a pointer on SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect*
   */
  SparseVect* getSparsePtr(unsigned int = 0) const;

  /** return the array of double values of the vector
   *  \exception SiconosVectorException
   *  \param unsigned int: vector position (only for block vector)
   *  \return double* : the pointer on the array
   */
  double* getArray(unsigned int = 0) const;

  /** get block starting at pos and of size block.size()
   *  \param an int, position of the first element of the required block
   *  \param a SiconosVector, in-out parameter.
   */
  void getBlock(unsigned int, SiconosVector&) const;

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
  const double normInf()const;

  /** return the Euclidian norm of the vector
   *  \return a double
   */
  const double norm2() const ;

  /** display data on standard output
   */
  void display(void) const;

  /** return the current object. This function is really usefull only for block vector
   * \return a pointer to a SiconosVector
   */
  inline SiconosVector* getVectorPtr(unsigned int)
  {
    return this;
  };

  /** set all values of the vector component to value.
   * \param a double
   */
  void fill(double);

  /** put data of the vector into a string
   */
  const std::string toString() const;

  //************************** VECTORS HANDLING AND OPERATORS *******************************

  /** return the element vector[i]
   *  \param an unsigned int i
   *  \return a double
   */
  const double getValue(unsigned int);

  /** set the element vector[i]
   *  \param an unsigned int i
   *  \param the value
   */
  void setValue(unsigned int, double);

  /** set the elements starting from position i with input vector
   *  \param an unsigned int i
   *  \param a SiconosVector
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

  /** get the element at position i in the vector
   *  \param an integer i
   *  \return a double
   */
  double& operator()(unsigned int);

  /** get the element at position i in the vector
   *  \param an integer i
   *  \return a double
   */
  const double operator()(unsigned int)const;

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SiconosVector*
   */
  inline SimpleVector* operator [](unsigned int)
  {
    return this;
  };

  /** get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a SiconosVector*
   */
  inline const SimpleVector* operator [](unsigned int) const
  {
    return this;
  };

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

  /** operator /=
   *  \param double, a scalar
   */
  SimpleVector& operator /= (double);

  /** operator /=
   *  \param int, a scalar
   */
  SimpleVector& operator /= (int);

  /** multiply the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   */
  SimpleVector& operator *= (double);

  /** multiply the current vector with an int
   *  \param an int
   *  \exception SiconosVectorException
   */
  SimpleVector& operator *= (int);

  /** Copy x into y */


  /**: A==B when (A-B).normInf()<tolerance
   * \param 2 SiconosVector
   * \return a boolean
   */
  friend bool operator ==(const SiconosVector&, const SiconosVector&);

  /** Addition of two vectors
   *  \param 2 SiconosMatrix
   *  \return a SimpleVector
   *  \exception SiconosVectorException, if the sizes are incompatible
   *  \exception SiconosVectorException, if the two vectors have different types, in this case use function add
   */
  friend SimpleVector operator + (const SimpleVector&, const SimpleVector&);

  /** Subtraction of two vectors
   *  \param 2 SiconosVector
   *  \return a SimpleVector
   *  \exception SiconosVectorException, if the sizes are incompatible
   *  \exception SiconosVectorException, if the two vectors have different types, in this case use function add
   */
  friend SimpleVector operator - (const SimpleVector&, const SimpleVector&);

  /** multiplication of a vector by a double
   *  \param a SiconosVector
   *  \param a double
   *  \return a SimpleVector
   */
  friend SimpleVector operator * (const SimpleVector&, double);

  /** multiplication of a vector by an int
   *  \param a SiconosVector
   *  \param an integer
   *  \return a SimpleVector
   */
  friend SimpleVector operator * (const SimpleVector&, int);

  /** multiplication of a vector by a double
   *  \param a SiconosVector
   *  \param a double
   *  \return a SimpleVector
   */
  friend SimpleVector operator * (double, const SimpleVector&);

  /** multiplication of a vector by an int
   *  \param a SiconosVector
   *  \param an integer
   *  \return a SimpleVector
   */
  friend SimpleVector operator * (int, const SimpleVector&);

  /** division of the vector by a double
   *  \param a SiconosVector
   *  \param a double
   *  \return a SimpleVector
   *  \exception SiconosVectorException, if the double d = 0
   */
  friend SimpleVector operator / (const SimpleVector&, double);

  /** division of the vector by an int
   *  \param a SiconosVector
   *  \param an integer
   *  \return a SimpleVector
   *  \exception SiconosVectorException, if the int d = 0
   */
  friend SimpleVector operator / (const SimpleVector&, int);

  /** compute the product m1 * trans(m2)
   *  \param 2 SiconosVectors
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix outer_prod(const SiconosVector&, const SiconosVector&);

  /** compute dot product m1.m2
   *  \param 2 SiconosVectors
   *  \return a double
   */
  friend const double inner_prod(const SiconosVector&, const SiconosVector&);

  /** computes this = x + y with atlas xpy .
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  void xpy(const SiconosVector &, const SiconosVector &);

  /** computes this = ax + by with atlas axpby .
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  void axpby(double, const SiconosVector&, double, const SiconosVector&);

  /** computes this = a*x with atlas scal.
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  void scal(double, const SiconosVector&);

  /** computes this = a*x + y with atlas axpy.
      \param a SiconosVector (x)
      \param a SiconosVector (y)
  */
  void axpy(double, const SiconosVector&, const SiconosVector&);

  /** computes y = a*x + b*y with atlas axpy.
      \param a SiconosVector (x), IN.
      \param a SiconosVector (y), IN-OUT.
  */
  friend void axpby(double, const SiconosVector&, double, SiconosVector&);

  /** computes y = a*x + y with atlas axpy.
      \param a SiconosVector (x), IN.
      \param a SiconosVector (y), IN-OUT.
  */
  friend void axpy(double, const SiconosVector&, SiconosVector&);

  /** Swap x and y contents, using atlas swap.*/
  friend void swap(SiconosVector&, SiconosVector&);

};
#endif
