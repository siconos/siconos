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

/** \class MySimpleVector
 *   \brief This class is an encapsulation of the Boost class managing vectors of double.
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.3.0.
 *   \date (Creation) 07/21/2006
 *
 *
 *
 *
 * MySimpleVector is used in the platform to store vectors (mathematical object).
 *
 */

#ifndef __MySimpleVector__
#define __MySimpleVector__

#include "MySiconosVector.h"

using namespace boost::numeric::ublas;
class MySimpleMatrix;

class MySimpleVector: public MySiconosVector
{
private:

  /** \var unsigned int num
   * \brief an unsigned int which make a correspondance with Boost Matrices: 1 -> DenseMat, 2 -> TriangMat, 3 -> SymMat
   */
  unsigned int num;

  /** \var MyVect vect (See MySiconosMatrix.h for more details on MyMat type);
   * \brief union of The Boost Matrices : DenseVect, SparseVect are encapsulated.
   */
  MyVect vect;

  /** \fn MySimpleVector (TYP = DENSE)
   *  \brief constructor with the type of the Boost vector
   *  \param TYP
   */
  MySimpleVector(TYP = DENSE);

public:
  /***************************** CONSTRUCTORS ****************************/

  /** \fn MySimpleVector (unsigned int, TYP = DENSE)
   *  \brief constructor with the type and the dimension of the Boost vector
   *  \param an unsigned int, dimension
   *  \param a TYP
   */
  MySimpleVector(unsigned int , TYP = DENSE);

  /** \fn MySimpleVector (const std::vector<double>&, unsigned int, TYP = DENSE)
   *  \brief constructor with the type of the boost vector, a std::vector of the values and the dimension of the vector
   *  \param a std::vector<double>
   *  \param an unsigned int
   *  \param a TYP
   */
  MySimpleVector(const std::vector<double> , unsigned int , TYP = DENSE);

  /** \fn MySimpleVector (const MySimpleVector&)
   *  \brief copy constructor
   *  \param MySimpleVector
   */
  MySimpleVector(const MySimpleVector&);

  /** \fn MySimpleVector (const MySiconosVector&)
   *  \brief copy constructor
   *  \param MySiconosVector
   */
  MySimpleVector(const MySiconosVector&);

  /** \fn MySimpleVector (const DenseVect&)
   *  \brief constructor with a DenseVect vector (see MySiconosMatrix.h for details)
   *  \param a DenseVect
   */
  MySimpleVector(const DenseVect&);

  /** \fn MySimpleVector (const SparseVect&)
   *  \brief constructor with a SparseVect vector (see MySiconosMatrix.h for details)
   *  \param a SparseVect
   */
  MySimpleVector(const SparseVect&);

  /** \fn MySimpleVector (std::string file, bool ascii)
   *  \brief constructor with an input file
   *  \param a std::string which contain the file path
   *  \param a boolean to indicate if the file is in ascii
   */
  MySimpleVector(const std::string&, bool = true);

  /** \fn ~MySimpleVector ()
   *  \brief destructor
   */
  ~MySimpleVector(void);

  /******************************** METHODS ******************************/

  /** \fn void zero();
   *  \brief sets all the values of the vector to 0.0
   */
  void zero(void);

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the vector
   *  \return a double
   */
  const double normInf(void)const;

  /** \fn unsigned unsigned int size() const
   *  \brief get the vector size, ie the total number of (double)
   *  elements in the vector
   *  \return unsigned int
   */
  unsigned int size(void)const;

  /** \fn  void resize (unsigned int nbcol, bool val = true)const
   *  \brief resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   *  \exception SiconosVectorException
   */
  void resize(unsigned int, bool = true);


  /** \fn void display();
   *  \brief display data on standard output
   */
  void display(void)const;

  //************************** VECTORS GETTERS AND SETTERS *******************************

  /** \fn DenseVect getDense()
   *  \brief get the attribute if it's type is DenseVect
   *  \return a DenseVect
   */
  const DenseVect getDense(void)const;

  /** \fn SparseVect getSparse()
   *  \brief get the attribute if it's type is SparseVect
   *  \return a SparseVect
   */
  const SparseVect getSparse(void)const;

  /** \fn DenseVect* getDensePtr()
   *  \brief get a pointer on DenseVect
   *  \return a DenseVect*
   */
  DenseVect* getDensePtr(void) const;

  /** \fn SparseVect* getSparsePtr()
   *  \brief get a pointer on SparseVect
   *  \return a SparseVect*
   */
  SparseVect* getSparsePtr(void)const;

  /** \fn unsigned int getNum() const
   *  \brief get the attribute num of current vector
   * \return an unsigned int.
   */
  unsigned int getNum(void)const;

  /** \fn void setNum(unsigned int n)
   *  \brief set the attribute num of current vector with n
   */
  void setNum(unsigned int);

  //************************** VECTORS HANDLING AND OPERATORS *******************************

  /** \fn double& operator ()(unsigned int i)
   *  \brief get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  double& operator()(unsigned int);

  /** \fn double operator ()(unsigned int i)const
   *  \brief get the element at position i in the vector
   *  \param an integer i
   *  \exception SiconosVectorException
   *  \return a double
   */
  double operator()(unsigned int)const;

  /** \fn operator * (const MySiconosVector &v)
   *  \brief multiply the current vector with the vector v
   *  \param a MySiconosVector
   *  \return a MySimpleVector, the result of the multiplication
   */
  double operator * (const MySiconosVector&);

  /** \fn operator = (const MySiconosVector&)
   *  \param MySiconosVector : the vector to be copied
   */
  const MySimpleVector& operator = (const MySiconosVector&);

  /** \fn operator = (const MySiconosVector&)
   *  \param MySimpleVector : the vector to be copied
   */
  const MySimpleVector& operator = (const MySimpleVector&);

  /** \fn operator += (const MySiconosVector&)
   *  \param MySiconosVector : a vector to add
   */
  const MySimpleVector& operator +=(const MySiconosVector&);

  /** \fn operator -= (const MySiconosVector&)
   *  \param MySiconosVector : a vector to subtract
   */
  const MySimpleVector& operator -=(const MySiconosVector&);

  /** \fn operator /= (double)
   *  \param double, a scalar
   */
  const MySimpleVector& operator /= (double);

  /** \fn operator /= (int)
   *  \param int, a scalar
   */
  const MySimpleVector& operator /= (int);

  /** \fn operator*=(double)
   *  \brief multiply the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   */
  const MySimpleVector& operator *= (double);

  /** \fn operator*=(int)
   *  \brief multiply the current vector with an int
   *  \param an int
   *  \exception SiconosVectorException
   */
  const MySimpleVector& operator *= (int);

  /** \fn operator ==
   * \brief: A==B when (A-B).normInf()<tolerance
   * \param 2 MySiconosVector
   * \return a boolean
   */
  friend  bool            operator ==(const MySiconosVector&, const MySiconosVector&);

  /** \fn operator + (const MySiconosVector& m1, const MySiconosVector& m2);
   *  \brief Addition of two vectors
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleVector
   *  \exception SiconosVectorException, if the sizes are incompatible
   *  \exception SiconosVectorException, if the two vectors have different types, in this case use function add
   */
  friend MySimpleVector operator + (const MySimpleVector&, const MySimpleVector&);

  /** \fn operator - (const MySiconosVector& m1, const MySiconosVector& m2);
   *  \brief Subtraction of two vectors
   *  \param 2 MySiconosVector
   *  \return a MySimpleVector
   *  \exception SiconosVectorException, if the sizes are incompatible
   *  \exception SiconosVectorException, if the two vectors have different types, in this case use function add
   */
  friend MySimpleVector operator - (const MySimpleVector&, const MySimpleVector&);

  /** \fn operator * (const MySiconosVector& m1, double d);
   *  \brief multiplication of a vector by a double
   *  \param a MySiconosVector
   *  \param a double
   *  \return a MySimpleVector
   */
  friend MySimpleVector operator * (const MySimpleVector&, double);

  /** \fn operator * (const MySiconosVector& m1, int d);
   *  \brief multiplication of a vector by an int
   *  \param a MySiconosVector
   *  \param an integer
   *  \return a MySimpleVector
   */
  friend MySimpleVector operator * (const MySimpleVector&, int);

  /** \fn operator * (double d, const MySiconosVector& m1);
   *  \brief multiplication of a vector by a double
   *  \param a MySiconosVector
   *  \param a double
   *  \return a MySimpleVector
   */
  friend MySimpleVector operator * (double, const MySimpleVector&);

  /** \fn operator * (int d, const MySiconosVector& m1);
   *  \brief multiplication of a vector by an int
   *  \param a MySiconosVector
   *  \param an integer
   *  \return a MySimpleVector
   */
  friend MySimpleVector operator * (int, const MySimpleVector&);

  /** \fn operator / (const MySiconosVector& m1, double d);
   *  \brief division of the vector by a double
   *  \param a MySiconosVector
   *  \param a double
   *  \return a MySimpleVector
   *  \exception SiconosVectorException, if the double d = 0
   */
  friend MySimpleVector operator / (const MySimpleVector&, double);

  /** \fn operator / (const MySiconosVector& m1, int d);
   *  \brief division of the vector by an int
   *  \param a MySiconosVector
   *  \param an integer
   *  \return a MySimpleVector
   *  \exception SiconosVectorException, if the int d = 0
   */
  friend MySimpleVector operator / (const MySimpleVector&, int);

  /** \fn add (const MySiconosVector& m1, const MySiconosVector& m2);
   *  \brief Addition of two vectors
   *  \param 2 MySiconosVector
   *  \return a MySimpleVector
   *  \exception SiconosVectorException, if the sizes are incompatible
   */
  friend MySimpleVector add(const MySiconosVector&, const MySiconosVector&);

  /** \fn operator sub (const MySiconosVector& m1, const MySiconosVector& m2);
   *  \brief Subtraction of two vectors
   *  \param 2 MySiconosVector
   *  \return a MySimpleVector
   *  \exception SiconosVectorException, if the sizes are incompatible
   */
  friend MySimpleVector sub(const MySiconosVector&, const MySiconosVector&);

  /** \fn inner_prod (const MySiconosVector& m1, const MySiconosVector& m2);
   *  \brief compute the product trans(m1) * m2
   *  \param 2 MySiconosMatrix
   *  \return a double
   */
  friend double inner_prod(const MySiconosVector&, const MySiconosVector&);

  /** \fn outer_prod (const MySiconosVector& m1, const MySiconosVector& m2);
   *  \brief compute the product m1 * trans(m2)
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleVector
   */
  friend MySimpleMatrix outer_prod(const MySiconosVector&, const MySiconosVector&);
};
#endif
