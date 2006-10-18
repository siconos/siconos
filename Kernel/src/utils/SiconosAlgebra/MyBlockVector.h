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

/** \class MyBlockVector
 *  \brief This class describes MyBlock Vectors, containers of several Siconos Vectors (that should be SimpleVectors)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *
 */

#ifndef MYBLOCKVECTOR_H
#define MYBLOCKVECTOR_H

#include "MySimpleVector.h"
class MyBlockVector : public MySiconosVector
{
private:
  // A container of pointers on SiconosVector (that are to be SimpleVector : no MyBlock of MyBlock allowed
  BlocksVect vect;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  std::deque<bool> isBlockAllocatedIn;

  // tabindex[i] = tabindex[i-1] + ni, ni being the size of svref[i].
  Index tabIndex;

public:

  /** \fn MyBlockVector()
   *  \brief default contructor
   */
  MyBlockVector();

  /** \fn MyBlockVector(const std::string&, const bool)
   *  \brief contructor from data by function read call
   *  \param a string
   *  \param a bool
   */
  MyBlockVector(const std::string&, const bool);

  /** \fn MyBlockVector(const MySiconosVector& v)
   *  \brief contructor with a MySiconosVector (copy)
   *  \param MySiconosVector& v
   */
  MyBlockVector(const MySiconosVector&);

  /** \fn MyBlockVector(const MyBlockVector& v)
   *  \brief copy contructor
   *  \param MyBlockVector& v
   */
  MyBlockVector(const MyBlockVector&);

  /** \fn MyBlockVector(MySiconosVector* v1, MySiconosVector* v2)
   *  \brief contructor with a 2 MySiconosVectors
   *  \param MySiconosVector* v1
   *  \param MySiconosVector* v2
   */
  MyBlockVector(MySiconosVector*, MySiconosVector*);

  /** \fn MyBlockVector(unsigned int i, unsigned int j)
   *  \brief constructor with the number of MyBlocks and their dimension (ie all MyBlocks have the same dim)
   *  \param unsigned int : number of MyBlocks
   *  \param unsigned int : dim of each MyBlock
   */
  MyBlockVector(unsigned int, unsigned int);

  /** \fn ~SiconosVector ()
   *  \brief destructor
   */
  ~MyBlockVector();

  /** \fn BlocksVect getAllVect() const
   *  \brief get vect, ie all the vectors of the object
   * \return a BlocksVect
   */
  inline BlocksVect getAllVect() const
  {
    return vect;
  }

  /** \fn unsigned int getNum() const
   *  \brief get the attribute num of current vector
   * \return an unsigned int.
   */
  unsigned int getNum() const;

  /** \fn DenseVect getDense(unsigned int = 0)
   *  \brief get the attribute if it's type is DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect
   */
  const DenseVect getDense(unsigned int = 0) const;

  /** \fn SparseVect getSparse(unsigned int = 0)
   *  \brief get the attribute if it's type is SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect
   */
  const SparseVect getSparse(unsigned int = 0) const;

  /** \fn DenseVect* getDensePtr(unsigned int = 0)
   *  \brief get a pointer on DenseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a DenseVect*
   */
  DenseVect* getDensePtr(unsigned int = 0) const;

  /** \fn SparseVect* getSparsePtr(unsigned int = 0)
   *  \brief get a pointer on SparseVect
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   *  \return a SparseVect*
   */
  SparseVect* getSparsePtr(unsigned int = 0) const;

  /** \fn bool check() const
   *  \brief return false if one of the block is a null pointer
   * \return a bool
   */
  bool check() const;

  /** \fn void zero();
   *  \brief sets all the values of the vector to 0.0
   *  \param unsigned int: position of the required vector (useless for SimpleVector, default = 0)
   */
  void zero();

  /** \fn unsigned int size() const
   *  \brief get the vector size, ie the total number of (double)
   *  elements in the vector
   *  \return unsigned int
   */
  inline unsigned int size() const
  {
    return tabIndex[tabIndex.size() - 1];
  };

  /** \fn unsigned int getNumberOfBlocks() const
   *  \brief get the number of SimpleVector-Blocks
   *  \return unsigned int
   */
  inline unsigned int getNumberOfBlocks() const
  {
    return tabIndex.size();
  };

  /** \fn  void resize (unsigned int nbcol, bool val = true)const
   *  \brief resize the vector with nbcol columns. The existing elements of the matrix are preseved when specified.
   *  \exception SiconosVectorException
   */
  void resize(unsigned int, bool = true);

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the vector
   *  \return a double
   */
  const double normInf() const;

  /** \fn double norm()
   *  \brief return the Euclidian norm of the vector
   *  \return a double
   */
  const double norm() const ;

  /** \fn void display();
   *  \brief display data on standard output
   */
  void display(void) const;

  /** \fn MySimpleVector getVector(unsigned int) const;
   *  \brief return i-eme MySiconosVector of vect
   * \return a MySimpleVector
   */
  MySimpleVector getVector(unsigned int) const;

  /** \fn MySiconosVector* getVectorPtr(unsigned int) const;
   *  \brief return i-eme MySiconosVector of vect
   * \return a pointer to a MySiconosVector
   */
  MySiconosVector* getVectorPtr(unsigned int);

  /** \fn void fill(double value);
   *  \brief set all values of the vector component to value.
   * \param a double
   */
  void fill(double);

  /** \fn Index getTabIndex() const
   *  \brief get the index tab
   * \return a standard vector of int
   */
  inline Index getTabIndex() const
  {
    return tabIndex;
  }

  /** \fn double& operator() (unsigned int index)
   *  \brief get the element at position i (warning: absolute position.)
   *  \param an unsigned integer i
   *  \return a reference to a double
   */
  double& operator()(unsigned int) ;

  /** \fn double operator() (unsigned int index) const
   *  \brief get the element at position i (warning: absolute position.)
   *  \param an unsigned integer i
   *  \return a double
   */
  double operator()(unsigned int) const;

  /** \fn MySiconosVector* operator[] (unsigned int i)
   *  \brief get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a MySiconosVector*
   */
  MySiconosVector* operator [](unsigned int) ;

  /** \fn MySiconosVector* operator[] (unsigned int i)
   *  \brief get the vector at position i(ie this for Simple and block i for BlockVector)
   *  \param an unsigned integer i
   *  \return a MySiconosVector*
   */
  const MySiconosVector* operator [](unsigned int) const;

  /** \fn operator = (const MySiconosVector&)
   *  \param MySiconosVector : the vector to be copied
   */
  MyBlockVector& operator =(const MySiconosVector&);

  /** \fn operator += (const MySiconosVector&)
   *  \param MySiconosVector : a vector to add
   */
  MyBlockVector& operator +=(const MySiconosVector&);

  /** \fn operator -= (const MySiconosVector&)
   *  \param MySiconosVector : a vector to subtract
   */
  MyBlockVector& operator -=(const MySiconosVector&);

  /** \fn operator /= (double)
   *  \param double, a scalar
   */
  MyBlockVector& operator /=(double);

  /** \fn operator /= (int)
   *  \param int, a scalar
   */
  MyBlockVector& operator /=(int);

  /** \fn operator*=(double)
   *  \brief multiply the current vector with a double
   *  \param a double
   *  \exception SiconosVectorException
   */
  MyBlockVector& operator *=(double);

  /** \fn operator*=(int)
   *  \brief multiply the current vector with an int
   *  \param an int
   *  \exception SiconosVectorException
   */
  MyBlockVector& operator *=(int);

  /** \fn bool add(const MySiconosVector& v)
   *  \brief add a subvector in this vector: allocation and copy
   *  \param MySiconosVector& v : the vector to add
   *  \exception SiconosVectorException
   */
  void add(const  MySiconosVector&) ;

  /** \fn bool addPtr(MySiconosVector* v)
   *  \brief add a pointer to a subvector in this vector: no reallocation nor copy.
   *  \param a pointer to MySiconosVector*
   *  \exception SiconosVectorException
   */
  void addPtr(MySiconosVector*) ;
};

#endif
