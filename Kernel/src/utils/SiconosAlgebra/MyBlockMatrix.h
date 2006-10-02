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
/** \class MyBlockMatrix
 *   \brief list of MySimpleMatrix * to work with "block-matrices"
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.3.0.
 *   \date (Creation) 21/07/2006
 *
 *  A block matrix is a map of pointers to MySimpleMatrix, with two indices as a key.
 *  This indices give row and column position of the block in the global matrix.
 *
 * If block 1 is n1xm1, block2 n2xm2, block3 n3xm3 ..., then:
 *  tabRow = [ n1 n1+n2 n1+n2+n3 ...]
 *  tabCol = [ m1 m1+m2 m1+m2+m3 ...]
 *
 */
#ifndef __MyBlockMatrix__
#define __MyBlockMatrix__

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <vector>
#include <deque>
#include <cassert>

#include "MySimpleMatrix.h"
/** \constant STDMAP
 *  \brief if STDMAP then we use an iterator of std::map, else we use an iterator1 and an iterator2 of map_std of Boost
 */
#define STDMAP 0

using namespace boost::numeric::ublas;

class MyBlockMatrix : public MySiconosMatrix
{
private:

  /** \var map
   * \brief a map of the blocks that compose the matrix
   */
  mapped map;

  /** \var isBlockAllocatedIn
   * \brief a list of bool, to check inside-class allocation for SimpleMatrix blocks.
   */
  std::deque<bool> isBlockAllocatedIn;

  /** \var tabRow
   * \brief list of blocks dimension (number of rows)
   */
  std::vector<unsigned int> tabRow;

  /** \var tabCol
   * \brief list of blocks dimension (number of columns)
   */
  std::vector<unsigned int> tabCol;

  /** \fn void makeTab(unsigned int row, unsigned int col)
   *  \brief build tabRow and tabCol
   *  \param unsigned int i: number of block-rows in this
   *  \param unsigned int i: number of block-columns in this
   */
  void makeTab(unsigned int, unsigned int);

  /** \fn addInTab(unsigned int i, unsigned int j)
   *  \brief change tabRow and tabCol, according to a new input matrix of dim i x j
   *  \param unsigned int i: number of rows in the new matrix
   *  \param unsigned int i: number of columns in the new matrix
   */
  void addInTab(unsigned int, unsigned int);

  /** \fn MyBlockMatrix ()
   *  \brief default constructor
   */
  MyBlockMatrix(void);

public:

  /** \fn MyBlockMatrix (const MySiconosMatrix&  m)
   *  \brief copy constructor
   *  \param MySiconosMatrix
   */
  MyBlockMatrix(const MySiconosMatrix&);

  /** \fn MyBlockMatrix (const MyBlockMatrix&  m)
   *  \brief copy constructor
   *  \param MyBlockMatrix
   */
  MyBlockMatrix(const MyBlockMatrix&);

  /** \fn MyBlockMatrix (const mapped&  m)
   *  \brief constructor with a map
   *  \param MyBlockMatrix
   */
  MyBlockMatrix(mapped&);

  /** \fn MyBlockMatrix (vector<MySiconosMatrix*>v, unsigned int row, unsigned int col)
   *  \brief constructor with a list of pointer to MySiconosMatrix (!links with pointer, no copy!)
   *  \param vector<MySiconosMatrix*>
   *  \param unsigned int: number of blocks in a row
   *  \param unsigned int: number of col in a row
   */
  MyBlockMatrix(const std::vector<MySiconosMatrix* >&, unsigned int, unsigned int);

  /** \fn ~MyBlockMatrix ()
   *  \brief destructor
   */
  ~MyBlockMatrix(void);

  /** \fn  unsigned int size1 (void)const
   *  \brief get the number of rows of the Block matrix
   *  \exception SiconosMatrixException
   *  \return the number of rows of the Block matrix
   */
  unsigned int  size1(void)const;

  /** \fn  unsigned int size2 (void)const
   *  \brief get the number of columns of the Block matrix
   *  \exception SiconosMatrixException
   *  \return the number of columns of the Block matrix
   */
  unsigned int  size2(void)const;

  /** \fn  void resize (unsigned int nbrow, unsigned int nbcol, unsigned int lower=0, unsigned int upper=0, bool val = true)const
   *  \brief resize the matrix with nbrow rows and nbcol columns, lower and upper are useful only for SparseMat.The existing elements of the Block matrix are preseved when specified.
   *  \exception SiconosMatrixException
   */
  void resize(unsigned int, unsigned int, unsigned int lower = 0, unsigned int upper = 0, bool = true);

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the Block matrix
   *  \return a double
   */
  const double normInf(void)const;

  /** \fn void zero();
   *  \brief sets all the values of the matrix to 0.0
   */
  void zero(void);

  /** \fn void eye();
   *  \brief set an identity matrix
   */
  void eye(void);

  /** \fn void display();
   *  \brief display data on standard output
   */
  void display(void)const;

  // ******************************* GETTERS/SETTERS **************************

  /** \fn void getBlock(unsigned int row, unsigned int col, MySiconosMatrix&)
   *  \brief get a block of a matrix, which position is defined by row and col
   *  \param 2 unsigned int for indexes and a MySiconosMatrix (in-out paramater)
   */
  void  getBlock(unsigned int, unsigned int, MySiconosMatrix&)const;

  /** \fn const std::deque<bool> getBlockAllocated()const
   *   \brief get std::deque of bool
   *   \useless for MySimpleMatrix
   *   \return a std::deque<bool>
   */
  const std::deque<bool> getBlockAllocated(void)const;

  /** \fn unsigned int getNum() const
   *  \brief get the attribute num of current matrix, useless for Block Matrix
   * \return an unsigned int.
   */
  unsigned int   getNum(void)const;

  /** \fn void getRow(unsigned int index, MySimpleVector& vOut) const
   *  \brief get row index of current matrix and save it unsigned into vOut
   *  \param unsigned int: index of required line
   *  \param ref to MySimpleVector: in-out parameter
   */
  void  getRow(unsigned int, MySimpleVector&) const;

  /** \fn void getCol(unsigned int index, MySimpleVector& vOut) const
   *  \brief get column index of current matrix and save it into vOut
   *  \param unsigned int: index of required column
   *  \param ref to MySimpleVector: in-out parameter
   */
  void  getCol(unsigned int, MySimpleVector&) const;

  /** \fn DenseMat* getDense(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a DenseMat
   */
  const DenseMat getDense(unsigned int = 0, unsigned int = 0) const;

  /** \fn TriangMat getTriang(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a TriangMat
   */
  const TriangMat getTriang(unsigned int = 0, unsigned int = 0) const;

  /** \fn SymMat getSym(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SymMat
   */
  const SymMat     getSym(unsigned int = 0, unsigned int = 0)const;

  /** \fn BandedMat getBanded(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a BandedMat
   */
  const BandedMat getBanded(unsigned int = 0, unsigned int = 0)const;

  /** \fn SparseMat getSparse(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SparseMat
   */
  const SparseMat getSparse(unsigned int = 0, unsigned int = 0)const;

  /** \fn DenseMat* getDensePtr(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get a pointer on DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a DenseMat*
   */
  const DenseMat*  getDensePtr(unsigned int = 0, unsigned int = 0)const;

  /** \fn TriangMat* getTriangPtr(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get a pointer on TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a TriangMat*
   */
  const TriangMat* getTriangPtr(unsigned int = 0, unsigned int = 0)const;

  /** \fn SymMat* getSymPtr(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get a pointer on SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SymMat*
   */
  const SymMat*    getSymPtr(unsigned int = 0, unsigned int = 0)const;

  /** \fn BandedMat* getBandedPtr(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get a pointer on BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a BandedMat*
   */
  const BandedMat* getBandedPtr(unsigned int = 0, unsigned int = 0)const;

  /** \fn SparseMat* getSparsePtr(unsigned int row = 0, unsigned int col = 0) const
   *  \brief get a pointer on SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SparseMat*
   */
  const SparseMat* getSparsePtr(unsigned int = 0, unsigned int = 0)const;

  /** \fn mapped getMap() const
   *  \brief get mapped matrix
   *  \useless for MySimpleMatrix
   *  \return a mapped
   */
  virtual const mapped getMap(void)const;

  /** \fn void setNum(unsigned int n)
   *  \brief set the attribute num of current matrix with n
   */
  void  setNum(unsigned int);

  /** \fn void setRow(unsigned int row, const MySimpleVector &v)
   *  \brief set line row of the current matrix with vector v
   *  \param an unsigned int and a MySimpleVector
   */
  void  setRow(unsigned int, const MySimpleVector&);

  /** \fn void setCol(unsigned int col, const MySimpleVector &v)
   *  \brief set column col of the current matrix with vector v
   *  \param an unsigned int and a MySimpleVector
   */
  void  setCol(unsigned int, const MySimpleVector&);

  //******************* MATRICES HANDLING AND OPERATORS ***********************


  /** \fn void blockMatrixCopy( MySiconosMatrix &blockMat, unsigned int, unsigned int)
   *  \brief copy the blockmatrix "blockMat" in the matrix "mat" at the position (xPos, yPos)
   *                    [0, 0, 0, 0].blockMatrixCopy([1], 0, 2) => mat = [0 0 1 0]
   *  \param MySiconosMatrix& : the block matrix to copy in the current matrix
   *  \param unsigned int : the line position to start the copy of the blockmatrix
   *  \param unsigned int : the column position to start the copy of the blockmatrix
   */
  void blockMatrixCopy(const MySiconosMatrix&, unsigned int, unsigned int);

  /** \fn double& operator() (unsigned int row, unsigned int col)
   *  \brief get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double& operator()(unsigned int, unsigned int);

  /** \fn double operator() (unsigned int row, unsigned int col)const
   *  \brief get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double operator()(unsigned int, unsigned int)const;

  /** \fn assignment operator
   *  \param MySiconosMatrix : the matrix to be copied
   */
  const MyBlockMatrix& operator = (const MySiconosMatrix&);

  /** \fn assignment operator
   *  \param MyBlockMatrix : the matrix to be copied
   */
  const MyBlockMatrix& operator = (const MyBlockMatrix&);

  /** \fn operator +=
   *  \param MySiconosMatrix : a matrix to add
   */
  const MyBlockMatrix& operator +=(const MySiconosMatrix&);

  /** \fn operator -=
   *  \param MySiconosMatrix : a matrix to subtract
   */
  const MyBlockMatrix& operator -=(const MySiconosMatrix&);

  /** \fn operator /=
   *  \param double, a scalar
   */
  const MyBlockMatrix& operator /=(double);

  /** \fn operator /=
   *  \param int, a scalar
   */
  const MyBlockMatrix& operator /=(int);

  /** \fn operator *=
   *  \param double, a scalar
   */
  const MyBlockMatrix& operator *=(double);

  /** \fn operator *=
   *  \param int, a scalar
   */
  const MyBlockMatrix& operator *=(int);

};

#endif
