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

/** \class MySiconosMatrix : abstract class
 *  \brief this class provides an interface for matrices in Siconos
 *  \author SICONOS Development Team - copyright INRIA
 *  \author NDIAYE Abdel-Aziz
 *  \date (creation) 07/21/2006
 *  Matrices can be either block or Simple.
 *  See Derived classes for details.
 */


#ifndef __MySiconosMatrix__
#define __MySiconosMatrix__

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <cassert>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "SiconosMatrixException.h"
#include "MySimpleVector.h"

using namespace boost::numeric::ublas;
//class MySimpleVector;
class MySimpleMatrix;

/**\brief DenseMat is a typedef of boost::ublas::numeric::matrix<double, row_major, std::vector<double> >
 */
typedef matrix<double, row_major, std::vector<double> > DenseMat;

/**\brief TriangMat is a typedef of boost::ublas::numeric::triangular_matrix<double, upper, row_major, std::vector<double> >
 */
typedef triangular_matrix<double, upper, row_major, std::vector<double> > TriangMat;

/**\brief SymMat is a typedef of boost::ublas::numeric::symmetric_matrix<double, upper, row_major, std::vector<double> >
 */
typedef symmetric_matrix<double, upper, row_major, std::vector<double> > SymMat;

/**\brief BandedMat is a typedef of boost::ublas::numeric::banded_matrix<double, row_major, std::vector<double> >
 */
typedef banded_matrix<double, row_major, std::vector<double> > BandedMat;

/**\brief SparseMat is a typedef of boost::ublas::numeric::mapped_matrix<double>
 */
typedef mapped_matrix<double> SparseMat;

/**\brief MyMat is an union of DenseMat pointer, TriangMat pointer BandedMat, SparseMat and SymMat pointer
 */
union MyMat
{
  DenseMat *Dense;
  TriangMat *Triang;
  SymMat *Sym;
  SparseMat *Sparse;
  BandedMat *Banded;
};

//enum TYP {DENSE=1, TRIANGULAR, SYMMETRIC, SPARSE, BANDED};

class MySiconosMatrix;

/**\brief mapped is a typedef of boost::ublas::numeric::mapped_matrix<MySiconosMatrix* >
 */
typedef mapped_matrix<MySiconosMatrix* > mapped;

class MySiconosMatrix
{
protected:
  /**\var isBlockMatrix
   * bool to check the type of the current matrix; true if block else false. */
  bool isBlockMatrix;

  /**\fn setIsBlockMatrix (bool)
   * \brief set the value of isBlockMatrix.
   * \param a bool to set isBlockMatrix (optional, default = false)
   */
  inline void setIsBlock(bool val)
  {
    isBlockMatrix = val;
  }

public:

  /**\fn bool isBlock ()
   * \brief true if the matrix is block else false.
   * \return a bool.*/
  inline bool isBlock(void)const
  {
    return isBlockMatrix;
  }

  /**\fn ~MySiconosMatrix ()
   * \brief Destructor. */
  virtual ~MySiconosMatrix(void) = 0;

  // ************ GETTERS/SETTERS ***************


  /** \fn unsigned int getNum() const = 0
   *  \brief get the attribute num of current matrix
   * \return an unsigned int.
   */
  virtual unsigned int getNum(void) const = 0;

  /** \fn DenseMat getDense(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a DenseMat
   */
  virtual const DenseMat getDense(unsigned int = 0, unsigned int = 0) const = 0;

  /** \fn TriangMat getTriang(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a TriangMat
   */
  virtual const TriangMat getTriang(unsigned int = 0, unsigned int = 0) const  = 0;

  /** \fn SymMat getSym(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SymMat
   */
  virtual const SymMat getSym(unsigned int = 0, unsigned int = 0) const  = 0;

  /** \fn BandedMat getBanded(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a BandedMat
   */
  virtual const BandedMat getBanded(unsigned int = 0, unsigned int = 0) const  = 0;

  /** \fn SparseMat getSparse(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SparseMat
   */
  virtual const SparseMat getSparse(unsigned int = 0, unsigned int = 0) const = 0;

  /** \fn DenseMat* getDensePtr(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get a pointer on DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a DenseMat*
   */
  virtual const DenseMat* getDensePtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** \fn TriangMat* getTriangPtr(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get a pointer on TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a TriangMat*
   */
  virtual const TriangMat* getTriangPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** \fn SymMat* getSymPtr(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get a pointer on SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SymMat*
   */
  virtual const SymMat* getSymPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** \fn BandedMat* getBandedPtr(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get a pointer on BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a BandedMat*
   */
  virtual const BandedMat* getBandedPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** \fn SparseMat* getSparsePtr(unsigned int row = 0, unsigned int col = 0) const = 0
   *  \brief get a pointer on SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SparseMat*
   */
  virtual const SparseMat* getSparsePtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** \fn mapped getMap() const  = 0
   *  \brief get mapped matrix
   *  \useless for MySimpleMatrix
   *  \return a mapped
   */
  virtual const mapped getMap(void) const = 0;

  /** \fn void getBlock(unsigned int numCow, unsigned int numCol, MySiconosMatrix& block) const  = 0
   *  \brief get block corresponding to lines given in numRow and columns in numCol
   *  \param unsigned int, row index
   *  \param unsigned int, col index
   *  \param a MySiconosMatrix (in-out paramater)
   */
  virtual void getBlock(unsigned int, unsigned int, MySiconosMatrix&) const = 0;

  /** \fn const std::deque<bool> getBlockAllocated() const = 0
   *   \brief get std::deque of bool
   *   \useless for MySimpleMatrix
   *   \return a std::deque<bool>
   */
  virtual const std::deque<bool> getBlockAllocated(void) const = 0;

  /** \fn void getRow(unsigned int index, MySimpleVector& vOut) const = 0
   *  \brief get row index of current matrix and save it into vOut
   *  \param unsigned int: index of required line
   *  \param ref to MySimpleVector: in-out parameter
   */
  virtual void getRow(unsigned int, MySimpleVector&) const = 0;

  /** \fn void getCol(unsigned int index, MySimpleVector& vOut) const = 0
   *  \brief get column index of current matrix and save it into vOut
   *  \param unsigned int: index of required column
   *  \param ref to MySimpleVector: in-out parameter
   */
  virtual void getCol(unsigned int, MySimpleVector&) const = 0;

  /** \fn void setNum(unsigned int n)
   *  \brief set the attribute num of current matrix with n
   */
  virtual void setNum(unsigned int) = 0;

  /** \fn void setRow(unsigned int row, const MySimpleVector &v)
   *  \brief set line row of the current matrix with vector v
   *  \param an unsigned int and a MySimpleVector
   */
  virtual void setRow(unsigned int, const MySimpleVector&) = 0;

  /** \fn void setCol(unsigned int col, const MySimpleVector &v)
   *  \brief set column col of the current matrix with vector v
   *  \param an unsigned int and a MySimpleVector
   */
  virtual void setCol(unsigned int, const MySimpleVector&) = 0;

  /** \fn  unsigned int size1 (void) const = 0
   *  \brief get the number of rows of the matrix
   *  \exception SiconosMatrixException
   *  \return the number of rows of the matrix
   */
  virtual unsigned int size1(void)const = 0;

  /** \fn  unsigned int size2 (void) const = 0
   *  \brief get the number of columns of the matrix
   *  \exception SiconosMatrixException
   *  \return the number of columns of the matrix
   */
  virtual unsigned int size2(void) const = 0;

  /** \fn  void resize (unsigned int nbrow, unsigned int nbcol, unsigned int lower=0, unsigned int upper=0,  bool val = true) = 0
   *  \brief resize the matrix with nbrow rows and nbcol columns, upper and lower are only useful for BandedMatrix .
   *   The existing elements of the matrix are preseved when specified.
   *  \exception SiconosMatrixException
   *  \param 2 unsigned int: number of rows and columns
   *  \param 2 unsigned int: for banded matrices
   */
  virtual void resize(unsigned int, unsigned int, unsigned int lower = 0, unsigned int upper = 0, bool = true) = 0;

  /** \fn const double normInf() const = 0;
   *  \brief compute the infinite norm of the matrix
   *  \return a double
   */
  virtual const double normInf(void) const = 0;

  /** \fn void display() = 0;
   *  \brief display data on standard output
   */
  virtual void display(void) const = 0;

  /** \fn void zero() = 0;
   *  \brief sets all the values of the matrix to 0.0
   */
  virtual void zero(void) = 0;

  /** \fn void eye() = 0;
   *  \brief set an identity matrix
   */
  virtual void eye(void) = 0;


  // --- MATRICES HANDLING AND OPERATORS ---

  /** \fn void blockMatrixCopy( MySiconosMatrix &blockMat, unsigned int, unsigned int) = 0
   *  \brief copy the blockmatrix "blockMat" in the matrix "mat" at the position (xPos, yPos)
   *                    [0, 0, 0, 0].blockMatrixCopy([1], 0, 2) => mat = [0 0 1 0]
   *  \param MySiconosMatrix& : the block matrix to copy in the current matrix
   *  \param unsigned int : the line position to start the copy of the blockmatrix
   *  \param unsigned int : the column position to start the copy of the blockmatrix
   */
  virtual void blockMatrixCopy(const MySiconosMatrix&, unsigned int, unsigned int) = 0;

  // Note: in the following functions, row and col are general;
  // that means that for a MySimpleMatrix m, m(i,j) is index (i,j) element but
  // for a MyBlockMatrix w that contains 2 MySiconosMatrix of size 3
  // w(1, 4) corresponds to the element (1,1) of the second matrix.


  /** \fn double& operator() (unsigned int row, unsigned int col) = 0
   *  \brief get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  virtual double& operator()(unsigned int , unsigned int) = 0;

  /** \fn double operator() (unsigned int row, unsigned int col)const
   *  \brief get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  virtual double operator()(unsigned int , unsigned int) const = 0;

  /** \fn assignment operator
   *  \param MySiconosMatrix : the matrix to be copied
   */
  virtual const MySiconosMatrix& operator  =(const MySiconosMatrix&) = 0;

  /** \fn operator /= (double) = 0
   *  \param double, a scalar
   */
  virtual const MySiconosMatrix& operator /=(double) = 0;

  /** \fn operator /= (int) = 0
   *  \param int, a scalar
   */
  virtual const MySiconosMatrix& operator /=(int) = 0;

  /** \fn operator += (const MySiconosMatrix&) = 0
   *  \param MySiconosMatrix : a matrix to add
   */
  virtual const MySiconosMatrix& operator +=(const MySiconosMatrix&) = 0;

  /** \fn operator -= (const MySiconosMatrix&) = 0
   *  \param MySiconosMatrix : a matrix to subtract
   */
  virtual const MySiconosMatrix& operator -=(const MySiconosMatrix&) = 0;

  /** \fn operator *= (double) = 0
   *  \param double, a scalar
   */
  virtual const MySiconosMatrix& operator *=(double) = 0;

  /** \fn operator *= (int) = 0
   *  \param int, a scalar
   */
  virtual const MySiconosMatrix& operator *=(int) = 0;

};

#endif
