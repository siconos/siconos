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

  /**\fn SetIsBlockMatrix (bool)
  * \brief Set the value of isBlockMatrix.
  * \param a bool to set isBlockMatrix (optional, default = false)
  */
  inline void       SetIsBlock(bool val)
  {
    isBlockMatrix = val;
  }

public:

  /**\fn bool isBlock ()
  * \brief true if the matrix is block else false.
  * \return a bool.*/
  inline  bool      isBlock(void)const
  {
    return isBlockMatrix;
  }

  /**\fn ~MySiconosMatrix ()
  * \brief Destructor. */
  virtual       ~MySiconosMatrix(void)          = 0;

  // ************ GETTERS/SETTERS ***************


  /** \fn int getNum() const
  *  \brief get the attribute num of current matrix
  * \return an int.
  */
  virtual int       GetNum(void)const         = 0;

  /** \fn DenseMat getDense(int row = 0, int col = 0)
  *  \brief get DenseMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a DenseMat
  */
  virtual const DenseMat    GetDense(int = 0, int = 0)const       = 0;

  /** \fn TriangMat getTriang(int row = 0, int col = 0)
  *  \brief get TriangMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a TriangMat
  */
  virtual const TriangMat   GetTriang(int = 0, int = 0)const       = 0;

  /** \fn SymMat getSym(int row = 0, int col = 0)
  *  \brief get SymMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a SymMat
  */
  virtual const SymMat    GetSym(int = 0, int = 0)const     = 0;

  /** \fn BandedMat getBanded(int row = 0, int col = 0)
  *  \brief get BandedMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a BandedMat
  */
  virtual const BandedMat     GetBanded(int = 0, int = 0)const      = 0;

  /** \fn SparseMat getSparse(int row = 0, int col = 0)
  *  \brief get SparseMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a SparseMat
  */
  virtual const SparseMat     GetSparse(int = 0, int = 0)const      = 0;

  /** \fn DenseMat* getDensePtr(int row = 0, int col = 0)
  *  \brief get a pointer on DenseMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a DenseMat*
  */
  virtual const DenseMat*   GetDensePtr(int = 0, int = 0)const       = 0;

  /** \fn TriangMat* getTriangPtr(int row = 0, int col = 0)
  *  \brief get a pointer on TriangMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a TriangMat*
  */
  virtual const TriangMat*  GetTriangPtr(int = 0, int = 0)const       = 0;

  /** \fn SymMat* getSymPtr(int row = 0, int col = 0)
  *  \brief get a pointer on SymMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a SymMat*
  */
  virtual const SymMat*     GetSymPtr(int = 0, int = 0)const     = 0;

  /** \fn BandedMat* getBandedPtr(int row = 0, int col = 0)
  *  \brief get a pointer on BandedMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a BandedMat*
  */
  virtual const BandedMat*    GetBandedPtr(int = 0, int = 0)const       = 0;

  /** \fn SparseMat* getSparsePtr(int row = 0, int col = 0)
  *  \brief get a pointer on SparseMat matrix
  *  \param an int, position of the block (row) - Useless for MySimpleMatrix
  *  \param an int, position of the block (column) - Useless for MySimpleMatrix
  *  \return a SparseMat*
  */
  virtual const SparseMat*    GetSparsePtr(int = 0, int = 0)const       = 0;

  /** \fn mapped getMap()
  *  \brief get mapped matrix
  *  \useless for MySimpleMatrix
  *  \return a mapped
  */
  virtual const mapped    GetMap(void)const       = 0;

  /** \fn void getBlock(int numCow, int numCol, MySiconosMatrix& block)
  *  \brief get block corresponding to lines given in numRow and columns in numCol
  *  \param 2 int for indexes and a MySiconosMatrix (in-out paramater)
  */
  virtual void      GetBlock(int, int, MySiconosMatrix&)const   = 0;

  /** \fn const std::deque<bool> GetBlockAllocated()const
  *   \brief get std::deque of bool
  *   \useless for MySimpleMatrix
  *   \return a std::deque<bool>
  */
  virtual const std::deque<bool>  GetBlockAllocated(void)const        = 0;

  /** \fn void getRow(int index, MySimpleVector& vOut) const
  *  \brief get row index of current matrix and save it into vOut
  *  \param int: index of required line
  *  \param ref to MySimpleVector: in-out parameter
  */
  virtual void      GetRow(int, MySimpleVector&)const     = 0;

  /** \fn void getCol(int index, MySimpleVector& vOut) const
  *  \brief get column index of current matrix and save it into vOut
  *  \param int: index of required column
  *  \param ref to MySimpleVector: in-out parameter
  */
  virtual void      GetCol(int, MySimpleVector&)const     = 0;

  /** \fn void SetNum(int n)
  *  \brief set the attribute num of current matrix with n
  */
  virtual void      SetNum(int)         = 0;

  /** \fn void setRow(int row, const MySimpleVector &v)
  *  \brief set line row of the current matrix with vector v
  *  \param an int and a MySimpleVector
  */
  virtual void      SetRow(int, const MySimpleVector&)    = 0;

  /** \fn void setCol(int col, const MySimpleVector &v)
  *  \brief set column col of the current matrix with vector v
  *  \param an int and a MySimpleVector
  */
  virtual void      SetCol(int, const MySimpleVector&)    = 0;

  /** \fn  int size1 (void)const
   *  \brief get the of rows of the matrix
   *  \exception SiconosMatrixException
   *  \return the number of rows of the matrix
   */
  virtual int     size1(void)const         = 0;

  /** \fn  int size2 (void)const
  *  \brief get the of columns of the matrix
  *  \exception SiconosMatrixException
  *  \return the number of columns of the matrix
  */
  virtual int     size2(void)const         = 0;

  /** \fn  void resize (int nbrow, int nbcol, int lower=0, int upper=0,  bool val = true)const
  *  \brief resize the matrix with nbrow rows and nbcol columns, upper and lower are only useful for BandedMatrix . The existing elements of the matrix are preseved when specified.
  *  \exception SiconosMatrixException
  */
  virtual void      resize(int, int, int lower = 0, int upper = 0, bool = true)      = 0;

  /** \fn const double normInf() const;
  *  \brief compute the infinite norm of the matrix
  *  \return a double
  */
  virtual const double    normInf(void)const         = 0;

  /** \fn void display();
  *  \brief display data on standard output
  */
  virtual void      display(void)const         = 0;

  /** \fn void zero();
  *  \brief sets all the values of the matrix to 0.0
  */
  virtual void      zero(void)          = 0;

  /** \fn void eye();
  *  \brief set an identity matrix
  */
  virtual void      eye(void)          = 0;


  // --- MATRICES HANDLING AND OPERATORS ---

  /** \fn void blockMatrixCopy( MySiconosMatrix &blockMat, int, int)
  *  \brief copy the blockmatrix "blockMat" in the matrix "mat" at the position (xPos, yPos)
  *                    [0, 0, 0, 0].blockMatrixCopy([1], 0, 2) => mat = [0 0 1 0]
  *  \param MySiconosMatrix& : the block matrix to copy in the current matrix
  *  \param int : the line position to start the copy of the blockmatrix
  *  \param int : the column position to start the copy of the blockmatrix
  */
  virtual void      BlockMatrixCopy(const MySiconosMatrix&, int, int)  = 0;

  // Note: in the following functions, row and col are general;
  // that means that for a MySimpleMatrix m, m(i,j) is index (i,j) element but
  // for a MyBlockMatrix w that contains 2 MySiconosMatrix of size 3
  // w(1, 4) corresponds to the element (1,1) of the second matrix.


  /** \fn double& operator() (int row, int col)
  *  \brief get or set the element matrix[i,j]
  *  \param an int i
  *  \param an int j
  *  \exception SiconosMatrixException
  *  \return the element matrix[i,j]
  */
  virtual double&     operator()(int , int)        = 0;

  /** \fn double operator() (int row, int col)const
  *  \brief get or set the element matrix[i,j]
  *  \param an int i
  *  \param an int j
  *  \exception SiconosMatrixException
  *  \return the element matrix[i,j]
  */
  virtual double      operator()(int , int)const     = 0;

  /** \fn assignment operator
  *  \param MySiconosMatrix : the matrix to be copied
  */
  virtual const MySiconosMatrix&  operator  = (const MySiconosMatrix&)    = 0;

  /** \fn operator /= (double)
  *  \param double, a scalar
  */
  virtual const MySiconosMatrix&  operator /= (double)        = 0;

  /** \fn operator /= (int)
  *  \param int, a scalar
  */
  virtual const MySiconosMatrix&  operator /= (int)         = 0;

  /** \fn operator += (const MySiconosMatrix&)
  *  \param MySiconosMatrix : a matrix to add
  */
  virtual const MySiconosMatrix&  operator += (const MySiconosMatrix&)    = 0;

  /** \fn operator -= (const MySiconosMatrix&)
  *  \param MySiconosMatrix : a matrix to subtract
  */
  virtual const MySiconosMatrix&  operator -= (const MySiconosMatrix&)    = 0;

  /** \fn operator *= (double)
  *  \param double, a scalar
  */
  virtual const MySiconosMatrix&  operator *= (double)        = 0;

  /** \fn operator *= (int)
  *  \param int, a scalar
  */
  virtual const MySiconosMatrix&  operator *= (int)         = 0;

};

#endif
