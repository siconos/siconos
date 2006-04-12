/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
/** \class SiconosMatrix: abstract class
 *   \brief this class provides an interface for matrices in Siconos
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.4.
 *   \date (Creation) 05/19/2004
 *
 *   Matrices can be either block or simple. See derived classes for details.
 *    Simple: a LaGenMatDouble from Lapack++
 *    Block: a map(stl) of Simple
 *
 */

#ifndef __SiconosMatrix__
#define __SiconosMatrix__

#include "SiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrixException.h"
#include "check.h"
#include <lapack++.h>
#include <string>
#include <fstream>
#include <iostream>

const bool printVerbose = true;

// This constant set the maximum allowed size for displaying a Matrix on standard output (see function display())
const unsigned int MAXSIZEFORDISPLAY = 10;
class SimpleVector;
class SimpleMatrix;

class SiconosMatrix
{
protected:

  /**\var isBlockMatrix
   *  bool to check this matrix type; true if block else false. */
  bool isBlockMatrix;

  /** \fn static void verbose(std::string msg)
   *  \brief print on the screen a message "printVerbose" static variable is true
   *  \param std::string : the message to print
   */
  static void verbose(const std::string& msg);

  /** \fn SiconosMatrix (const bool & = false)
   *  \brief default contructor
   *  \param a bool to set isBlockMatrix (optional, default = false)
   */
  SiconosMatrix(const bool & = false);

public:

  /** \fn ~SiconosMatrix ()
   *  \brief destructor
   */
  virtual ~SiconosMatrix();

  /** \fn  unsigned int size (const unsigned int&)
   *  \brief get the dimension of the matrix
   *  \param the dimension to get (0 to get the number of rows or 1 to get the number of columns)
   *  \exception SiconosMatrixException
   *  \return the a dimension of the matrix
   */
  virtual unsigned int size(const unsigned int&) const = 0;

  /** \fn bool isSquare()
   *  \brief determines if the matrix is square
   *  \return true if the matrix is square
   */
  virtual bool isSquare() const = 0;

  /** \fn bool isInversed()
   *  \brief determines if the matrix has been inversed
   *  \return true if the matrix is inversed
   */
  virtual bool isInversed() const  = 0;

  /** \fn bool isfactorized()
   *  \brief determines if the matrix has been factorized
   *  \return true if the matrix is factorized
   */
  virtual bool isFactorized() const  = 0;

  /** \fn bool isBlock()
   *  \brief true if the matrix is block else false
   *  \return a bool
   */
  inline bool isBlock() const
  {
    return isBlockMatrix;
  };

  // --- GETTERS/SETTERS ---

  /** \fn LaGenMatDouble getLaGenMatDouble(const int& row = 0, const int& col = 0)
   *  \brief get LaGenMatDouble matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a LaGenMatDouble
   */
  virtual const LaGenMatDouble getLaGenMatDouble(const unsigned int& = 0, const unsigned int& = 0) const = 0;

  /** \fn LaGenMatDouble* getLaGenMatDoubleRef(const int& row = 0, const int& col = 0)
   *  \brief get LaGenMatDouble matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a LaGenMatDouble
   */
  virtual const LaGenMatDouble* getLaGenMatDoubleRef(const unsigned int& = 0, const unsigned int& = 0) const = 0;

  /** \fn setValue(const unsigned int& row, const unsigned int& col, const double& d)
   * \brief set the element matrix[row, col]
   * \param an integer col
   * \param an integer row
   * \param a double d : the value.
   */
  virtual void setValue(const unsigned int & row, const unsigned int& col, const double& d) = 0;

  /** \fn setValue(const LaGenMatDouble&, const unsigned int& row = 0, const unsigned int& col = 0)
   * \brief set mat of a SimpleMatrix, of mat-block(row,col) of a BlockMatrix
   * \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   * \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   * \param a LaGenMatDouble
   */
  virtual void setValue(const LaGenMatDouble&, const unsigned int& = 0, const unsigned int& = 0) = 0;

  /** \fn const double getValue(const int& row, const int& col) const
   * \brief get the element matrix[row, col]
   * \param an integer col
   * \param an integer row
   * \return a double d : the value.
   */
  inline const double getValue(const unsigned int& row, const unsigned int& col)
  {
    return (*this)(row, col);
  }

  /** \fn void setRow(const unsigned int& row, const SiconosVector &v)
   *  \brief set line row of the current matrix with vector v
   *  \param an int and a SiconosVector
   */
  virtual void setRow(const unsigned int& , const SiconosVector &) = 0;

  /** \fn void getRow(const int& index, const SimpleVector& vOut) const
   *  \brief get row index of current matrix and save it into vOut
   *  \param int: index of required line
   *  \param ref to SimpleVector: in-out parameter
   */
  virtual void getRow(const unsigned int& i, const SimpleVector&) const = 0;

  /** \fn void getCol(const int& index, const SimpleVector& vOut) const
   *  \brief get column index of the current matrix and save it into vOut
   *  \param int: index of required column
   *  \param ref to SimpleVector: in-out parameter
   */
  virtual void getCol(const unsigned int&, const SimpleVector&) const = 0;

  /** \fn void setCol(const unsigned int& col, const SiconosVector &v)
   *  \brief set column col of the current matrix with vector v
   *  \param an int and a SiconosVector
   */
  virtual void setCol(const unsigned int& , const SiconosVector &) = 0;

  /** \fn void getBlock(const std::vector<unsigned int>& indexList, SiconosMatrix&)
   *  \brief get a block of a matrix, which position is defined by index_list
   *  \param vector<unsigned int> for indexes and a SiconosMatrix (in-out paramater)
   */
  virtual void getBlock(const std::vector<unsigned int>&, SiconosMatrix&) const = 0;

  /** \fn void getBlock(const vector<unsigned int>& indexRow, const vector<unsigned int>& indexCol, SiconosMatrix& block)
   *  \brief get block corresponding to lines given in indexRow and columns in indexCol
   *  \param 2 vector<unsigned int> for indexes and a SiconosMatrix (in-out paramater)
   */
  virtual void getBlock(const std::vector<unsigned int>& , const std::vector<unsigned int>&, SiconosMatrix&) const = 0;

  /** \fn SiconosMatrix* getBlockPtr(const unsigned int& row, const unsigned int& col)
   *  \brief get block at position row-col if BlockMatrix, else if SimpleMatrix return this
   *  \param unsigned int row
   *  \param unsigned int col
   */
  virtual SiconosMatrix* getBlockPtr(const unsigned int& = 0, const unsigned int& = 0) = 0;

  /** \fn double* getArray(const unsigned int& i = 0, const unsigned int&  j = 0)
   *  \brief return the adress of the array of double values of the matrix ( for block(i,j) if this is a block matrix)
   *  \param: row position for the required block
   *  \param: col position for the required block
   *  \return double* : the pointer on the double array
   */
  virtual double* getArray(const unsigned int& = 0, const unsigned int& = 0) = 0;

  /** \fn bool read(std::string fileName, std::string mode = BINARY)
   *  \brief read the matrix in a file
   *  \param std::string fileName : the input file
   *  \param std::string mode : ASCII or BINARY (binary is the default mode)
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  virtual bool read(const std::string& , const std::string& = "binary") = 0;

  /** \fn bool write(std::string fileName, std::string mode = BINARY)
   *  \brief write the matrix in a file
   *  \param std::string fileName : the output file
   *  \param std::string mode : ASCII or BINARY (binary is the default mode
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  virtual bool write(const std::string& fileName, const std::string& mode = "binary") const  = 0;

  /** \fn bool rawWrite(std::string fileName, std::string mode = BINARY)
   *  \brief write the matrix in a file without dimensions output
   *  \param std::string fileName : the output file
   *  \param std::string mode : ASCII or BINARY (binary is the default mode
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  virtual bool rawWrite(const std::string& fileName, const std::string& mode = "binary") const  = 0;

  /** \fn void zero();
   *  \brief sets all the values of the matrix to 0.0
   */
  virtual void zero() = 0;

  /** \fn void eye();
   *  \brief set an identity matrix
   */
  virtual void eye() = 0;

  /** \fn void display();
   *  \brief display data on standard output
   */
  virtual void display() const = 0;

  // --- MATRICES HANDLING AND OPERATORS ---

  /** \fn SimpleMatrix multTranspose(const SiconosMatrix B)
   *  \brief compute A*Bt
   *  \param SiconosMatrix B
   *  \return SimpleMatrix : the result of the multiplication
   */
  virtual SimpleMatrix multTranspose(const SiconosMatrix &) = 0;

  /** \fn void blockMatrixCopy( SiconosMatrix &blockMat, const int&, const int&)
   *  \brief copy the blockmatrix "blockMat" in the matrix "mat" at the position (xPos, yPos)
   *      blockMatrixCopy([1], [0 0 0 0], 0, 2) => mat = [0 0 1 0]
   *  \param SiconosMatrix& : the block matrix to copy in the current matrix
   *  \param int : the line position to start the copy of the blockmatrix
   *  \param int : the column position to start the copy of the blockmatrix
   */
  virtual void blockMatrixCopy(const SiconosMatrix &, const unsigned int&, const unsigned int&) = 0;

  // Io-stream operators
  friend std::istream& operator >> (std::istream& i, SiconosMatrix& m);
  friend std::ostream& operator << (std::ostream& o, SiconosMatrix& m);

  /** \fn operator() (const int& row, const int& col)= 0
   *  \brief get or set the element matrix[i,j]
   *  \param an integer i
   *  \param an integer j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  virtual double& operator()(const int& row, const int& col) = 0;

  /** \fn operator () (const unsigned int& row, const unsigned int& col) = 0;
   *  \brief get or set the element matrix[i,j]
   *  \param an integer i
   *  \param an integer j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  virtual double& operator()(const unsigned int& row, const unsigned int& col) = 0;
  virtual double& operator()(const unsigned int& row, const unsigned int& col) const = 0;

  /** \fn assignment operator
   *  \param SiconosMatrix : the matrix to be copied
   */
  virtual SiconosMatrix& operator = (const SiconosMatrix& m) = 0;

  /** \fn operator +=
  *  \param SiconosMatrix : a matrix to add
  */
  virtual SiconosMatrix& operator+=(const SiconosMatrix  &)  = 0;

  /** \fn operator -=
  *  \param SiconosMatrix : a matrix to subtract
  */
  virtual SiconosMatrix& operator-=(const SiconosMatrix  &)  = 0;

  /** \fn operator *=
  *  \param double, a scalar
  */
  virtual SiconosMatrix& operator*=(const double  &)  = 0;

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the matrix
   *  \return a double
   */
  virtual const double normInf() const = 0;

  // --- COMPUTING WITH MATRICES  ---

  /** \fn SiconosMatrix linearSolve(const SiconosMatrix & B,SiconosMatrix &X )
   *  \brief compute the vector x in  Ax=B
   *  \param SiconosMatrix & B
   *  \param SiconosMatrix X
   */
  virtual void linearSolve(const SiconosMatrix & , SiconosMatrix&) = 0;

  /** \fn  PLUFactorizationInPlace(void);
   *  \brief Compute the LU factorization with Partial pivoting.  The result is returned in this (InPlace)
   */
  virtual void PLUFactorizationInPlace(void) = 0;

  /** \fn  SiconosMatrix  PLUInverseInPlace(void);
   *  \brief  use the LU factorization with Partial pivoting of the matrix this (or the matrix this itself if it is triangular)
   *       to compute the inverse matrix.  The result is returned in this (InPlace)
   */
  virtual void  PLUInverseInPlace(void) = 0;

  /** \fn SiconosMatrix  PLUForwardBackward(SiconosMatrix &B);
   *  \brief use the LU factorization with Partial pivoting of the matrix this (or the matrix this itself if it is triangular)
   *  to solve A x = B by forward or backward substitution;
   *  \param the RHS matrix b  which contains the result x
   */
  virtual void  PLUForwardBackwardInPlace(SiconosMatrix &B) = 0 ;

  /** \fn SiconosVector  PLUForwardBackward(SiconosVector &B);
   *  \brief use the LU factorization with Partial pivoting of the matrix this (or the matrix this itself if it is triangular)
   *  to solve A x = B by forward or backward substitution;
   *  \param the RHS vector b which contains the result x
   */
  virtual void   PLUForwardBackwardInPlace(SiconosVector &B) = 0;
};

#endif // __SiconosMatrix__
