/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
/** \class SimpleMatrix
 *   \brief This class is an encapsulation of the Lapack++ class managing vmatrices of double.
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.3.
 *   \date (Creation) 05/19/2004
 *
 *
 *
 *
 * SimpleMatrix is used in the platform to store matrices (mathematical object).
 *
 */

#ifndef __SimpleMatrix__
#define __SimpleMatrix__

#include "SiconosMatrix.h"

class SimpleVector;

class SimpleMatrix: public SiconosMatrix
{
private:

  /** \var LaGenMatDouble mat;
   * \brief The Matrix Lapack ++, LaGenMatDouble, encapsulated.
   */
  LaGenMatDouble mat;

  /** \var LaVectorLongInt* ipiv;
   * \brief Pointer to a  LaVectorLongInt for the partial pivoting in PLUFactorization.
   */
  LaVectorLongInt* ipiv;

  /** \var  bool isPLUFactorized;
   * \brief  Boolean = true if the Matrix is PLU Factorized
   */
  bool isPLUFactorized;

  /** \var  bool isPLUInversed;
   * \brief  Boolean = true if the Matrix is Inversed in Place
   */
  bool isPLUInversed;

  /** \fn SimpleMatrix ()
   *  \brief contructor
   *  \return SimpleMatrix
   */
  SimpleMatrix();

public:

  /** \fn SimpleMatrix (const SiconosMatrix&  m)
   *  \brief copy contructor
   *  \param SiconosMatrix
   */
  SimpleMatrix(const SiconosMatrix&);

  /** \fn SimpleMatrix (const SimpleMatrix&  m)
   *  \brief copy contructor
   *  \param SimpleMatrix
   */
  SimpleMatrix(const SimpleMatrix&);

  /** \fn SimpleMatrix (const unsigned int& row , const unsigned int& col)
   *  \brief contructor with the dimension
   *  \param an integer to set the number of row
   *  \param an integer to set the number of column
   */
  SimpleMatrix(const unsigned int&, const unsigned int&);

  /** \fn SimpleMatrix (const LaGenMatDouble m)
   *  \brief contructor with a LaGenMatDouble
   *  \param a LaGenMatDouble
   */
  SimpleMatrix(const LaGenMatDouble&);

  /** \fn SimpleMatrix (const LaVectorDouble v, unsigned int& row, unsigned int& col);
   *  \brief contructor with a LaVectorDouble
   *  \param a LaVectorDouble
   *  \param an integer to set the number of row
   *  \param an integer to set the number of column
   */
  SimpleMatrix(const LaVectorDouble&, const unsigned int&, const unsigned int&);

  /** \fn SimpleMatrix (std::string file, bool ascii)
   *  \brief contructor with an input file
   *  \param a std::string which contain the file path
   *  \param a boolean to indicate if the file is in ascii
   */
  SimpleMatrix(const std::string&, const bool&);

  /** \fn ~SimpleMatrix ()
   *  \brief destructor
   */
  ~SimpleMatrix();

  /** \fn  unsigned int size (const unsigned int&)
   *  \brief get the dimension of the matrix
   *  \param the dimension to get (0 to get the number of rows or 1 to get the number of columns)
   *  \exception SiconosMatrixException
   *  \return the a dimension of the matrix
   */
  unsigned int size(const unsigned int&) const;

  /** \fn bool isSquare()
   *  \brief determines if the matrix is square
   *  \return true if the matrix is square
   */
  bool isSquare() const;

  /** \fn bool isInversed()
   *  \brief determines if the matrix has been inversed
   *  \return true if the matrix is inversed
   */
  inline bool isInversed() const
  {
    return isPLUInversed;
  }

  /** \fn bool isfactorized()
   *  \brief determines if the matrix has been factorized
   *  \return true if the matrix is factorized
   */
  inline bool isFactorized() const
  {
    return isPLUFactorized;
  }

  // --- GETTERS/SETTERS ---

  /** \fn LaGenMatDouble getLaGenMatDouble(const int& row = 0, const int& col = 0)
   *  \brief get LaGenMatDouble matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a LaGenMatDouble
   */
  const LaGenMatDouble getLaGenMatDouble(const unsigned int& = 0, const unsigned int& = 0) const;

  /** \fn setValue(const unsigned int& row, const unsigned int& col, const double& d)
   * \brief set the element matrix[row, col]
   * \param an integer col
   * \param an integer row
   * \param a double d : the value.
   */
  inline void setValue(const unsigned int & row, const unsigned int& col, const double& d)
  {
    (*this)(row, col) = d;
  }

  /** \fn setValue(const LaGenMatDouble&, const unsigned int& row = 0, const unsigned int& col = 0)
   * \brief set mat of a SimpleMatrix, of mat-block(row,col) of a BlockMatrix
   * \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   * \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   * \param a LaGenMatDouble
   */
  void setValue(const LaGenMatDouble&, const unsigned int& = 0, const unsigned int& = 0);

  /** \fn void setRow(const unsigned int& row, const SiconosVector &v)
   *  \brief set line row of the current matrix with vector v
   *  \param an int and a SiconosVector
   */
  void setRow(const unsigned int& , const SiconosVector &);

  /** \fn void getRow(const int& index, const SimpleVector& vOut) const
   *  \brief get row index of current matrix and save it into vOut
   *  \param int: index of required line
   *  \param ref to SimpleVector: in-out parameter
   */
  void getRow(const unsigned int& i, const SimpleVector&) const;

  /** \fn void getCol(const int& index, const SimpleVector& vOut) const
   *  \brief get column index of the current matrix and save it into vOut
   *  \param int: index of required column
   *  \param ref to SimpleVector: in-out parameter
   */
  void getCol(const unsigned int&, const SimpleVector&) const;

  /** \fn void setCol(const unsigned int& col, const SiconosVector &v)
   *  \brief set column col of the current matrix with vector v
   *  \param an int and a SiconosVector
   */
  void setCol(const unsigned int& , const SiconosVector &);

  /** \fn void getBlock(const std::vector<unsigned int>& indexList, SiconosMatrix&)
   *  \brief get a block of a matrix, which position is defined by index_list
   *  \param vector<unsigned int> for indexes and a SiconosMatrix (in-out paramater)
   */
  void getBlock(const std::vector<unsigned int>&, SiconosMatrix&) const;

  /** \fn void getBlock(const vector<unsigned int>& indexRow, const vector<unsigned int>& indexCol, SiconosMatrix& block)
   *  \brief get block corresponding to lines given in indexRow and columns in indexCol
   *  \param 2 vector<unsigned int> for indexes and a SiconosMatrix (in-out paramater)
   */
  void getBlock(const std::vector<unsigned int>& , const std::vector<unsigned int>&, SiconosMatrix&) const;

  /** \fn double* getArray(const unsigned int& i = 0, const unsigned int&  j = 0)
   *  \brief return the adress of the array of double values of the matrix
   *  \param: row position for the required block ->useless for SimpleMatrix
   *  \param: col position for the required block ->useless for SimpleMatrix
   *  \return double* : the pointer on the double array
   */
  inline double* getArray(const unsigned int& = 0, const unsigned int& = 0)
  {
    return mat.addr();
  }

  /** \fn bool read(std::string fileName, std::string mode = BINARY)
   *  \brief read the matrix in a file
   *  \param std::string fileName : the input file
   *  \param std::string mode : ASCII or BINARY (binary is the default mode)
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool read(const std::string& , const std::string& = "binary");

  /** \fn bool write(std::string fileName, std::string mode = BINARY)
   *  \brief write the matrix in a file
   *  \param std::string fileName : the output file
   *  \param std::string mode : ASCII or BINARY (binary is the default mode
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool write(const std::string& fileName, const std::string& mode = "binary") const ;

  /** \fn bool rawWrite(std::string fileName, std::string mode = BINARY)
   *  \brief write the matrix in a file without dimensions output
   *  \param std::string fileName : the output file
   *  \param std::string mode : ASCII or BINARY (binary is the default mode
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool rawWrite(const std::string& fileName, const std::string& mode = "binary") const ;

  /** \fn void zero();
   *  \brief sets all the values of the matrix to 0.0
   */
  void zero();

  /** \fn void eye();
   *  \brief set an identity matrix
   */
  void eye();

  /** \fn void display();
   *  \brief display data on standard output
   */
  void display() const;

  // --- MATRICES HANDLING AND OPERATORS ---

  /** \fn SimpleMatrix multTranspose(const SiconosMatrix B)
   *  \brief compute A*Bt
   *  \param SiconosMatrix B
   *  \return SimpleMatrix : the result of the multiplication
   */
  SimpleMatrix multTranspose(const SiconosMatrix &);

  /** \fn void blockMatrixCopy( SiconosMatrix &blockMat, const int&, const int&)
   *  \brief copy the blockmatrix "blockMat" in the matrix "mat" at the position (xPos, yPos)
   *      blockMatrixCopy([1], [0 0 0 0], 0, 2) => mat = [0 0 1 0]
   *  \param SiconosMatrix& : the block matrix to copy in the current matrix
   *  \param int : the line position to start the copy of the blockmatrix
   *  \param int : the column position to start the copy of the blockmatrix
   */
  void blockMatrixCopy(const SiconosMatrix &, const unsigned int&, const unsigned int&);

  /** \fn operator (int row, int col)
   *  \brief get or set the element matrix[i,j]
   *  \param an integer i
   *  \param an integer j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double& operator()(const int& row, const int& col);

  /** \fn operator (int row, int col)
   *  \brief get or set the element matrix[i,j]
   *  \param an integer i
   *  \param an integer j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double& operator()(const unsigned int& row, const unsigned int& col);
  double& operator()(const unsigned int& row, const unsigned int& col) const;

  /** \fn assignment operator
   *  \param SiconosMatrix : the matrix to be copied
   */
  SimpleMatrix& operator = (const SiconosMatrix&);

  /** \fn assignment operator
   *  \param SimpleMatrix : the matrix to be copied
   */
  SimpleMatrix& operator = (const SimpleMatrix&);

  /** \fn operator +=
  *  \param  SimpleMatrix: a matrix to add
  */
  SimpleMatrix& operator+=(const SiconosMatrix  &) ;

  /** \fn operator -=
  *  \param SimpleMatrix : a matrix to subtract
  */
  SimpleMatrix& operator-=(const SiconosMatrix  &);

  /** \fn operator *=
  *  \param double, a scalar
  */
  SimpleMatrix& operator*=(const double  &);

  /** \fn operator ==
   * \brief: A==B when (A-B).normInf()<tolerance
   * \param a SiconosMatrix
   */
  friend bool operator==(const SiconosMatrix&, const SiconosMatrix&);

  /** \fn operator * (const SimpleMatrix& m1, const SimpleMatrix& m2);
   *  \brief multiplication of two matrices
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   */
  friend SimpleMatrix operator * (const SimpleMatrix&, const SimpleMatrix&);

  /** \fn operator * (const SimpleMatrix& m1, const double& d);
   *  \brief multiplication of a matrix by a double
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix operator * (const SiconosMatrix&, const double&);

  /** \fn operator * (const double& d, const SimpleMatrix& m1);
   *  \brief multiplication of a matrix by a double
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix operator * (const double&, const SiconosMatrix&);

  /** \fn operator / (const SimpleMatrix& m1, const double d);
   *  \brief division of the matrix by a double
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the double d = 0
   */
  friend SimpleMatrix operator / (const SiconosMatrix&, const double&);

  /** \fn operator + (const SiconosMatrix& m1, const SiconosMatrix& m2);
   *  \brief Addition of two matrices
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   */
  friend SimpleMatrix operator + (const SiconosMatrix&, const SiconosMatrix&);

  /** \fn operator - (const SiconosMatrix& m1, const SiconosMatrix& m2);
   *  \brief subtraction of two matrices
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   */
  friend SimpleMatrix operator - (const SiconosMatrix&, const SiconosMatrix&);

  /** \fn pow(const SimpleMatrix& m1,  const unsigned int&);
   *  \brief compute the power of the matrix (!)
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the power < 0
   */
  friend SimpleMatrix pow(const SimpleMatrix&, const unsigned int&);

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the matrix
   *  \return a double
   */
  const double normInf() const;

  // --- COMPUTING WITH MATRICES  ---

  /** \fn void linearSolve(const SiconosMatrix & B, SimpleMatrix & X)
   *  \brief compute the vector x in  Ax=B
   *  \param SimpleMatrix & B
   */
  void linearSolve(const SiconosMatrix & , SiconosMatrix &);

  /** \fn  PLUFactorizationInPlace(void);
   *  \brief Compute the LU factorization with Partial pivoting.  The result is returned in this (InPlace)
   */
  void PLUFactorizationInPlace();

  /** \fn  SiconosMatrix  PLUInverseInPlace(void);
   *  \brief  use the LU factorization with Partial pivoting of the matrix this (or the matrix this itself if it is triangular)
   *       to compute the inverse matrix.  The result is returned in this (InPlace)
   */
  void  PLUInverseInPlace();

  /** \fn SiconosMatrix  PLUForwardBackward(SiconosMatrix &B);
   *  \brief use the LU factorization with Partial pivoting of the matrix this (or the matrix this itself if it is triangular)
   *  to solve A x = B by forward or backward substitution;
   *  \param the RHS matrix b  which contains the result x
   */
  void  PLUForwardBackwardInPlace(SiconosMatrix &B) ;

  /** \fn SiconosVector  PLUForwardBackward(SiconosVector &B);
   *  \brief use the LU factorization with Partial pivoting of the matrix this (or the matrix this itself if it is triangular)
   *  to solve A x = B by forward or backward substitution;
   *  \param the RHS vector b which contains the result x
   */
  void   PLUForwardBackwardInPlace(SiconosVector &B);

  /** \fn SiconosMatrix BlockMatrixAssemble(vector<SiconosMatrix*>);
   *  \brief build a matrix from n matrices
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix BlockMatrixAssemble(const std::vector<SiconosMatrix*>&);

};

#endif // __SimpleMatrix__
