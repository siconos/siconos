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

/*! \file SimpleMatrix.hpp
 */

#ifndef __SimpleMatrix__
#define __SimpleMatrix__

#include "SimpleMatrixFriends.hpp"

/**Input parameter for copy and transpose constructor.*/
const std::string transpose = "transpose";

/** A STL vector of int */
typedef std::vector<int> VInt;


/**  Matrix (embedded various types of Boost matrices of double)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 07/21/2006
 *
 * SimpleMatrix is used in the platform to store matrices (mathematical object) of double.
 *
 * Possible types: Siconos::DENSE (default), TRIANGULAR, SYMMETRIC, SPARSE, BANDED, ZERO, Siconos::IDENTITY.
 *
 * \todo: review resize function for Banded, Symetric and Triangular. Error in tests.
 *
 * Notes:
 *  - to handle with sparse matrices see http://freenet-homepage.de/guwi17/ublas/matrix_sparse_usage.html#Q2, for operations improvments.
 *  - See SandBox/Algebra for comparison between ublas (direct call) and Siconos perf.
 *  - Different way to compute matrix-vector or matrix-matrix products are proposed (prod, axpy_prod, gem...) based either on ublas or boost numeric bindings.
 *   See SandBox/Algebra/TestOperators to know which the most performant op. on your system.
 *   axpy_prod is only efficient for sparse or for large objects. For small matrices and vectors it is slower.
 *
 *  See also Siconos Algebra's page in Users Guide, \ref UMsiconosAlgebra.
 *
 *
 */
class SimpleMatrix: public SiconosMatrix
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SimpleMatrix);

  /** Union of The Boost Matrices : DenseMat, TriangMat, SymMat ...
      (See SiconosMatrix.h for more details on MATRIX_UBLAS_TYPE);
  */
  MATRIX_UBLAS_TYPE mat;


private:
  /** VInt _ipiv;
   * The pivot indices obtained from DGETRF (PLUFactorizationInPlace)
   */
  SP::VInt _ipiv;

  /** bool _isPLUFactorized;
   *  Boolean = true if the Matrix has been PLU Factorized in place.
   */
  bool _isPLUFactorized;

  /** bool _isQRFactorized;
   *  Boolean = true if the Matrix has been QR Factorized in place.
   */
  bool _isQRFactorized;

  /** bool _isPLUInversed;
   *  Boolean = true if the Matrix has been Inversed in Place
   */
  bool _isPLUInversed;

  /**  computes res = subA*x +res, subA being a submatrix of A (rows from startRow to startRow+sizeY and columns between startCol and startCol+sizeX).
       If x is a block vector, it call the present function for all blocks.
       \param a pointer to SiconosMatrix A
       \param an int, sub-block position (startRow)
       \param an int, sub-block position (startCol)
       \param a pointer to a SiconosVector, x
       \param a DenseVect, res.
  */
  friend void private_addprod(const SiconosMatrix&, unsigned int, unsigned int, const SiconosVector&, SiconosVector&);

  /**  computes res = a*subA*x +res, subA being a submatrix of A (rows from startRow to startRow+sizeY and columns between startCol and startCol+sizeX).
       If x is a block vector, it call the present function for all blocks.
       \param a double, a
       \param a pointer to SiconosMatrix A
       \param an int, sub-block position (startRow)
       \param an int, sub-block position (startCol)
       \param a pointer to a SiconosVector, x
       \param a DenseVect, res.
  */

  /**  computes res = subA*x +res, subA being a submatrix of trans(A) (rows from startRow to startRow+sizeY and columns between startCol and startCol+sizeX).
       If x is a block vector, it call the present function for all blocks.
       \param a pointer to a SiconosVector, x
       \param a pointer to SiconosMatrix A
       \param an int, sub-block position (startRow)
       \param an int, sub-block position (startCol)
       \param a DenseVect, res.
  */
  friend void private_addprod(SPC::SiconosVector, SPC::SiconosMatrix, unsigned int, unsigned int, SP::SiconosVector);

  /**  computes y = subA*x (init =true) or += subA * x (init = false), subA being a submatrix of A (all columns, and rows between start and start+sizeY).
       If x is a block vector, it call the present function for all blocks.
       \param a pointer to SiconosMatrix A
       \param an int, sub-block position (start)
       \param a pointer to a SiconosVector, x
       \param a pointer to a SiconosVector, y
       \param init, a bool
  */
  friend void private_prod(const SiconosMatrix& A, unsigned int startRow, const SiconosVector& x, SiconosVector& y, bool init);
  friend void private_prod(const SiconosMatrix& A, unsigned int, BlockVector&, const SiconosVector&, bool);
  friend void private_prod(SPC::SiconosMatrix, unsigned int, SPC::SiconosVector , SP::BlockVector, bool);
  friend void private_prod(SPC::SiconosMatrix, unsigned int, SPC::BlockVector , SP::BlockVector, bool);


  friend void private_addprod(double, SPC::SiconosMatrix, unsigned int, unsigned int, SPC::SiconosVector, SP::SiconosVector);


  /**  computes y = a*subA*x (init =true) or += a*subA * x (init = false), subA being a submatrix of A (all columns, and rows between start and start+sizeY).
       If x is a block vector, it call the present function for all blocks.
       \param a double, a
       \param a pointer to SiconosMatrix A
       \param an int, sub-block position (start)
       \param a pointer to a SiconosVector, x
       \param a pointer to a SiconosVector, y
       \param init, a bool
  */
  friend void private_prod(double, SPC::SiconosMatrix, unsigned int, SPC::SiconosVector, SP::SiconosVector, bool);

  /**  computes y = subA*x (init =true) or += subA * x (init = false), subA being a submatrix of trans(A) (all columns, and rows between start and start+sizeY).
       If x is a block vector, it call the present function for all blocks.
       \param a pointer to a SiconosVector, x
       \param a pointer to SiconosMatrix A
       \param an int, sub-block position (start)
       \param a pointer to a SiconosVector, y
       \param init, a bool
  */
  friend void private_prod(SPC::SiconosVector, SPC::SiconosMatrix, unsigned int, SP::SiconosVector, bool);
  friend void private_prod(SPC::BlockVector, SPC::SiconosMatrix, unsigned int, SP::SiconosVector, bool);
  friend void private_prod(SPC::BlockVector, SPC::SiconosMatrix, unsigned int, SP::BlockVector, bool);
  friend void private_prod(SPC::SiconosVector, SPC::SiconosMatrix, unsigned int, SP::BlockVector, bool);

public:
  /** Default constructor
   */
  SimpleMatrix();
  SimpleMatrix(int i);

  /** constructor with the type and the dimension of the Boost matrix
   *  \param row number of rows.
   *  \param col number of columns.
   *  \param typ the type of matrix
   *  \param upper if Siconos::UBLAS_TYPE==SPARSE, number of non-zero terms, if Siconos::UBLAS_TYPE == BANDED, number of diags. under the main diagonal
   *  \param lower if Siconos::UBLAS_TYPE == BANDED, number of diags. over the main diagonal
   */
  SimpleMatrix(unsigned int row, unsigned int col, Siconos::UBLAS_TYPE typ = Siconos::DENSE, unsigned int upper = 1, unsigned int lower = 1);

  /** constructor with the the dimensions of the Boost matrix, a default value and the type.
   *  \param row number of rows.
   *  \param col number of columns.
   *  \param double a, so that *this = [a a a ...]
   *  \param typ the type of matrix
   *  \param upper if Siconos::UBLAS_TYPE==SPARSE, number of non-zero terms, if Siconos::UBLAS_TYPE == BANDED, number of diags. under the main diagonal
   *  \param lower if Siconos::UBLAS_TYPE == BANDED, number of diags. over the main diagonal
   */
  SimpleMatrix(unsigned int row, unsigned int col, double inputValue, Siconos::UBLAS_TYPE typ = Siconos::DENSE, unsigned int upper = 1, unsigned int lower = 1);

  /** copy constructor
   *  \param smat the matrix to copy
   */
  SimpleMatrix(const SimpleMatrix& smat);

  /** copy constructor
   *  \param m the matrix to copy
   */
  SimpleMatrix(const SiconosMatrix& m);

  /** constructor with a DenseMat matrix (see SiconosMatrix.h for details)
   *  \param m a DenseMat
   */
  SimpleMatrix(const DenseMat& m);

  /** constructor with a TriangMat matrix (see SiconosMatrix.h for details)
   *  \param m a TriangMat
   */
  SimpleMatrix(const TriangMat& m);

  /** constructor with a SymMat matrix (see SiconosMatrix.h for details)
   *  \param m a SymMat
   */
  SimpleMatrix(const SymMat& m);

  /** constructor with a BandedMat matrix (see SiconosMatrix.h for details)
   *  \param m a BandedMat
   */
  SimpleMatrix(const BandedMat& m);

  /** constructor with a SparseMat matrix (see SiconosMatrix.h for details)
   *  \param m a SparseMat
   */
  SimpleMatrix(const SparseMat& m);

  /** constructor with a ZeroMat matrix (see SiconosMatrix.h for details)
   *  \param m a ZeroMat
   */
  SimpleMatrix(const ZeroMat& m);

  /** constructor with a IdentityMat matrix (see SiconosMatrix.h for details)
   *  \param m a IdentityMat
   */
  SimpleMatrix(const IdentityMat& m);

  /** constructor with an input file
   *  \param file the input file path
   *  \param ascii a boolean to indicate if the file is in ascii
   */
  SimpleMatrix(const std::string& file, bool ascii = true);

  /** destructor
   */
  ~SimpleMatrix();
  //************************** GETTERS/SETTERS  **************************

  /** determines if the matrix has been inversed
   *  \return true if the matrix is inversed
   */
  inline bool isPLUInversed() const
  {
    return _isPLUInversed;
  }

  /** determines if the matrix has been factorized
   *  \return true if the matrix is factorized
   */
  inline bool isPLUFactorized() const
  {
    return _isPLUFactorized;
  }

  /** determines if the matrix has been factorized
   *  \return true if the matrix is factorized
   */
  inline bool isQRFactorized() const
  {
    return _isQRFactorized;
  }

  inline SP::VInt ipiv() const
  {
    return _ipiv;
  }

  /** get DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a DenseMat
   */
  const DenseMat getDense(unsigned int = 0, unsigned int = 0) const;

  /** get TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat
   */
  const TriangMat getTriang(unsigned int = 0, unsigned int = 0) const;

  /** get SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat
   */
  const SymMat getSym(unsigned int = 0, unsigned int = 0) const;

  /** get BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat
   */
  const BandedMat getBanded(unsigned int = 0, unsigned int = 0) const;

  /** get SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat
   */
  const SparseMat getSparse(unsigned int = 0, unsigned int = 0) const;

  /** get ZeroMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a ZeroMat
   */
  const ZeroMat getZero(unsigned int = 0, unsigned int = 0) const;

  /** get  getIdentity matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return an IdentityMat
   */
  const IdentityMat getIdentity(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a DenseMat*
   */
  DenseMat* dense(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat*
   */
  TriangMat* triang(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat*
   */
  SymMat* sym(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat*
   */
  BandedMat* banded(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat*
   */
  SparseMat* sparse(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on ZeroMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a ZeroMat*
   */
  ZeroMat* zero(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on Identity matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return an IdentityMat*
   */
  IdentityMat* identity(unsigned int = 0, unsigned int = 0) const;

  /** return the adress of the array of double values of the matrix
   *  \param: row position for the required block ->useless for SimpleMatrix
   *  \param: col position for the required block ->useless for SimpleMatrix
   *  \return double* : the pointer on the double array
   */
  double* getArray(unsigned int = 0, unsigned int = 0) const;

  /** sets all the values of the matrix to 0.0
   */
  void zero();

  /** Initialize the matrix with random values
   */
  void randomize();

  /** Initialize a symmetric matrix with random values
   */
  void randomize_sym();

  /** set an identity matrix
   */
  void eye();

  /** copy the matrix data to the array given in parameter'
   * Works only for dense matrices !
   * \param data array where the matrix is copied
   * \return the size of the matrix
   */
  unsigned copyData(double* data) const;

  /** get the number of rows or columns of the matrix
   *  \param index 0 for rows, 1 for columns
   *  \return the size
   */
  unsigned int size(unsigned int index) const;

  /** resize the matrix with nbrow rows and nbcol columns The existing elements of the matrix are preseved when specified.
   *  \exception SiconosMatrixException
   */
  void resize(unsigned int, unsigned int, unsigned int = 0, unsigned int = 0, bool = true);

  /** compute the infinite norm of the matrix
   *  \return a double
   */
  double normInf() const;

  /** Compute the normInf for each column*/
  void normInfByColumn(SP::SiconosVector vIn) const;

  /** compute the determinant of the matrix (use LU factorization)
   *  \return a double
   */
  double det() const;

  /** display data on standard output
   */
  void display() const;


  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double& operator()(unsigned int , unsigned int);

  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double operator()(unsigned int , unsigned int) const;

  /** return the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \return a double
   */
  double getValue(unsigned int, unsigned int) const;

  /** set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \param the value
   */
  void setValue(unsigned int, unsigned int, double);

  /** set block starting at row "posRow" (first argument) and column "posCol" (second argument) with m (last arg)
   *  \param posRow an int, row-position of the first element of the required block
   *  \param posCol an int, col-position of the first element of the required block
   *  \param a ref to a SiconosMatrix  m
   */
  void setBlock(unsigned int, unsigned int, const SiconosMatrix&);

  friend void setBlock(SPC::SiconosMatrix , SP::SiconosMatrix , const Index&, const Index&);

  /** get block at position row-col, ie return this in SimpleMatrix case
   *  \param unsigned int row
   *  \param unsigned int col
   */
  inline SP::SiconosMatrix block(unsigned int = 0, unsigned int = 0)
  {
    return shared_from_this();
  };

  /** get block at position row-col, ie return this in SimpleMatrix case
   *  \param unsigned int row
   *  \param unsigned int col
   */
  inline SPC::SiconosMatrix block(unsigned int = 0, unsigned int = 0) const
  {
    return shared_from_this();
  };

  /** get row index of current matrix and save it into vOut
   *  \param index row we want to get
   *  \param[out] vOut SiconosVector that will contain the desired row
   */
  void getRow(unsigned int index, SiconosVector& vOut) const;

  /** get column index of current matrix and save it into vOut
   *  \param index column we want to get
   *  \param[out] vOut SiconosVector that will contain the desired column
   */
  void getCol(unsigned int index, SiconosVector& vOut) const;

  /** set line row of the current matrix with vector v
   *  \param index row we want to set
   *  \param vIn SiconosVector containing the new row
   */
  void setRow(unsigned int index, const SiconosVector& vIn);

  /** set column col of the current matrix with vector v
   *  \param index column we want to set
   *  \param vIn a SiconosVector containing the new column
   */
  void setCol(unsigned int index, const SiconosVector& vIn);

  /** get column number index of current matrix, starting from element at position pos and save it into vOut
   *  \param index index of required column
   *  \param pos index of the first required element in the column
   *  \param[out] vOut a SP::SiconosVector
   */
  void getSubCol(unsigned int index, unsigned int pos, SP::SiconosVector vOut) const;

  /** get row number index of current matrix, starting from element at position pos and save it into vOut
   *  \param index index of the required row
   *  \param pos index of the first required element in the row
   *  \param[out] vOut a SP::SiconosVector that will contain the sub row
   */
  void getSubRow(unsigned int index, unsigned int pos, SP::SiconosVector vOut) const;

  /** set column number index of current matrix, starting from element at position pos, with vIn
   *  \param index index of required column
   *  \param pos index of the first required element in the column
   *  \param vIn a vector
   */
  void setSubCol(unsigned int index, unsigned int pos, SP::SiconosVector vIn);

  /** set row number index of current matrix, starting from element at position pos, with vIn
   *  \param index index of required row
   *  \param pos index of the first required element in the row
   *  \param vIn a vector
   */
  void setSubRow(unsigned int index, unsigned int pos, SP::SiconosVector vIn);

  /** add the input matrix to the elements starting from position i (row) and j (col).
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \param a SiconosMatrix
   */
  void addBlock(unsigned int, unsigned int, const SiconosMatrix&);

  /** subtract the input matrix to the elements starting from position i (row) and j (col).
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \param a SiconosMatrix
   */
  void subBlock(unsigned int, unsigned int, const SiconosMatrix&);

  /** transpose in place: x->trans() is x = transpose of x.
   */
  void trans();

  /** transpose a matrix: x->trans(m) is x = transpose of m.
   *  \param mat the matrix to be transposed.
   */
  void trans(const SiconosMatrix& mat);

  /** assignment
   *  \param m the matrix to be copied
   */
  SimpleMatrix& operator = (const SiconosMatrix& m);

  /** assignment
   *  \param m the matrix to be copied
   */
  SimpleMatrix& operator = (const SimpleMatrix& m);

  /** assignment to a DenseMat
   *  \param the matrix to be copied
   */
  SimpleMatrix& operator = (const DenseMat& m);

  /** operator +=
   *  \param m a matrix to add
   */
  SimpleMatrix& operator +=(const SiconosMatrix& m);

  /** operator -=
   *  \param m a matrix to subtract
   */
  SimpleMatrix& operator -=(const SiconosMatrix& m);

  /** computes an LU factorization of a general M-by-N matrix using partial pivoting with row interchanges.
   *  The result is returned in this (InPlace). Based on Blas dgetrf function.
   */
  void PLUFactorizationInPlace();

  /**  compute inverse of this thanks to LU factorization with Partial pivoting.
   * This method inverts U and then computes inv(A) by solving the system
   *  inv(A)*L = inv(U) for inv(A). The result is returned in this (InPlace). Based on Blas dgetri function.
   */
  void PLUInverseInPlace();

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
   *   by PLUFactorizationInPlace. Based on Blas dgetrs function.
   *  \param[in,out] B on input the RHS matrix b; on output the result x
   */
  void PLUForwardBackwardInPlace(SiconosMatrix& B);

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
   *   by PLUFactorizationInPlace.  Based on Blas dgetrs function.
   *  \param[in,out] B on input the RHS matrix b; on output the result x
   */
  void PLUForwardBackwardInPlace(SiconosVector& B);

  /** solves a system of linear equations A * X = B  (A=this)
      with a general N-by-N matrix A using the Least squares method
   *  \param[in,out] B on input the RHS matrix b; on output the result x
   */
  void SolveByLeastSquares(SiconosMatrix& B);

  /** solves a system of linear equations A * X = B  (A=this)
       with a general N-by-N matrix A using the Least squares method
   *  \param[in,out] B on input the RHS matrix b; on output the result x
    */
  void SolveByLeastSquares(SiconosVector& B);

  /** set to false all LU indicators. Useful in case of
      assignment for example.
  */

  void resetLU();

  /** set to false all QR indicators. Useful in case of
      assignment for example.
  */
  void resetQR();

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();
  
  friend const SimpleMatrix operator * (const SiconosMatrix&, double);
  //  friend const SimpleMatrix operator * (const SimpleMatrix&,const SimpleMatrix&);

  friend  SP::SimpleMatrix operator * (const SP::SimpleMatrix, const SP::SimpleMatrix);

  friend  void operator +=(SP::SiconosMatrix, SP::SimpleMatrix);

  friend  SimpleMatrix operator * (double , const SiconosMatrix&);

  friend const SimpleMatrix operator /(const SiconosMatrix&, double);

  friend const SimpleMatrix operator +(const SiconosMatrix&, const SiconosMatrix&);

  friend SP::SimpleMatrix operator +(const SP::SimpleMatrix, const SP::SimpleMatrix);

  friend void add(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&);

  friend const SimpleMatrix operator -(const SiconosMatrix&, const SiconosMatrix&);

  friend void sub(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&);

  friend bool operator == (const SiconosMatrix&, const SiconosMatrix&);

  friend bool operator!= (const SiconosMatrix&, const SiconosMatrix&);

  friend const SimpleMatrix pow(const SimpleMatrix&, unsigned int);

  friend const SimpleMatrix prod(const SiconosMatrix&, const SiconosMatrix&);

  friend void prod(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&, bool);

  friend void axpy_prod(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&, bool);

  friend const SiconosVector prod(const SiconosMatrix&, const SiconosVector&);

  friend void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init);

  friend void prod(const SiconosMatrix&, const BlockVector&, SiconosVector&, bool);

  friend void prod(const SiconosMatrix&, const SiconosVector&, BlockVector&, bool);

  friend void prod(double, const SiconosMatrix&, const SiconosVector&, SiconosVector&, bool);

  friend void prod(const SiconosVector&, const SiconosMatrix&, SiconosVector&, bool);

  friend void prod(const SiconosVector&, const SiconosMatrix&, BlockVector&, bool);

  friend void subprod(const SiconosMatrix&, const SiconosVector&, SiconosVector&, const Index&, bool);

  friend void axpy_prod(const SiconosMatrix&, const SiconosVector&, SiconosVector&, bool);

  friend void gemvtranspose(double, const SiconosMatrix&, const SiconosVector&, double, SiconosVector&);

  friend void gemv(double, const SiconosMatrix&, const SiconosVector&, double, SiconosVector&);

  friend void gemmtranspose(double, const SiconosMatrix&, const SiconosMatrix&, double, SiconosMatrix&);

  friend void gemm(double, const SiconosMatrix&, const SiconosMatrix&, double, SiconosMatrix&);

  friend void scal(double, const SiconosMatrix&, SiconosMatrix&, bool);

  friend void invertMatrix(const SimpleMatrix&, SimpleMatrix&);

  /** returns a vector of maximum relative error for each column
   * \param data the matrix filled with simulation results
   * \param data the matrix filled with the reference values
   * \return a pointer filled with the maximum relative error for each value in data
   */
  friend SP::SiconosVector compareMatrices(const SimpleMatrix& data, const SimpleMatrix& ref);

};



#endif

