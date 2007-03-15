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

/*! \file SimpleMatrix.h
 */

#ifndef __SimpleMatrix__
#define __SimpleMatrix__

#include "SiconosMatrix.h"

/**Input parameter for copy and transpose constructor.*/
const std::string transpose = "transpose";

/**  Matrix (embedded various types of Boost matrices of double)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.0.1.
 *   \date (Creation) 07/21/2006
 *
 *
 * SimpleMatrix is used in the platform to store matrices (mathematical object).
 *
 * \todo: review resize function for Banded, Symetric and Triangular. Error in tests.
 *
 */
class SimpleMatrix: public SiconosMatrix
{
private:

  /**  unsigned int num
   * an unsigned int which make a correspondance with Boost Matrices: 1 -> DenseMat, 2 -> TriangMat, 3 -> SymMat, 4->SparseMat, 5->BandedMat, 6->zeroMat, 7->IdentityMat
   */
  unsigned int num;

  /** Mat mat (See SiconosMatrix.h for more details on Mat type);
   * union of The Boost Matrices : DenseMat, TriangMat, SymMat are encapsulated.
   */
  Mat mat;

  /** constructor with the type of the Boost matrix
   *  \param TYP
   */
  SimpleMatrix(TYP = DENSE);

  /** std::vector<int> ipiv;
   * The pivot indices obtained from DGETRF (PLUFactorizationInPlace)
   */
  std::vector<int> ipiv;

  /** bool isPLUFactorized;
   *  Boolean = true if the Matrix is PLU Factorized
   */
  bool isPLUFactorized;

  /** bool isPLUInversed;
   *  Boolean = true if the Matrix is Inversed in Place
   */
  bool isPLUInversed;

public:
  /***************************** CONSTRUCTORS ****************************/

  /** copy constructor
   *  \param SimpleMatrix
   */
  SimpleMatrix(const SimpleMatrix&);

  /** copy constructor
   *  \param SiconosMatrix
   */
  SimpleMatrix(const SiconosMatrix&);

  /** copy and transpose constructor
   *  \param a string: "trans"
   *  \param SiconosMatrix
   */
  SimpleMatrix(const std::string&, const SiconosMatrix&);

  /** constructor with the type and the dimension of the Boost matrix
   *  \param 2 unsigned int (number of rows and columns)
   *  \param TYP
   *  \param unsigned int, if TYP==SPARSE, number of non-zero terms, if TYP == BANDED, number of diags. under the main diagonal
   *  \param unsigned int, if TYP == BANDED, number of diags. over the main diagonal
   */
  SimpleMatrix(unsigned int, unsigned int, TYP = DENSE, unsigned int = 1, unsigned int = 1);

  /** constructor with a DenseMat matrix (see SiconosMatrix.h for details)
   *  \param a DenseMat
   */
  SimpleMatrix(const DenseMat&);

  /** constructor with a TriangMat matrix (see SiconosMatrix.h for details)
   *  \param a TriangMat
   */
  SimpleMatrix(const TriangMat&);

  /** constructor with a SymMat matrix (see SiconosMatrix.h for details)
   *  \param a SymMat
   */
  SimpleMatrix(const SymMat&);

  /** constructor with a BandedMat matrix (see SiconosMatrix.h for details)
   *  \param a BandedMat
   */
  SimpleMatrix(const BandedMat&);

  /** constructor with a SparseMat matrix (see SiconosMatrix.h for details)
   *  \param a SparseMat
   */
  SimpleMatrix(const SparseMat&);

  /** constructor with a ZeroMat matrix (see SiconosMatrix.h for details)
   *  \param a ZeroMat
   */
  SimpleMatrix(const ZeroMat&);

  /** constructor with a IdentityMat matrix (see SiconosMatrix.h for details)
   *  \param a IdentityMat
   */
  SimpleMatrix(const IdentityMat&);

  /** constructor with the type of the boost matrix, a vector of the values and the dimensions
   *  of the matrix, the integers upper and lower are useful only for BandedMat
   *  \param a TYP
   *  \param a std::vector<double>
   *  \param unsigned int
   */
  SimpleMatrix(const std::vector<double>& , unsigned int, unsigned int = 0, TYP = DENSE, unsigned int = 0, unsigned int = 0);

  /** constructor with an input file
   *  \param a std::string which contain the file path
   *  \param a boolean to indicate if the file is in ascii
   */
  SimpleMatrix(const std::string&, bool = true);

  /** destructor
   */
  ~SimpleMatrix(void);
  //************************** GETTERS/SETTERS  **************************

  /** determines if the matrix is square
   *  \return true if the matrix is square
   */
  inline bool isSquare() const
  {
    return (size(0) == size(1));
  } ;

  /** determines if the matrix has been inversed
   *  \return true if the matrix is inversed
   */
  inline bool isInversed() const
  {
    return false;
  }

  /** determines if the matrix has been factorized
   *  \return true if the matrix is factorized
   */
  inline bool isFactorized() const
  {
    return false;
  }

  /** get the attribute num of current matrix, useless for Block Matrix
   * \return an unsigned int.
   */
  inline unsigned int getNum(void) const
  {
    return num;
  };

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
  DenseMat* getDensePtr(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat*
   */
  TriangMat* getTriangPtr(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat*
   */
  SymMat* getSymPtr(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat*
   */
  BandedMat* getBandedPtr(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat*
   */
  SparseMat* getSparsePtr(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on ZeroMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a ZeroMat*
  */
  ZeroMat* getZeroPtr(unsigned int = 0, unsigned int = 0) const;

  /** get a pointer on Identity matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return an IdentityMat*
  */
  IdentityMat* getIdentityPtr(unsigned int = 0, unsigned int = 0) const;

  /** get BlocksMat matrix
   *  \useless for SimpleMatrix
   *  \return a BlocksMat
   */
  const BlocksMat getAllBlocks(void) const;

  /** get block corresponding to lines given in numRow and columns in numCol
   *  \param 2 unsigned int for indexes and a SiconosMatrix (in-out paramater)
   */
  void getBlock(unsigned int, unsigned int, SiconosMatrix&) const;

  /** get block at position row-col, ie return this in SimpleMatrix case
   *  \param unsigned int row
   *  \param unsigned int col
   */
  inline SiconosMatrix* getBlockPtr(unsigned int = 0, unsigned int = 0)
  {
    return this;
  };

  /** get std::deque of bool
   *   \useless for SimpleMatrix
   *   \return a std::deque<bool>
   */
  const std::deque<bool> getBlockAllocated(void) const;

  /** get row index of current matrix and save it into vOut
   *  \param unsigned int: index of required line
   *  \param ref to SimpleVector: in-out parameter
   */
  void getRow(unsigned int, SimpleVector&) const;

  /** get column index of current matrix and save it into vOut
   *  \param unsigned int: index of required column
   *  \param ref to SimpleVector: in-out parameter
   */
  void getCol(unsigned int, SimpleVector&) const;

  /** set line row of the current matrix with vector v
   *  \param an unsigned int and a SimpleVector
   */
  void setRow(unsigned int, const SimpleVector&);

  /** set column col of the current matrix with vector v
   *  \param an unsigned int and a SimpleVector
   */
  void setCol(unsigned int, const SimpleVector&);

  /** return the adress of the array of double values of the matrix
   *  \param: row position for the required block ->useless for SimpleMatrix
   *  \param: col position for the required block ->useless for SimpleMatrix
   *  \return double* : the pointer on the double array
   */
  double* getArray(unsigned int = 0, unsigned int = 0) const;

  /** sets all the values of the matrix to 0.0
   */
  void zero(void);

  /** set an identity matrix
   */
  void eye(void);

  /** Computes dim according to the matrix type.
   */
  void computeDim();

  /** resize the matrix with nbrow rows and nbcol columns The existing elements of the matrix are preseved when specified.
   *  \exception SiconosMatrixException
   */
  void resize(unsigned int, unsigned int, unsigned int = 0, unsigned int = 0, bool = true);

  /** compute the infinite norm of the matrix
   *  \return a double
   */
  const double normInf(void) const;

  /** transpose in place: x->trans() is x = transpose of x.
   */
  void trans();

  /** transpose a matrix: x->trans(m) is x = transpose of m.
   *  \param a SiconosMatrix: the matrix to be transposed.
   */
  void trans(const SiconosMatrix&);

  /** display data on standard output
   */
  void display(void) const;

  // --- MATRICES HANDLING AND OPERATORS ---

  /** copy the matrix "blockMat" into the matrix "mat" at the position (xPos, yPos)
   *  \param SiconosMatrix& : the block matrix to be copied in the current matrix
   *  \param unsigned int : the line position to start the copy of the blockmatrix
   *  \param unsigned int : the column position to start the copy of the blockmatrix
   */
  void matrixCopy(const SiconosMatrix&, unsigned int, unsigned int);

  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double& operator()(unsigned int , unsigned int);

  /** return the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \return a double
   */
  double getValue(unsigned int, unsigned int);

  /** set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \param the value
   */
  void setValue(unsigned int, unsigned int, double);

  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double operator()(unsigned int , unsigned int) const;

  /** operator =
   *  \param SiconosMatrix : the matrix to be copied
   */
  SimpleMatrix& operator = (const SiconosMatrix&);

  /** operator =
   *  \param SimpleMatrix : the matrix to be copied
   */
  SimpleMatrix& operator = (const SimpleMatrix&);

  /** operator /=
   *  \param double, a scalar
   */
  SimpleMatrix& operator /= (double);

  /** operator /=
   *  \param int, a scalar
   */
  SimpleMatrix& operator /= (int);

  /** operator +=
   *  \param SiconosMatrix : a matrix to add
   */
  SimpleMatrix& operator +=(const SiconosMatrix&);

  /** operator -=
   *  \param SiconosMatrix : a matrix to subtract
   */
  SimpleMatrix& operator -=(const SiconosMatrix&);

  /** operator *=
   *  \param double, a scalar
   */
  SimpleMatrix& operator *= (double);

  /** operator *=
   *  \param int, a scalar
   */
  SimpleMatrix& operator *= (int);

  /** computes an LU factorization of a general M-by-N matrix using partial pivoting with row interchanges.
   *  The result is returned in this (InPlace). Based on Blas dgetrf function.
   */
  void PLUFactorizationInPlace(void);

  /**  compute inverse of this thanks to LU factorization with Partial pivoting. This method inverts U and then computes inv(A) by solving the system
   *  inv(A)*L = inv(U) for inv(A). The result is returned in this (InPlace). Based on Blas dgetri function.
   */
  void PLUInverseInPlace(void);

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
   *   by PLUFactorizationInPlace. Based on Blas dgetrs function.
   *  \param input: the RHS matrix b - output: the result x
   */
  void PLUForwardBackwardInPlace(SiconosMatrix &B);

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
   *   by PLUFactorizationInPlace.  Based on Blas dgetrs function.
   *  \param input: the RHS matrix b - output: the result x
   */
  void PLUForwardBackwardInPlace(SiconosVector &B);

  /** set to false all LU indicators. Useful in case of
      assignment for example.
  */
  void resetLU();

  /**: A==B when (A-B).normInf()<tolerance
   * \param 2 SiconosMatrix
   * \return a boolean
   */
  friend bool operator == (const SiconosMatrix&, const SiconosMatrix&);

  /** Addition of two matrices
   *  \param 2 SiconosMatrix
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   *  \exception SiconosMatrixException, if the two matrices have different types, in this case use function add
   */
  friend SimpleMatrix operator +(const SiconosMatrix&, const SiconosMatrix&);

  /** subtraction of two matrices
   *  \param 2 SiconosMatrix
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   *  \exception SiconosMatrixException, if the two matrices have different types, in this case use function sub
   */
  friend SimpleMatrix operator - (const SiconosMatrix&, const SiconosMatrix&);

  /** multiplication of two matrices
   *  \param 2 SiconosMatrix
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the two matrices have different types, in this case use function prod
   */
  friend SimpleMatrix operator *(const SiconosMatrix&, const SiconosMatrix&);

  /** multiplication of a matrix by a double
   *  \param a SiconosMatrix
   *  \param a double
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix operator * (const SiconosMatrix&, double);


  /** multiplication of a matrix by an int
   *  \param a SiconosMatrix
   *  \param an int
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix operator *(const SiconosMatrix&, int);

  /** multiplication of a matrix by a double
   *  \param a double
   *  \param a SiconosMatrix
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix operator * (double , const SiconosMatrix&);


  /** multiplication of a matrix by an int
   *  \param an int
   *  \param a SiconosMatrix
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix operator *(int, const SiconosMatrix&);

  /** division of the matrix by a double
   *  \param a SiconosMatrix
   *  \param a double
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the double d = 0
   */
  friend SimpleMatrix operator /(const SiconosMatrix&, double);

  /** division of the matrix by an int
   *  \param a SiconosMatrix
   *  \param an int
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the int d = 0
   */
  friend SimpleMatrix operator / (const SiconosMatrix&, int);

  /** Addition of two matrices
   *  \param 2 SiconosMatrix
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   */
  friend SimpleMatrix add(const SiconosMatrix&, const SiconosMatrix&);

  /** subtraction of two matrices
   *  \param 2 SiconosMatrix
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   */
  friend SimpleMatrix sub(const SiconosMatrix&, const SiconosMatrix&);

  /** multiplication of two matrices
   *  \param 2 SiconosMatrix
   *  \return a SimpleMatrix
   */
  friend SimpleMatrix prod(const SiconosMatrix&, const SiconosMatrix&);

  /** compute the power of the matrix (!)
   *  \return a SimpleMatrix
   *  \exception SiconosMatrixException, if the power < 0
   */
  friend SimpleMatrix pow(const SimpleMatrix&, unsigned int);

  /** compute the product matrix-vector
   *  \return a SimpleVector
   */
  friend SimpleVector prod(const SiconosMatrix&, const SiconosVector&);

  /** gemv(transA, a, A, x, b, y) computes y = a*op(A)*x + b*y
      with transX = CblasNoTrans (op(X) = X), CblasTrans (op(X) = transpose(X)), CblasConjTrans (op(X) = conj(X))
      This function encapsulates atlas::gemv.
      \param CBLAS_TRANSPOSE, op for A
      \param a double, a (in)
      \param a SiconosMatrix, A (in)
      \param a SiconosVector, x (in)
      \param a double, b (in)
      \param a SiconosVector, y (in-out)
  */
  friend void gemv(CBLAS_TRANSPOSE, double, const SiconosMatrix&, const SiconosVector&, double, SiconosVector&);

  /** gemv(a, A, x, b, y) computes y = a*A*x+ b*y
      This function encapsulates atlas::gemm.
      \param a double, a (in)
      \param a SiconosMatrix, A (in)
      \param a SiconosVector, x (in)
      \param a double, b (in)
      \param a SiconosVector, y (in-out)
  */
  friend void gemv(double, const SiconosMatrix&, const SiconosVector&, double, SiconosVector&);

  /** prod(A, x, y) computes y = A*x
      \param a SiconosMatrix, A (in)
      \param a SiconosVector, x (in)
      \param a SiconosVector, y (in-out)
  */
  friend void prod(const SiconosMatrix&, const SiconosVector&, SiconosVector&);

  /** gemm(transA, transB, a, A, B, b, C) computes C = a*op(A)*op(B) + b*C
      with transX = CblasNoTrans (op(X) = X), CblasTrans (op(X) = transpose(X)), CblasConjTrans (op(X) = conj(X))
      This function encapsulates atlas::gemm.
      \param CBLAS_TRANSPOSE, op for A
      \param CBLAS_TRANSPOSE, op for B
      \param a double, a (in)
      \param a SiconosMatrix, A (in)
      \param a SiconosMatrix, B (in)
      \param a double, b (in)
      \param a SiconosMatrix, C (in-out)
  */
  friend void gemm(const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE, double, const SiconosMatrix&, const SiconosMatrix&, double, SiconosMatrix&);

  /** gemm(a, A, B, b, C) computes C = a*A*B+ b*C
      This function encapsulates atlas::gemm.
      \param a double, a (in)
      \param a SiconosMatrix, A (in)
      \param a SiconosMatrix, B (in)
      \param a double, b (in)
      \param a SiconosMatrix, C (in-out)
  */
  friend void gemm(double, const SiconosMatrix&, const SiconosMatrix&, double, SiconosMatrix&);

  /** prod(A, B, C) computes C = A*B
      This function encapsulates atlas::gemm
      \param a SiconosMatrix, A (in)
      \param a SiconosMatrix, B (in)
      \param a SiconosMatrix, C (in-out)
  */
  friend void prod(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&);

};
#endif
