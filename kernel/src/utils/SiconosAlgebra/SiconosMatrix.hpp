/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file SiconosMatrix.hpp
  Interface for matrices handling.
*/

#ifndef __SiconosMatrix__
#define __SiconosMatrix__

#include <stddef.h>                    // for size_t
#include <iosfwd>                      // for ostream
#include "CSparseMatrix.h"             // for CSparseMatrix
#include "SiconosAlgebraTypeDef.hpp"   // for DenseMat, BandedMat, IdentityMat
#include "SiconosFwd.hpp"              // for SiconosMatrix
#include "SiconosSerialization.hpp"    // for ACCEPT_SERIALIZATION
#include "SiconosVisitor.hpp"          // for VIRTUAL_ACCEPT_VISITORS
class BlockVector;


/** Union of DenseMat pointer, TriangMat pointer BandedMat, SparseMat, SymMat, Zero and Identity mat pointers.
 */
union MATRIX_UBLAS_TYPE
{
  DenseMat *Dense;    // num = 1
  TriangMat *Triang;  // num = 2
  SymMat *Sym;        // num = 3
  SparseMat *Sparse;  // num = 4
  BandedMat *Banded;  // num = 5
  ZeroMat *Zero;      // num = 6
  IdentityMat *Identity; // num = 7
  SparseCoordinateMat *SparseCoordinate; // num = 8
};
/** A STL vector of int */
typedef std::vector<int> VInt;
TYPEDEF_SPTR(VInt)


/**
   Abstract class to provide interface for matrices handling
   
   Matrices can be either block or Simple.
   See Derived classes for details.
   
   In Siconos, a "matrix" can be either a SimpleMatrix or a BlockMatrix, ie a container of several pointers to SiconosMatrix
   
   You can find an overview on how to build and use vectors and matrices in siconos users guide .
   
*/
class SiconosMatrix //: public std::enable_shared_from_this<SiconosMatrix>
{
protected:
  ACCEPT_SERIALIZATION(SiconosMatrix);

  /** A number to specify the type of the matrix: (block or ublas-type)
   *  0-> BlockMatrix, 1 -> DenseMat, 2 -> TriangMat, 3 -> SymMat, 4->SparseMat, 5->BandedMat, 6->zeroMat, 7->IdentityMat
   */
  siconos::UBLAS_TYPE _num;

  /** bool _isSymmetric;
   *  Boolean = true if the Matrix is symmetric
   */
  bool _isSymmetric;

  /** bool _isPositiveDefinite;
   *  Boolean = true if the Matrix is positive definite
   */
  bool _isPositiveDefinite;

  /** default constructor */
  SiconosMatrix() {};

  /** basic constructor
   *
   *  \param type unsigned int type-number of the vector
   */
  SiconosMatrix(siconos::UBLAS_TYPE type);

  /** computes y = subA*x (init =true) or += subA * x (init = false), subA being a submatrix of A (all columns, and rows between start and start+sizeY).
   *  If x is a block vector, it call the present function for all blocks.
   *
   *  \param A a pointer to SiconosMatrix
   *  \param startRow an int, sub-block position
   *  \param x a pointer to a SiconosVector
   *  \param y a pointer to a SiconosVector
   *  \param init a bool
   */
  void private_prod(unsigned int startRow, const SiconosVector& x, SiconosVector& y, bool init) const;


  /** computes res = subA*x +res, subA being a submatrix of A (rows from startRow to startRow+sizeY and columns between startCol and startCol+sizeX).
   *  If x is a block vector, it call the present function for all blocks.
   *
   *  \param A a pointer to SiconosMatrix
   *  \param startRow an int, sub-block position
   *  \param startCol an int, sub-block position
   *  \param x a pointer to a SiconosVector
   *  \param res a DenseVect
   */
  void private_addprod(unsigned int startRow, unsigned int startCol, const SiconosVector& x, SiconosVector& res) const;

public:

  /** Destructor. */
  virtual ~SiconosMatrix() {};

  /** true if the matrix is block else false.
   *
   *  \return a bool
   */
  inline bool isBlock(void) const
  {
    if (_num == siconos::BLOCK) return true ;
    else return false;
  }

  /** determines if the matrix has been inversed in place
   *
   *  \return true if the matrix is inversed
   */
  inline virtual bool isPLUInversed() const
  {
    return false;
  };

  /** true if the matrix is symmetric (the flag is just returned)
   *
   *  \return true if the matrix is symmetric
   */
  inline bool isSymmetric() const
  {
    return _isSymmetric;
  }

  /** set the flag _isSymmetric */
  inline void setIsSymmetric(bool b)
  {
    _isSymmetric= b;
  }

  /** true if the matrix is definite positive (the flag is just returned)
   *
   *  \return true if the matrix is
   */
  inline bool isPositiveDefinite() const
  {
    return _isPositiveDefinite;
  }

  /** set the flag _isPositiveDefinite */
  inline void setIsPositiveDefinite(bool b)
  {
    _isPositiveDefinite= b ;
  }

  /** determines if the matrix is symmetric up to a given tolerance
   *
   *  \return true if the matrix is inversed
   */
  virtual bool checkSymmetry(double tol) const =0;

  /** determines if the matrix has been PLU factorized
   *
   *  \return true if the matrix is factorized
   */
  inline virtual bool isPLUFactorized() const
  {
    return false;
  };

  /** determines if the matrix has been PLU factorized in place
   *
   *  \return true if the matrix is factorized
   */
  inline virtual bool isPLUFactorizedInPlace() const
  {
    return false;
  };

  /** determines if the matrix has been Cholesky factorized
   *
   *  \return true if the matrix is factorized
   */
  inline virtual bool isCholeskyFactorized() const
  {
    return false;
  };

  /** determines if the matrix has been QR factorized
   *
   *  \return true if the matrix is factorized
   */
  inline bool isQRFactorized() const
  {
    return false;
  }

  inline virtual SP::VInt ipiv() const
  {
    SP::VInt dummy;
    return dummy;
  }

  /** determines if the matrix has been factorized
   *
   *  \return true if the matrix is factorized
   */
  inline virtual bool isFactorized() const
  {
    return (isPLUFactorized() || isPLUFactorizedInPlace() || isCholeskyFactorized() || isQRFactorized());
  };


  /** get the number of rows or columns of the matrix
   *
   *  \param index 0 for rows, 1 for columns
   *  \return an int
   */
  virtual unsigned int size(unsigned int index) const = 0;

  /** get the attribute num of current matrix
   *
   *  \return an unsigned int.
   */
  inline siconos::UBLAS_TYPE num() const
  {
    return _num;
  };

  /** get the number of block (i=0, row, i=1 col)
   *
   *  \param i unsigned int(i=0, row, i=1 col)
   *  \return an unsigned int. 1 as default for SimpleMatrix.
   */
  inline virtual unsigned int numberOfBlocks(unsigned int i) const
  {
    return 1;
  };

  /** reserved to BlockMatrix - get the index tab for rows
   *
   *  \return a pointer to a standard vector of int
   */
  virtual const SP::Index tabRow() const ;

  /** reserved to BlockMatrix - get the index tab of columns
   *
   *  \return a pointer to a standard vector of int
   */
  virtual const SP::Index tabCol() const ;

  /** get DenseMat matrix
   *
   *  \param row an unsigned int position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int position of the block (column) - Useless for SimpleMatrix
   *  \return a DenseMat
   */
  virtual const DenseMat getDense(unsigned int row = 0, unsigned int col = 0) const =  0;

  /** get TriangMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat
   */
  virtual const TriangMat getTriang(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get SymMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat
   */
  virtual const SymMat getSym(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get BandedMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat
   */
  virtual const BandedMat getBanded(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get SparseMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat
   */
  virtual const SparseMat getSparse(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get SparseCoordinateMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseCoordinateMat
   */
  virtual const SparseCoordinateMat getSparseCoordinate(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get ZeroMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a ZeroMat
   */
  virtual const ZeroMat getZero(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get  getIdentity matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return an IdentityMat
   */
  virtual const IdentityMat getIdentity(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on DenseMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a DenseMat*
   */
  virtual  DenseMat* dense(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on TriangMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat*
   */
  virtual TriangMat* triang(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on SymMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat*
   */
  virtual SymMat* sym(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on BandedMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat*
   */
  virtual BandedMat* banded(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on SparseMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat*
   */
  virtual SparseMat* sparse(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on SparseCoordinateMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseCoordinateMat*
   */
  virtual SparseCoordinateMat* sparseCoordinate(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on ZeroMat matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a ZeroMat*
   */
  virtual ZeroMat* zero_mat(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** get a pointer on Identity matrix
   *
   *  \param row an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param col an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return an IdentityMat*
   */
  virtual IdentityMat* identity(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** return the address of the array of double values of the matrix
   *  ( for block(i,j) if this is a block matrix)
   *
   *  \param row position for the required block
   *  \param col position for the required block
   *  \return double* : the pointer on the double array
   */
  virtual double* getArray(unsigned int row = 0, unsigned int col = 0) const = 0;

  /** sets all the values of the matrix to 0.0
   */
  virtual void zero() = 0;

  /** Initialize the matrix with random values
   */
  virtual void randomize() = 0;

  /** set an identity matrix
   */
  virtual void eye() = 0;

  /** resize the matrix with nbrow rows and nbcol columns, upper and lower are only useful for BandedMatrix .
   *  The existing elements of the matrix are preseved when specified.
   *
   *  \param nbrow
   *  \param nbcol
   *  \param lower,upper for banded matrices
   *  \param preserve
   */
  virtual void resize(unsigned int nbrow, unsigned int nbcol,
                      unsigned int lower = 0, unsigned int upper = 0, bool preserve = true) = 0;

  /** compute the infinite norm of the matrix
   *
   *  \return a double
   */
  virtual double normInf() const = 0;

  /** display data on standard output
   */
  virtual void display() const = 0;

  /** display data on standard output
   */
  virtual void displayExpert(bool brief = true ) const = 0;

  /** put data of the matrix into a std::string
   *
   *  \return std::string
   */
  virtual std::string toString() const = 0;

  // Note: in the following functions, row and col are general;
  // that means that for a SimpleMatrix m, m(i,j) is index (i,j) element but
  // for a BlockMatrix w that contains 2 SiconosMatrix of size 3
  // w(1, 4) corresponds to the element (1,1) of the second matrix.
  /** get or set the element matrix[i,j]
   *
   *  \param i an unsigned int i
   *  \param j an unsigned int j
   *  \return the element matrix[i,j]
   */
  virtual double& operator()(unsigned int i, unsigned int j) = 0;

  /** get or set the element matrix[i,j]
   *
   *  \param i an unsigned int i
   *  \param j an unsigned int j
   *  \return the element matrix[i,j]
   */
  virtual double operator()(unsigned int i, unsigned int j) const = 0;

  /** return the element matrix[i,j]
   *
   *  \param i an unsigned int i
   *  \param j an unsigned int j
   *  \return a double
   */
  virtual double getValue(unsigned int i, unsigned int j) const = 0;

  /** set the element matrix[i,j]
   *
   *  \param i an unsigned int i
   *  \param j an unsigned int j
   *  \param value
   */
  virtual void setValue(unsigned int i, unsigned int j, double value) = 0;

  /** get block at position row-col if BlockMatrix, else if SimpleMatrix return this
   *
   *  \param row unsigned int row
   *  \param col unsigned int col
   *  \return SP::SiconosMatrix
   */
  virtual SP::SiconosMatrix block(unsigned int row = 0, unsigned int col = 0)
  {
    THROW_EXCEPTION("must be implemented");
  };


  /** get block at position row-col if BlockMatrix, else if SimpleMatrix return this
   *
   *  \param row unsigned int row
   *  \param col unsigned int col
   *  \return SPC::SiconosMatrix
   */
  virtual SPC::SiconosMatrix block(unsigned int row = 0, unsigned int col = 0) const
  {
    THROW_EXCEPTION("must be implemented");
  };

  /** get row index of current matrix and save it into vOut
   *
   *  \param index row we want to get
   *  \param[out] vOut SiconosVector that will contain the desired row
   */
  virtual void getRow(unsigned int index, SiconosVector& vOut) const = 0;

  /** get column index of current matrix and save it into vOut
   *
   *  \param index column we want to get
   *  \param[out] vOut SiconosVector that will contain the desired column
   */
  virtual void getCol(unsigned int index, SiconosVector& vOut) const = 0;

  /** set line row of the current matrix with vector v
   *
   *  \param index row we want to set
   *  \param vIn SiconosVector containing the new row
   */
  virtual void setRow(unsigned int index, const SiconosVector& vIn) = 0;

  /** set column col of the current matrix with vector v
   *
   *  \param index column we want to set
   *  \param vIn a SiconosVector containing the new column
   */
  virtual void setCol(unsigned int index, const SiconosVector& vIn) = 0;

  /** transpose in place: x->trans() is x = transpose of x.
   */
  virtual void trans() = 0;

  /** transpose a matrix: x->trans(m) is x = transpose of m.
   *
   *  \param m the matrix to be transposed.
   */
  virtual void trans(const SiconosMatrix& m) = 0;

  /** operator =
   *
   *  \param m the matrix to be copied
   *  \return SiconosMatrix&
   */
  virtual SiconosMatrix& operator  = (const SiconosMatrix& m) = 0;

  /** operator = to a DenseMat
   *
   *  \param m the DenseMat to be copied
   *  \return SiconosMatrix&
   */
  virtual SiconosMatrix& operator  = (const DenseMat& m) = 0;

  /** operator +=
   *
   *  \param m a matrix to add
   *  \return SiconosMatrix&
   */
  virtual SiconosMatrix& operator +=(const SiconosMatrix& m) = 0;

  /** operator -=
   *
   *  \param m a matrix to subtract
   *  \return SiconosMatrix&
   */
  virtual SiconosMatrix& operator -=(const SiconosMatrix& m) = 0;

  virtual void updateNumericsMatrix() =0;

  virtual NumericsMatrix * numericsMatrix() const
  {
    return nullptr;
  };

  
  /** computes a LU factorization of a general M-by-N matrix
   *  with partial pivoting and row interchanges.
   *  The result is returned in this (InPlace).
   *  Based on Blas dgetrf function for dense matrix and
   *  ublas cholesky decomposition for sparse matrix
   *  (work only for a symmetric matrix and very slow because it uses
   *  matric accessor)
   *  use preferably PLUFactorize()
   */
  virtual void PLUFactorizationInPlace() = 0;

  /** computes a factorization of a general M-by-N matrix
   *  The implementation is based on an internal NumericsMatrix
   */
  virtual void Factorize() = 0;

  /** compute inverse of this thanks to LU factorization with partial pivoting.
   *  This method inverts U and then computes inv(A) by solving the system
   *  inv(A)*L = inv(U) for inv(A).
   *  The result is returned in this (InPlace).
   *  Based on Blas dgetri function for dense function
   */
  virtual void  PLUInverseInPlace() = 0;

  /** solves a system of linear equations A * X = B  (A=this)
   *  for a general N-by-N matrix A using the LU factorization computed
   *  by PLUFactorizationInPlace.
   *  Based on Blas dgetrs function for dense matrix.
   *
   *  \param[in,out] B on input the RHS matrix b; on output the result x
   */
  virtual void  PLUForwardBackwardInPlace(SiconosMatrix &B) = 0;

  /** solves a system of linear equations A * X = B  (A=this)
   *  for a general N-by-N matrix A using the LU factorization computed
   *  by PLUFactorize.
   *
   *  \param[in,out] B on input the RHS matrix b; on output the result x
   */
  virtual void  Solve(SiconosMatrix &B) = 0;

  /** solves a system of linear equations A * X = B  (A=this)
   *  for a general N-by-N matrix A using the LU factorization computed
   *  by PLUFactorizationInPlace.
   *  Based on Blas dgetrs function for dense matrix.
   *
   *  \param[in,out] B on input the RHS matrix b; on output the result x
   */
  virtual void   PLUForwardBackwardInPlace(SiconosVector &B) = 0;

  /** solves a system of linear equations A * X = B  (A=this)
   *  for a general N-by-N matrix A using the LU factorization computed
   *  by PLUFactorize.
   *
   *  \param[in,out] B on input the RHS matrix b; on output the result x
   */
  virtual void   Solve(SiconosVector &B) = 0;

  /** set to false all LU indicators. Useful in case of
      assignment for example.
  */
  virtual void resetLU()
  {
    THROW_EXCEPTION(" SiconosMatrix::resetLU not yet implemented for BlockMatrix.");
  };

  /** set to false all factorization indicators. Useful in case of
      assignment for example.
  */
  virtual void resetFactorizationFlags()
  {
    THROW_EXCEPTION(" SiconosMatrix::resetFactorizationFlags not yet implemented for BlockMatrix.");
  };

  /** return the number of non-zero in the matrix
   *
   *  \param tol the tolerance to consider a number zero (not used if the matrix is sparse)
   *  \return the number of non-zeros
   */
  virtual size_t nnz(double tol = 1e-14);

  /** Fill CSparseMatrix compresses column sparse matrix
   *
   *  \param csc the compressed column sparse matrix
   *  \param row_off
   *  \param col_off
   *  \param tol the tolerance under which a number is considered as equal to zero
   *  \return true if function worked.
   *  \warning not clear that it works for an empty csr matrix with row_off =0  and col_off =0
   */
  bool fillCSC(CSparseMatrix* csc, size_t row_off, size_t col_off, double tol = 1e-14);

  /** Fill CSparseMatrix compresses column sparse matrix
   *
   *  \param csc the compressed column sparse matrix
   *  \param tol the tolerance under which a number is considered as equal to zero
   *  \return true if function worked.
   */
  bool fillCSC(CSparseMatrix* csc, double tol = 1e-14);

  bool fromCSC(CSparseMatrix* csc);
  
  /** return the number of non-zero in the matrix
   *
   *  \param csc the compressed column sparse matrix
   *  \param row_off
   *  \param col_off
   *  \param tol the tolerance to consider a number zero (not used if the matrix is sparse)
   *  \return the number of non-zeros
   */
  bool fillTriplet(CSparseMatrix* csc, size_t row_off, size_t col_off, double tol = 1e-14);

  VIRTUAL_ACCEPT_VISITORS(SiconosMatrix);

  /** \defgroup SiconosMatrixFriends

      List of friend functions of the SimpleMatrix class

      Declared in SimpleMatrixFriends.hpp.
      Implemented in SimpleMatrixFriends.cpp.

   */

  // grant access to private_prod
  friend void prod(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y, bool init);
  // grant access to private_addprod
  friend void prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, bool init);

  /** send data of the matrix to an ostream
   *
   *  \param os An output stream
   *  \param sm a SiconosMatrix
   *  \return The same output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const SiconosMatrix& sm);

  /** multiply the current matrix with a scalar
   *
   *  \param m the matrix to operate on
   *  \param s the scalar
   *  \return SiconosMatrix&
   */
  friend SiconosMatrix& operator *=(SiconosMatrix& m, const double& s);

  /** divide the current matrix with a scalar
   *
   *  \param m the matrix to operate on
   *  \param s the scalar
   *  \return SiconosMatrix&
   */
  friend SiconosMatrix& operator /=(SiconosMatrix& m, const double& s);

 
};

#endif
