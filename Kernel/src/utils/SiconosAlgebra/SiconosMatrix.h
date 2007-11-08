/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

/*! \file SiconosMatrix.h
  \brief Interface for matrices handling.

*/

#ifndef __SiconosMatrix__
#define __SiconosMatrix__

#include "SiconosAlgebra.h"
#include "SiconosMatrixException.h"

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
};

/** Abstract class to provide interface for matrices handling
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \author NDIAYE Abdel-Aziz
 *  \date (creation) 07/21/2006
 *  Matrices can be either block or Simple.
 *  See Derived classes for details.
 *
 *  In Siconos, a "matrix" can be either a SimpleMatrix or a BlockMatrix, ie a container of several pointers to SiconosMatrix
 *
 * You can find an overview on how to build and use vectors and matrices in \ref GS_SicAlgebra .
 *
 */
class SiconosMatrix
{
protected:

  /** Number of rows (Warning: total number of scalar elements, not number of blocks) */
  unsigned int dimRow;

  /** Number of columns (Warning: total number of scalar elements, not number of blocks) */
  unsigned int dimCol;

  /** A number to specify the type of the matrix: (block or ublas-type)
   * 0-> BlockMatrix, 1 -> DenseMat, 2 -> TriangMat, 3 -> SymMat, 4->SparseMat, 5->BandedMat, 6->zeroMat, 7->IdentityMat
   */
  unsigned int num;

  /** default constructor */
  SiconosMatrix();

  /** basic constructor
      \param unsigned int, type-number of the vector
  */
  SiconosMatrix(unsigned int);

  /** basic constructor
      \param unsigned int, type-number of the vector
      \param unsigned int, number of rows
      \param unsigned int, number of columns
  */
  SiconosMatrix(unsigned int, unsigned int, unsigned int);

public:

  /** Destructor. */
  virtual ~SiconosMatrix() {};

  /** true if the matrix is block else false.
   * \return a bool.*/
  inline bool isBlock(void) const
  {
    if (num == 0) return true ;
    else return false;
  }

  /** determines if the matrix is square
   *  \return a bool, true if dimRow == dimCol
   */
  inline bool isSquare() const
  {
    return (dimRow == dimCol);
  };

  /** determines if the matrix has been inversed in place
   *  \return true if the matrix is inversed
   */
  inline virtual bool isInversed() const
  {
    return false;
  };

  /** determines if the matrix has been factorized in place
   *  \return true if the matrix is factorized
   */
  inline virtual bool isFactorized() const
  {
    return false;
  };

  /** get the number of rows or columns of the matrix
   *  \param : unsigned int, 0 for rows, 1 for columns
   *  \return an int
   */
  inline unsigned int size(unsigned int index) const
  {
    if (index == 0) return dimRow;
    else return dimCol;
  };

  /** get the attribute num of current matrix
   * \return an unsigned int.
   */
  inline unsigned int getNum() const
  {
    return num;
  };

  /** get the number of block (i=0, row, i=1 col)
   *  \param unsigned int(i=0, row, i=1 col)
   *  \return an unsigned int. 1 as default for SimpleMatrix.
   */
  inline virtual unsigned int getNumberOfBlocks(unsigned int) const
  {
    return 1;
  };

  /** reserved to BlockMatrix - get the index tab for rows
   * \return a pointer to a standard vector of int
   */
  virtual const Index * getTabRowPtr() const ;

  /** reserved to BlockMatrix - get the index tab of columns
   * \return a pointer to a standard vector of int
   */
  virtual const Index * getTabColPtr() const ;

  /** get an iterator pointing at the beginning of the block matrix
  *  \return a BlockIterator1
  */
  virtual BlockIterator1 begin() ;

  /** get an iterator pointing at the end of the block matrix
  *  \return a BlockIterator1
  */
  virtual BlockIterator1 end() ;

  /** get an iterator pointing at the beginning of the block matrix
  *  \return a BlockIterator1
  */
  virtual ConstBlockIterator1 begin() const;

  /** get an iterator pointing at the end of the block matrix
  *  \return a BlockIterator1
  */
  virtual ConstBlockIterator1 end() const;

  /** get DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a DenseMat
   */
  virtual const DenseMat getDense(unsigned int = 0, unsigned int = 0) const =  0;

  /** get TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat
   */
  virtual const TriangMat getTriang(unsigned int = 0, unsigned int = 0) const = 0;

  /** get SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat
   */
  virtual const SymMat getSym(unsigned int = 0, unsigned int = 0) const = 0;

  /** get BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat
   */
  virtual const BandedMat getBanded(unsigned int = 0, unsigned int = 0) const = 0;

  /** get SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat
   */
  virtual const SparseMat getSparse(unsigned int = 0, unsigned int = 0) const = 0;

  /** get ZeroMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a ZeroMat
   */
  virtual const ZeroMat getZero(unsigned int = 0, unsigned int = 0) const = 0;

  /** get  getIdentity matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return an IdentityMat
   */
  virtual const IdentityMat getIdentity(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on DenseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a DenseMat*
   */
  virtual  DenseMat* getDensePtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat*
   */
  virtual TriangMat* getTriangPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat*
   */
  virtual SymMat* getSymPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat*
   */
  virtual BandedMat* getBandedPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat*
   */
  virtual SparseMat* getSparsePtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on ZeroMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a ZeroMat*
   */
  virtual ZeroMat* getZeroPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get a pointer on Identity matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return an IdentityMat*
   */
  virtual IdentityMat* getIdentityPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** return the adress of the array of double values of the matrix ( for block(i,j) if this is a block matrix)
   *  \param: row position for the required block
   *  \param: col position for the required block
   *  \return double* : the pointer on the double array
   */
  virtual double* getArray(unsigned int = 0, unsigned int = 0) const = 0;

  /** sets all the values of the matrix to 0.0
   */
  virtual void zero() = 0;

  /** set an identity matrix
   */
  virtual void eye() = 0;

  /** resize the matrix with nbrow rows and nbcol columns, upper and lower are only useful for BandedMatrix .
   *   The existing elements of the matrix are preseved when specified.
   *  \param 2 unsigned int: number of rows and columns
   *  \param 2 unsigned int: for banded matrices
   */
  virtual void resize(unsigned int, unsigned int, unsigned int = 0, unsigned int = 0, bool = true) = 0;

  /** compute the infinite norm of the matrix
   *  \return a double
   */
  virtual const double normInf() const = 0;

  /** display data on standard output
   */
  virtual void display() const = 0;

  // Note: in the following functions, row and col are general;
  // that means that for a SimpleMatrix m, m(i,j) is index (i,j) element but
  // for a BlockMatrix w that contains 2 SiconosMatrix of size 3
  // w(1, 4) corresponds to the element (1,1) of the second matrix.
  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \return the element matrix[i,j]
   */
  virtual double& operator()(unsigned int , unsigned int) = 0;

  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \return the element matrix[i,j]
   */
  virtual const double operator()(unsigned int , unsigned int) const = 0;

  /** return the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \return a double
   */
  virtual const double getValue(unsigned int, unsigned int) const = 0;

  /** set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \param the value
   */
  virtual void setValue(unsigned int, unsigned int, double) = 0;

  /** get block at position row-col if BlockMatrix, else if SimpleMatrix return this
   *  \param unsigned int row
   *  \param unsigned int col
   */
  virtual SiconosMatrix* getBlockPtr(unsigned int = 0, unsigned int = 0) = 0;

  /** get block at position row-col if BlockMatrix, else if SimpleMatrix return this
   *  \param unsigned int row
   *  \param unsigned int col
   */
  virtual const SiconosMatrix* getBlockPtr(unsigned int = 0, unsigned int = 0) const = 0;

  /** get row index of current matrix and save it into vOut
   *  \param unsigned int: index of required line
   *  \param ref to SiconosVector: in-out parameter
   */
  virtual void getRow(unsigned int, SiconosVector&) const = 0;

  /** get column index of current matrix and save it into vOut
   *  \param unsigned int: index of required column
   *  \param ref to SiconosVector: in-out parameter
   */
  virtual void getCol(unsigned int, SiconosVector&) const = 0;

  /** set line row of the current matrix with vector v
   *  \param an unsigned int and a SiconosVector
   */
  virtual void setRow(unsigned int, const SiconosVector&) = 0;

  /** set column col of the current matrix with vector v
   *  \param an unsigned int and a SiconosVector
   */
  virtual void setCol(unsigned int, const SiconosVector&) = 0;

  /** transpose in place: x->trans() is x = transpose of x.
   */
  virtual void trans() = 0;

  /** transpose a matrix: x->trans(m) is x = transpose of m.
   *  \param a SiconosMatrix: the matrix to be transposed.
   */
  virtual void trans(const SiconosMatrix&) = 0;

  /** operator =
   *  \param SiconosMatrix : the matrix to be copied
   */
  virtual SiconosMatrix& operator  = (const SiconosMatrix&) = 0;

  /** operator = to a DenseMat
   *  \param the ublas-matrix to be copied
   */
  virtual SiconosMatrix& operator  = (const DenseMat&) = 0;

  /** operator +=
   *  \param SiconosMatrix : a matrix to add
   */
  virtual SiconosMatrix& operator +=(const SiconosMatrix&) = 0;

  /** operator -=
   *  \param SiconosMatrix : a matrix to subtract
   */
  virtual SiconosMatrix& operator -=(const SiconosMatrix&) = 0;

  /** multiply the current matrix with a scalar
   *  \param template, double, int ...
   */
  template <class T> SiconosMatrix& operator *=(const T& s)
  {
    if (num == 0) // BlockMatrix
    {
      BlocksMat::iterator1 it;
      BlocksMat::iterator2 it2;
      for (it = begin(); it != end(); ++it)
      {
        for (it2 = it.begin(); it2 != it.end(); ++it2)
          (**it2) *= s;
      }
    }
    else if (num == 1)
      *getDensePtr() *= s;
    else if (num == 2)
      *getTriangPtr() *= s;
    else if (num == 3)
      *getSymPtr() *= s;
    else if (num == 4)
      *getSparsePtr() *= s;
    else if (num == 5)
      *getBandedPtr() *= s;
    else if (num == 6) {} // nothing!
    else //if(num == 7)
      SiconosMatrixException::selfThrow(" SiconosMatrix *= (double) : invalid type of matrix");

    return *this;
  }

  /** divide the current matrix with a scalar
   *  \param template, double, int ...
   */
  template <class T> SiconosMatrix& operator /=(const T& s)
  {
    if (num == 0) // BlockMatrix
    {
      BlocksMat::iterator1 it;
      BlocksMat::iterator2 it2;
      for (it = begin(); it != end(); ++it)
      {
        for (it2 = it.begin(); it2 != it.end(); ++it2)
          (**it2) /= s;
      }
    }
    else if (num == 1)
      *getDensePtr() /= s;
    else if (num == 2)
      *getTriangPtr() /= s;
    else if (num == 3)
      *getSymPtr() /= s;
    else if (num == 4)
      *getSparsePtr() /= s;
    else if (num == 5)
      *getBandedPtr() /= s;
    else if (num == 6) {} // nothing!
    else //if(num == 7)
      SiconosMatrixException::selfThrow(" SiconosMatrix *= (double) : invalid type of matrix");

    return *this;
  }

  /** computes a LU factorization of a general M-by-N matrix using partial pivoting with row interchanges.
   *  The result is returned in this (InPlace). Based on Blas dgetrf function.
   */
  virtual void PLUFactorizationInPlace() = 0;

  /**  compute inverse of this thanks to LU factorization with Partial pivoting. This method inverts U and then computes inv(A) by solving the system
   *  inv(A)*L = inv(U) for inv(A). The result is returned in this (InPlace). Based on Blas dgetri function.
   */
  virtual void  PLUInverseInPlace() = 0;

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
   *   by PLUFactorizationInPlace. Based on Blas dgetrs function.
   *  \param input: the RHS matrix b - output: the result x
   */
  virtual void  PLUForwardBackwardInPlace(SiconosMatrix &B) = 0;

  /** solves a system of linear equations A * X = B  (A=this) with a general N-by-N matrix A using the LU factorization computed
   *   by PLUFactorizationInPlace.  Based on Blas dgetrs function.
   *  \param input: the RHS matrix b - output: the result x
   */
  virtual void   PLUForwardBackwardInPlace(SiconosVector &B) = 0;

  /** set to false all LU indicators. Useful in case of
      assignment for example.
  */
  virtual void resetLU()
  {
    SiconosMatrixException::selfThrow(" SiconosMatrix::resetLU not yet implemented for BlockMatrix.");
  };

  /** Compares two (block) matrices: true if they have the same number of blocks and if
      blocks which are facing each other have the same size;
      always true if one of the two is a SimpleMatrix.
      \param a SiconosMatrix*
      \param a SiconosMatrix*
  */
  friend const bool isComparableTo(const  SiconosMatrix*, const  SiconosMatrix*);

};

#endif
