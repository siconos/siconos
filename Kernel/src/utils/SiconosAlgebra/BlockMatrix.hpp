/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file BlockMatrix.hpp
  \brief Object to handle block-matrices.

*/

#ifndef __BlockMatrix__
#define __BlockMatrix__

#include "SiconosMatrix.hpp"

/** Object to handle block-matrices (ie lists of SiconosMatrix*)
 *
 * \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 21/07/2006
 *
 *  A BlockMatrix is a boost::ublas::compressed_matrix of SiconosMatrix*.
 * The blocks positions are given by two Index objects, tabRow and tabCol.
 *
 * If block 1 is n1xm1, block2 n2xm2, block3 n3xm3 ..., then:\n
 *  tabRow = [ n1 n1+n2 n1+n2+n3 ...] \n
 *  tabCol = [ m1 m1+m2 m1+m2+m3 ...] \n
 *
 */
class BlockMatrix : public SiconosMatrix
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BlockMatrix);


  /** A container of pointers to SiconosMatrix
   */
  SP::BlocksMat _mat;

  /** list of blocks dimension - tabRow[i] = tabRow[i-1] + ni, ni being the number of rows of block i.
   */
  SP::Index _tabRow;

  /** list of blocks dimension - tabCol[i] = tabCol[i-1] + ni, ni being the number of columns of block i.
   */
  SP::Index _tabCol;

  /** default constructor
   */
  BlockMatrix() {};

public:

  /** copy constructor
   *  \param SiconosMatrix
   */
  BlockMatrix(const SiconosMatrix&);

  /** copy constructor
   *  \param BlockMatrix
   */
  BlockMatrix(const BlockMatrix&);

  /** constructor with a list of pointer to SiconosMatrix (!links with pointer, no copy!)
   *  \param vector<SiconosMatrix*>
   *  \param unsigned int: number of blocks in a row
   *  \param unsigned int: number of col in a row
   */
  BlockMatrix(const std::vector<SP::SiconosMatrix>&, unsigned int, unsigned int);

  /** contructor with a list of 4 pointer to SiconosMatrix (!links with pointer, no copy!)
   *  \param SP::SiconosMatrix  m1, block (0,0)
   *  \param SP::SiconosMatrix  m2, block (0,1)
   *  \param SP::SiconosMatrix  m3, block (1,0)
   *  \param SP::SiconosMatrix  m4, block (1,1)
   */
  BlockMatrix(SP::SiconosMatrix, SP::SiconosMatrix, SP::SiconosMatrix, SP::SiconosMatrix);

  /** destructor
   */
  ~BlockMatrix(void);

  /** get the number of block (i=0, row, i=1 col)
   *  \param unsigned int(i=0, row, i=1 col)
   *  \return an unsigned int
   */
  unsigned int getNumberOfBlocks(unsigned int) const;

  /** get an iterator pointing at the beginning of the block matrix
   *  \return a BlockIterator1
   */
  BlockIterator1 begin() ;

  /** get an iterator pointing at the end of the block matrix
   *  \return a BlockIterator1
   */
  BlockIterator1 end() ;

  /** get an iterator pointing at the beginning of the block matrix
   *  \return a BlockIterator1
   */
  ConstBlockIterator1 begin() const;

  /** get an iterator pointing at the end of the block matrix
   *  \return a BlockIterator1
   */
  ConstBlockIterator1 end() const;

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
  const SymMat getSym(unsigned int = 0, unsigned int = 0)const;

  /** get BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat
   */
  const BandedMat getBanded(unsigned int = 0, unsigned int = 0)const;

  /** get SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat
   */
  const SparseMat getSparse(unsigned int = 0, unsigned int = 0)const;

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
  DenseMat* dense(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on TriangMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a TriangMat*
   */
  TriangMat* triang(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on SymMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SymMat*
   */
  SymMat* sym(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on BandedMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a BandedMat*
   */
  BandedMat* banded(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on SparseMat matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a SparseMat*
   */
  SparseMat* sparse(unsigned int = 0, unsigned int = 0)const;

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

  /** set an identity matrix
   */
  void eye();

  /** resize the matrix with nbrow rows and nbcol columns, lower and upper are useful only for SparseMat.The existing elements of the Block matrix are preseved when specified.
   */
  void resize(unsigned int, unsigned int, unsigned int lower = 0, unsigned int upper = 0, bool = true);

  /** compute the infinite norm of the Block matrix
   *  \return a double
   */
  double normInf() const;

  /** display data on standard output
   */
  void display() const;

  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double& operator()(unsigned int, unsigned int);

  /** get or set the element matrix[i,j]
   *  \param an unsigned int i
   *  \param an unsigned int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double operator()(unsigned int, unsigned int) const;

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

  /** transpose in place: x->trans() is x = transpose of x.
   */
  void trans();

  /** transpose a matrix: x->trans(m) is x = transpose of m.
   *  \param a SiconosMatrix: the matrix to be transposed.
   */
  void trans(const SiconosMatrix&);

  /** get the vector tabRow
   *  \return a vector of int
   */
  inline Index getTabRow() const
  {
    return *_tabRow;
  };

  /** get the vector tabCol
   *  \return a vector of int
   */
  inline Index getTabCol() const
  {
    return *_tabCol;
  };

  /** get the vector tabRow
   *  \return a pointer to vector of int
   */
  inline const SP::Index tabRow() const
  {
    return _tabRow;
  };

  /** get the vector tabCol
   *  \return a pointer to vector of int
   */
  inline const SP::Index tabCol() const
  {
    return _tabCol;
  };

  /** get block at position row-col
   *  \param unsigned int row
   *  \param unsigned int col
   */
  inline SP::SiconosMatrix block(unsigned int row = 0, unsigned int col = 0)
  {
    return (*_mat)(row, col);
  };

  /** get block at position row-col
   *  \param unsigned int row
   *  \param unsigned int col
   */
  inline SPC::SiconosMatrix block(unsigned int row = 0, unsigned int col = 0) const
  {
    return cpp11ns::shared_ptr<SiconosMatrix>((*_mat)(row, col));
  };

  /** get row index of current matrix and save it unsigned into vOut
   *  \param unsigned int: index of required line
   *  \param ref to SiconosVector: in-out parameter
   */
  void  getRow(unsigned int, SiconosVector&) const;

  /** set line row of the current matrix with vector v
   *  \param an unsigned int and a SiconosVector
   */
  void  setRow(unsigned int, const SiconosVector&);

  /** get column index of current matrix and save it into vOut
   *  \param unsigned int: index of required column
   *  \param ref to SiconosVector: in-out parameter
   */
  void  getCol(unsigned int, SiconosVector&) const;

  /** set column col of the current matrix with vector v
   *  \param an unsigned int and a SiconosVector
   */
  void  setCol(unsigned int, const SiconosVector&);

  /** add a part of the input matrix (starting from (i,j) pos) to the current matrix
   *  \param an unsigned int i (in-out)
   *  \param an unsigned int j (in-out)
   *  \param a SiconosMatrix (in-out)
   */
  void addSimple(unsigned int&, unsigned int&, const SiconosMatrix&);

  /** subtract a part of the input matrix (starting from (i,j) pos) to the current matrix
   *  \param an unsigned int i (in-out)
   *  \param an unsigned int j (in-out)
   *  \param a SiconosMatrix (in-out)
   */
  void subSimple(unsigned int&, unsigned int&, const SiconosMatrix&);

  /** assignment
   *  \param SiconosMatrix : the matrix to be copied
   */
  BlockMatrix& operator = (const SiconosMatrix&);

  /** assignment
   *  \param BlockMatrix : the matrix to be copied
   */
  BlockMatrix& operator = (const BlockMatrix&);

  /** assignment
   *  \param  DenseMat: the matrix to be copied
   */
  BlockMatrix& operator = (const DenseMat&);

  /**operator +=
   *  \param SiconosMatrix : a matrix to add
   */
  BlockMatrix& operator +=(const SiconosMatrix&);

  /**operator -=
   *  \param SiconosMatrix : a matrix to subtract
   */
  BlockMatrix& operator -=(const SiconosMatrix&);

  /** computes an LU factorization of a general M-by-N matrix using partial pivoting with row interchanges.
   *  The result is returned in this (InPlace). Based on Blas dgetrf function.
   */
  void PLUFactorizationInPlace();

  /**  compute inverse of this thanks to LU factorization with Partial pivoting. This method inverts U and then computes inv(A) by solving the system
   *  inv(A)*L = inv(U) for inv(A). The result is returned in this (InPlace). Based on Blas dgetri function.
   */
  void PLUInverseInPlace();

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

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();
};

//DEFINE_SPTR(BlockMatrix)

#endif
