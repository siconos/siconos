/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

/*! \file BlockMatrix.h

*/

#ifndef __BlockMatrix__
#define __BlockMatrix__

#include "SiconosMatrix.h"

/** STDMAP
 * if STDMAP then we use an iterator of std::map, else we use an iterator1 and an iterator2 of map_std of Boost
 */
#define STDMAP 0

/** Object to handle block-matrices (ie lists of SimpleMatrix*)
 *
 * \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date (Creation) 21/07/2006
 *
 *  A block matrix is a map of pointers to SimpleMatrix, with two indices as a key.
 *  This indices give row and column position of the block in the global matrix.
 *
 * If block 1 is n1xm1, block2 n2xm2, block3 n3xm3 ..., then:\n
 *  tabRow = [ n1 n1+n2 n1+n2+n3 ...] \n
 *  tabCol = [ m1 m1+m2 m1+m2+m3 ...] \n
 *
 */
class BlockMatrix : public SiconosMatrix
{
private:

  /** a map of the blocks that compose the matrix
   */
  BlocksMat map;

  /** a list of bool, to check inside-class allocation for SimpleMatrix blocks.
   */
  std::deque<bool> isBlockAllocatedIn;

  /** list of blocks dimension - tabRow[i] = tabRow[i-1] + ni, ni being the number of rows of block i.
   */
  Index tabRow;

  /** list of blocks dimension - tabCol[i] = tabCol[i-1] + ni, ni being the number of columns of block i.
   */
  Index tabCol;

  /** build tabRow and tabCol
  *  \param unsigned int i: number of block-rows in this
  *  \param unsigned int i: number of block-columns in this
  */
  void makeTab(unsigned int, unsigned int);

  /** change tabRow and tabCol, according to a new input matrix of dim i x j
  *  \param unsigned int i: number of rows in the new matrix
  *  \param unsigned int i: number of columns in the new matrix
  */
  void addInTab(unsigned int, unsigned int);

  /** default constructor
  */
  BlockMatrix();

public:

  /** copy constructor
  *  \param SiconosMatrix
  */
  BlockMatrix(const SiconosMatrix&);

  /** copy constructor
  *  \param BlockMatrix
  */
  BlockMatrix(const BlockMatrix&);

  /** constructor with a map
  *  \param BlockMatrix
  */
  BlockMatrix(BlocksMat&);

  /** constructor with a list of pointer to SiconosMatrix (!links with pointer, no copy!)
  *  \param vector<SiconosMatrix*>
  *  \param unsigned int: number of blocks in a row
  *  \param unsigned int: number of col in a row
  */
  BlockMatrix(const std::vector<SiconosMatrix* >&, unsigned int, unsigned int);

  /** contructor with a list of 4 pointer to SiconosMatrix (!links with pointer, no copy!)
  *  \param SiconosMatrix * m1, block (0,0)
  *  \param SiconosMatrix * m2, block (0,1)
  *  \param SiconosMatrix * m3, block (1,0)
  *  \param SiconosMatrix * m4, block (1,1)
  */
  BlockMatrix(SiconosMatrix*, SiconosMatrix*, SiconosMatrix*, SiconosMatrix*);

  /** destructor
  */
  ~BlockMatrix(void);

  /** Computes dim according to the matrix type.
   */
  void computeDim();

  /** get the number of blocks in a row
   *  \return an unsigned int
   */
  inline const unsigned int numberOfBlockInARow(void) const
  {
    return tabRow.size();
  };

  /** get the number of blocks in a column
  *  \return an unsigned int
  */
  inline const unsigned int numberOfBlockInACol(void) const
  {
    return tabCol.size();
  };

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

  /** resize the matrix with nbrow rows and nbcol columns, lower and upper are useful only for SparseMat.The existing elements of the Block matrix are preseved when specified.
  *  \exception SiconosMatrixException
  */
  void resize(unsigned int, unsigned int, unsigned int lower = 0, unsigned int upper = 0, bool = true);

  /** determines if the matrix is square
  *  \return true if the matrix is square
  */
  inline bool isSquare() const
  {
    return (size(0) == size(1));
  };

  /** determines if the matrix has been inversed
  *  \return true if the matrix is inversed
  */
  inline bool isInversed() const
  {
    return false;
  };

  /** determines if the matrix has been factorized
  *  \return true if the matrix is factorized
  */
  inline bool isFactorized() const
  {
    return false;
  };

  /** return the adress of the array of double values of the matrix
  *  \param: row position for the required block ->useless for SimpleMatrix
  *  \param: col position for the required block ->useless for SimpleMatrix
  *  \return double* : the pointer on the double array
  */
  double* getArray(unsigned int = 0, unsigned int = 0) const;

  /** compute the infinite norm of the Block matrix
  *  \return a double
  */
  const double normInf(void)const;

  /** transpose in place: x->trans() is x = transpose of x.
   */
  void trans();

  /** transpose a matrix: x->trans(m) is x = transpose of m.
   *  \param a SiconosMatrix: the matrix to be transposed.
   */
  void trans(const SiconosMatrix&);

  /** sets all the values of the matrix to 0.0
  */
  void zero(void);

  /** set an identity matrix
  */
  void eye(void);

  /** display data on standard output
  */
  void display(void)const;

  // ******************************* GETTERS/SETTERS **************************

  /** get the vector tabRow
  *  \return a vector of int
  */
  inline Index getTabRow() const
  {
    return tabRow;
  };

  /** get the vector tabCol
  *  \return a vector of int
  */
  inline Index getTabCol() const
  {
    return tabCol;
  };

  /** get block(row,col)
  *  \param 2 unsigned int for indexes and a SiconosMatrix (in-out paramater)
  */
  void  getBlock(unsigned int, unsigned int, SiconosMatrix&)const;

  /** get block at position row-col
  *  \param unsigned int row
  *  \param unsigned int col
  */
  SiconosMatrix* getBlockPtr(unsigned int = 0, unsigned int = 0);

  /** get std::deque of bool
  *   \useless for SimpleMatrix
  *   \return a std::deque<bool>
  */
  const std::deque<bool> getBlockAllocated(void)const;

  /** get the attribute num of current matrix, useless for Block Matrix
  * \return an unsigned int.
  */
  unsigned int getNum()const;

  /** get the number of block (i=0, row, i=1 col)
  *  \param unsigned int(i=0, row, i=1 col)
  *  \return an unsigned int
  */
  unsigned int getNumberOfBlocks(unsigned int) const;

  /** get row index of current matrix and save it unsigned into vOut
  *  \param unsigned int: index of required line
  *  \param ref to SimpleVector: in-out parameter
  */
  void  getRow(unsigned int, SimpleVector&) const;

  /** set line row of the current matrix with vector v
  *  \param an unsigned int and a SimpleVector
  */
  void  setRow(unsigned int, const SimpleVector&);

  /** get column index of current matrix and save it into vOut
  *  \param unsigned int: index of required column
  *  \param ref to SimpleVector: in-out parameter
  */
  void  getCol(unsigned int, SimpleVector&) const;

  /** set column col of the current matrix with vector v
  *  \param an unsigned int and a SimpleVector
  */
  void  setCol(unsigned int, const SimpleVector&);

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
  DenseMat* getDensePtr(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on TriangMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a TriangMat*
  */
  TriangMat* getTriangPtr(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on SymMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a SymMat*
  */
  SymMat* getSymPtr(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on BandedMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a BandedMat*
  */
  BandedMat* getBandedPtr(unsigned int = 0, unsigned int = 0)const;

  /** get a pointer on SparseMat matrix
  *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
  *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
  *  \return a SparseMat*
  */
  SparseMat* getSparsePtr(unsigned int = 0, unsigned int = 0)const;

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

  /** get the objects that holds all the blocks.
  *  \return a BlocksMat
  */
  const BlocksMat getAllBlocks(void)const;

  //******************* MATRICES HANDLING AND OPERATORS ***********************


  /** copy the the matrix M in the matrix "mat" at the position (x,y) (! block coordinates)
  *  \param SiconosMatrix& : the matrix to be copied in the current matrix
  *  \param unsigned int : the block-line position to start the copy of the blockmatrix
  *  \param unsigned int : the block-column position to start the copy of the blockmatrix
  */
  void matrixCopy(const SiconosMatrix&, unsigned int, unsigned int);

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
  double operator()(unsigned int, unsigned int)const;

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

  /**operator =
   *  \param SiconosMatrix : the matrix to be copied
   */
  BlockMatrix& operator = (const SiconosMatrix&);

  /**operator =
   *  \param BlockMatrix : the matrix to be copied
   */
  //  BlockMatrix& operator = (const BlockMatrix&);

  /**operator +=
   *  \param SiconosMatrix : a matrix to add
   */
  BlockMatrix& operator +=(const SiconosMatrix&);

  /**operator -=
   *  \param SiconosMatrix : a matrix to subtract
   */
  BlockMatrix& operator -=(const SiconosMatrix&);

  /**operator /=
   *  \param double, a scalar
   */
  BlockMatrix& operator /=(double);

  /**operator /=
   *  \param int, a scalar
   */
  BlockMatrix& operator /=(int);

  /**operator *=
   *  \param double, a scalar
   */
  BlockMatrix& operator *=(double);

  /**operator *=
   *  \param int, a scalar
   */
  BlockMatrix& operator *=(int);

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
};

#endif
