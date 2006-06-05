/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
/** \class BlockMatrix
 *   \brief list of SimpleMatrix * to work with "block-matrices"
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.2.0.
 *   \date (Creation) 05/19/2004
 *
 *  A block matrix is a map of pointers to SimpleMatrix, with a vector of two indices as a key.
 *  This vector gives row and column position of the block in the global matrix.
 *
 * If block 1 is n1xm1, block2 n2xm2, block3 n3xm3 ..., then:
 *  tabRow = [ n1 n1+n2 n1+n2+n3 ...]
 *  tabCol = [ m1 m1+m2 m1+m2+m3 ...]
 *
 */

#ifndef __BlockMatrix__
#define __BlockMatrix__

#include "SimpleMatrix.h"
#include<vector>
#include<map>
#include<deque>

class SimpleVector;
class SimpleMatrix;

typedef std::map< std::vector<unsigned int>, SiconosMatrix*> blocksMap;

class BlockMatrix: public SiconosMatrix
{
private:

  /** \var matrixOfBlocks
   * \brief a map of the blocks that compose the matrix
   */
  blocksMap matrixOfBlocks;

  /** \var isBlockAllocatedIn
   * \brief a list of bool, to check inside-class allocation for SimpleMatrix blocks.
   */
  std::deque<bool> isBlockAllocatedIn;

  /** \var tabRow
   * \brief list of blocks dimension (number of rows)
   */
  std::vector<unsigned int> tabRow;

  /** \var tabRow
   * \brief list of blocks dimension (number of columns)
   */
  std::vector<unsigned int> tabCol;

  /** \fn BlockMatrix ()
   *  \brief contructor
   */
  BlockMatrix();

  /** \fn addInTab(const unsigned int& i, const unsigned int& j)
   *  \brief change tabRow and tabCol, according to a new input matrix of dim ixj
   *  \param unsigned int i: number of rows in the new matrix
   *  \param unsigned int i: number of columns in the new matrix */
  void addInTab(const unsigned int&, const unsigned int&);

  /** \fn void BlockMatrix::makeTab(const unsigned int& row, const unsigned int& col)
   *  \brief build tabRow and tabCol
   *  \param unsigned int i: number of block-rows in this
   *  \param unsigned int i: number of block-columns in this */
  void makeTab(const unsigned int& row, const unsigned int& col);

public:

  /** \fn BlockMatrix (const SiconosMatrix&  m)
   *  \brief copy contructor
   *  \param SiconosMatrix
   */
  BlockMatrix(const SiconosMatrix&);

  /** \fn BlockMatrix (const BlockMatrix&  m)
   *  \brief copy contructor
   *  \param BlockMatrix
   */
  BlockMatrix(const BlockMatrix&);

  /** \fn BlockMatrix (const vector<LaGenMatDouble>& v, const unsigned int& row, const unsigned int& col)
   *  \brief contructor with a list of LaGenMatDouble
   *  \param vector<LaGenMatDouble>
   *  \param unsigned int: number of blocks in a row
   *  \param unsigned int: number of col in a row
   */
  BlockMatrix(const std::vector<LaGenMatDouble>&, const unsigned int&, const unsigned int&);

  /** \fn BlockMatrix (vector<SimpleMatrix*>v, const unsigned int& row, const unsigned int& col)
   *  \brief contructor with a list of pointer to SimpleMatrix (!links with pointer, no copy!)
   *  \param vector<SimpleMatrix*>
   *  \param unsigned int: number of blocks in a row
   *  \param unsigned int: number of col in a row
   */
  BlockMatrix(std::vector<SiconosMatrix*>, const unsigned int&, const unsigned int&);

  /** \fn BlockMatrix(SiconosMatrix* m1, SiconosMatrix* m2, SiconosMatrix* m3, SiconosMatrix* m4);
   *  \brief contructor with a list of 4 pointer to SiconosMatrix (!links with pointer, no copy!)
   *  \param SiconosMatrix * m1, block (0,0)
   *  \param SiconosMatrix * m2, block (0,1)
   *  \param SiconosMatrix * m3, block (1,0)
   *  \param SiconosMatrix * m4, block (1,1)
   */
  BlockMatrix(SiconosMatrix*, SiconosMatrix*, SiconosMatrix*, SiconosMatrix*);

  /** \fn ~BlockMatrix ()
   *  \brief destructor
   */
  ~BlockMatrix();

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
  bool isInversed() const;

  /** \fn bool isfactorized()
   *  \brief determines if the matrix has been factorized
   *  \return true if the matrix is factorized
   */
  bool isFactorized() const;

  // --- GETTERS/SETTERS ---

  /** \fn LaGenMatDouble getLaGenMatDouble(const int& row = 0, const int& col = 0)
   *  \brief get LaGenMatDouble matrix
   *  \param an unsigned int, position of the block (row) - Useless for BlockMatrix
   *  \param an unsigned int, position of the block (column) - Useless for BlockMatrix
   *  \return a LaGenMatDouble
   */
  const LaGenMatDouble getLaGenMatDouble(const unsigned int& = 0, const unsigned int& = 0) const;

  /** \fn LaGenMatDouble* getLaGenMatDoubleRef(const int& row = 0, const int& col = 0)
   *  \brief get LaGenMatDouble matrix
   *  \param an unsigned int, position of the block (row) - Useless for SimpleMatrix
   *  \param an unsigned int, position of the block (column) - Useless for SimpleMatrix
   *  \return a LaGenMatDouble
   */
  const LaGenMatDouble* getLaGenMatDoubleRef(const unsigned int& = 0, const unsigned int& = 0) const;

  /** \fn inline blocksMap getListOfBlocks() const {return matrixOfBlocks;};
   *  \brief get the map of blocks
   *  \return a blocksMap
   */
  inline blocksMap getListOfBlocks() const
  {
    return matrixOfBlocks;
  };

  /** \fn inline std::vector<unsigned int> getTabRow() const
   *  \brief get the vector tabRow
   *  \return a vector of int
   */
  inline std::vector<unsigned int> getTabRow() const
  {
    return tabRow;
  };

  /** \fn inline std::vector<unsigned int> getTabCol() const
   *  \brief get the vector tabCol
   *  \return a vector of int
   */
  inline std::vector<unsigned int> getTabCol() const
  {
    return tabCol;
  };

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
   * \brief set mat of a BlockMatrix, of mat-block(row,col) of a BlockMatrix
   * \param an unsigned int, position of the block (row) - Useless for BlockMatrix
   * \param an unsigned int, position of the block (column) - Useless for BlockMatrix
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

  /** \fn SiconosMatrix* getBlockPtr(const unsigned int& row, const unsigned int& col)
   *  \brief get block at position row-col
   *  \param unsigned int row
   *  \param unsigned int col
   */
  SiconosMatrix* getBlockPtr(const unsigned int& = 0, const unsigned int& = 0);

  /** \fn double* getArray(const unsigned int& i,const unsigned int&  j)
   *  \brief return the adress of the array of double values of the matrix for block(i,j)
   *  \param: row position for the required block
   *  \param: col position for the required block
   *  \return double* : the pointer on the double array
   */
  double* getArray(const unsigned int&, const unsigned int&);

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
   *  \return BlockMatrix : the result of the multiplication
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
  BlockMatrix& operator = (const SiconosMatrix&);

  /** \fn assignment operator
   *  \param BlockMatrix : the matrix to be copied
   */
  BlockMatrix& operator = (const BlockMatrix&);

  /** \fn operator +=
  *  \param  BlockMatrix: a matrix to add
  */
  BlockMatrix& operator+=(const SiconosMatrix  &) ;

  /** \fn operator -=
  *  \param BlockMatrix : a matrix to subtract
  */
  BlockMatrix& operator-=(const SiconosMatrix  &);

  /** \fn operator *=
  *  \param double, a scalar
  */
  BlockMatrix& operator*=(const double  &);

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the matrix
   *  \return a double
   */
  const double normInf() const;

  // --- COMPUTING WITH MATRICES  ---

  /** \fn void linearSolve(const SiconosMatrix & B, BlockMatrix & X)
   *  \brief compute the vector x in  Ax=B
   *  \param BlockMatrix & B
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

};

#endif // __BlockMatrix__
