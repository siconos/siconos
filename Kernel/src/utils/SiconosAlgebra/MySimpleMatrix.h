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
/** \class MySimpleMatrix
 *   \brief This class is an encapsulation of the Boost class managing matrices of double.
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.3.0.
 *   \date (Creation) 07/21/2006
 *
 *
 *
 *
 * MySimpleMatrix is used in the platform to store matrices (mathematical object).
 *
 */

#ifndef __MySimpleMatrix__
#define __MySimpleMatrix__

#include "MySiconosMatrix.h"

using namespace boost::numeric::ublas;
class MySimpleMatrix: public MySiconosMatrix
{
private:

  /** \var int num
   * \brief an int which make a correspondance with Boost Matrices: 1 -> DenseMat, 2 -> TriangMat, 3 -> SymMat, 4->SparseMat, 5->BandedMat
   */
  int num;

  /** \var MyMat mat (See MySiconosMatrix.h for more details on MyMat type);
   * \brief union of The Boost Matrices : DenseMat, TriangMat, SymMat are encapsulated.
   */
  MyMat mat;

  /** \fn MySimpleMatrix ()
   *  \brief default constructor
   *  \return MySimpleMatrix
   */
  MySimpleMatrix(void);

public:
  /***************************** CONSTRUCTORS ****************************/

  /** \fn MySimpleMatrix (const MySimpleMatrix&)
   *  \brief copy constructor
   *  \param MySimpleMatrix
   */
  MySimpleMatrix(const MySimpleMatrix&);

  /** \fn MySimpleMatrix (const MySiconosMatrix&)
   *  \brief copy constructor
   *  \param MySiconosMatrix
   */
  MySimpleMatrix(const MySiconosMatrix&);

  /** \fn MySimpleMatrix (TYP)
   *  \brief constructor with the type of the Boost matrix
   *  \param TYP
   */
  MySimpleMatrix(TYP);

  /** \fn MySimpleMatrix (TYP, int, int)
   *  \brief constructor with the type and the dimension of the Boost matrix
   *  \param TYP
   *  \param 2 int
   */
  MySimpleMatrix(TYP, int, int);

  /** \fn MySimpleMatrix (const DenseMat&)
   *  \brief constructor with a DenseMat matrix (see MySiconosMatrix.h for details)
   *  \param a DenseMat
   */
  MySimpleMatrix(const DenseMat&);

  /** \fn MySimpleMatrix (const TriangMat&)
   *  \brief constructor with a TriangMat matrix (see MySiconosMatrix.h for details)
   *  \param a TriangMat
   */
  MySimpleMatrix(const TriangMat&);

  /** \fn MySimpleMatrix (const SymMat&)
   *  \brief constructor with a SymMat matrix (see MySiconosMatrix.h for details)
   *  \param a SymMat
   */
  MySimpleMatrix(const SymMat&);

  /** \fn MySimpleMatrix (const BandedMat&)
   *  \brief constructor with a BandedMat matrix (see MySiconosMatrix.h for details)
   *  \param a BandedMat
   */
  MySimpleMatrix(const BandedMat&);

  /** \fn MySimpleMatrix (const SparseMat&)
   *  \brief constructor with a SparseMat matrix (see MySiconosMatrix.h for details)
   *  \param a SparseMat
   */
  MySimpleMatrix(const SparseMat&);


  /** \fn MySimpleMatrix (TYP, const std::vector<double>&, int, int=0)
   *  \brief constructor with the type of the boost matrix, a vector of the values and the dimensions
   *  of the matrix, the integers upper and lower are useful only for BandedMat
   *  \param a TYP
   *  \param a std::vector<double>
   *  \param 4 int
   */
  MySimpleMatrix(TYP , const std::vector<double>& , int row, int col = 0, int lower = 0, int upper = 0);

  /** \fn MySimpleMatrix (std::string file, bool ascii)
   *  \brief constructor with an input file
   *  \param a std::string which contain the file path
   *  \param a boolean to indicate if the file is in ascii
   */
  MySimpleMatrix(const std::string&, bool = true);

  /** \fn ~MySimpleMatrix ()
   *  \brief destructor
   */
  ~MySimpleMatrix(void);
  //************************** GETTERS/SETTERS  **************************

  /** \fn int getNum() const
   *  \brief get the attribute num of current matrix, useless for Block Matrix
   * \return an int.
   */
  int getNum(void) const;

  /** \fn DenseMat* getDense(int row = 0, int col = 0)
   *  \brief get DenseMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a DenseMat
   */
  const DenseMat getDense(int = 0, int = 0) const;

  /** \fn TriangMat getTriang(int row = 0, int col = 0)
   *  \brief get TriangMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a TriangMat
   */
  const TriangMat getTriang(int = 0, int = 0) const;

  /** \fn SymMat getSym(int row = 0, int col = 0)
   *  \brief get SymMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SymMat
   */
  const SymMat getSym(int = 0, int = 0) const;

  /** \fn BandedMat getBanded(int row = 0, int col = 0)
   *  \brief get BandedMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a BandedMat
   */
  const BandedMat getBanded(int = 0, int = 0) const;


  /** \fn SparseMat getSparse(int row = 0, int col = 0)
   *  \brief get SparseMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SparseMat
   */
  const SparseMat getSparse(int = 0, int = 0) const;

  /** \fn DenseMat* getDensePtr(int row = 0, int col = 0)
   *  \brief get a pointer on DenseMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a DenseMat*
   */
  const DenseMat* getDensePtr(int = 0, int = 0) const;

  /** \fn TriangMat* getTriangPtr(int row = 0, int col = 0)
   *  \brief get a pointer on TriangMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a TriangMat*
   */
  const TriangMat* getTriangPtr(int = 0, int = 0) const;

  /** \fn SymMat* getSymPtr(int row = 0, int col = 0)
   *  \brief get a pointer on SymMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SymMat*
   */
  const SymMat* getSymPtr(int = 0, int = 0) const;

  /** \fn BandedMat* getBandedPtr(int row = 0, int col = 0)
   *  \brief get a pointer on BandedMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a BandedMat*
   */
  const BandedMat* getBandedPtr(int = 0, int = 0) const;

  /** \fn SparseMat* getSparsePtr(int row = 0, int col = 0)
   *  \brief get a pointer on SparseMat matrix
   *  \param an int, position of the block (row) - Useless for MySimpleMatrix
   *  \param an int, position of the block (column) - Useless for MySimpleMatrix
   *  \return a SparseMat*
   */
  const SparseMat* getSparsePtr(int = 0, int = 0) const;

  /** \fn mapped getMap()
   *  \brief get mapped matrix
   *  \useless for MySimpleMatrix
   *  \return a mapped
   */
  const mapped getMap(void) const;

  /** \fn void getBlock(int numCow, int numCol, MySiconosMatrix& block)
   *  \brief get block corresponding to lines given in numRow and columns in numCol
   *  \param 2 int for indexes and a MySiconosMatrix (in-out paramater)
   */
  void getBlock(int, int, MySiconosMatrix&) const;

  /** \fn const std::deque<bool> getBlockAllocated()const
   *   \brief get std::deque of bool
   *   \useless for MySimpleMatrix
   *   \return a std::deque<bool>
   */
  const std::deque<bool> getBlockAllocated(void) const;

  /** \fn void getRow(int index, MySimpleVector& vOut) const
   *  \brief get row index of current matrix and save it into vOut
   *  \param int: index of required line
   *  \param ref to MySimpleVector: in-out parameter
   */
  void getRow(int, MySimpleVector&) const;

  /** \fn void getCol(int index, MySimpleVector& vOut) const
   *  \brief get column index of current matrix and save it into vOut
   *  \param int: index of required column
   *  \param ref to MySimpleVector: in-out parameter
   */
  void getCol(int, MySimpleVector&) const;

  /** \fn void setNum(int n)
   *  \brief set the attribute num of current matrix with n
   */
  void setNum(int);

  /** \fn void setRow(int row, const MySimpleVector &v)
   *  \brief set line row of the current matrix with vector v
   *  \param an int and a MySimpleVector
   */
  void setRow(int, const MySimpleVector&);

  /** \fn void setCol(int col, const MySimpleVector &v)
   *  \brief set column col of the current matrix with vector v
   *  \param an int and a MySimpleVector
   */
  void setCol(int, const MySimpleVector&);

  /** \fn void zero();
   *  \brief sets all the values of the matrix to 0.0
   */
  void zero(void);

  /** \fn void eye();
   *  \brief set an identity matrix
   */
  void eye(void);

  /** \fn  int size1 (void)const
   *  \brief get the number of rows of the matrix
   *  \exception SiconosMatrixException
   *  \return the number of rows of the matrix
   */
  int size1(void) const;

  /** \fn  int size2 (void)const
   *  \brief get the number of columns of the matrix
   *  \exception SiconosMatrixException
   *  \return the number of columns of the matrix
   */
  int size2(void) const;

  /** \fn  void resize (int nbrow, int nbcol, bool val = true)const
   *  \brief resize the matrix with nbrow rows and nbcol columns The existing elements of the matrix are preseved when specified.
   *  \exception SiconosMatrixException
   */
  void resize(int, int, int = 0, int = 0, bool = true);

  /** \fn const double normInf() const;
   *  \brief compute the infinite norm of the matrix
   *  \return a double
   */
  const double normInf(void) const;

  /** \fn void display();
   *  \brief display data on standard output
   */
  void display(void) const;

  // --- MATRICES HANDLING AND OPERATORS ---

  /** \fn void blockMatrixCopy( MySiconosMatrix &blockMat, int, int)
   *  \brief copy the blockmatrix "blockMat" in the matrix "mat" at the position (xPos, yPos)
   *                    [0, 0, 0, 0].blockMatrixCopy([1], 0, 2) => mat = [0 0 1 0]
   *  \param MySiconosMatrix& : the block matrix to copy in the current matrix
   *  \param int : the line position to start the copy of the blockmatrix
   *  \param int : the column position to start the copy of the blockmatrix
   */
  void blockMatrixCopy(const MySiconosMatrix&, int, int);

  /** \fn double& operator() (int row, int col)
   *  \brief get or set the element matrix[i,j]
   *  \param an int i
   *  \param an int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double& operator()(int , int);

  /** \fn double operator() (int row, int col)const
   *  \brief get or set the element matrix[i,j]
   *  \param an int i
   *  \param an int j
   *  \exception SiconosMatrixException
   *  \return the element matrix[i,j]
   */
  double operator()(int , int) const;

  /** \fn assignment operator
   *  \param MySiconosMatrix : the matrix to be copied
   */
  const MySimpleMatrix& operator = (const MySiconosMatrix&);

  /** \fn assignment operator
   *  \param MySimpleMatrix : the matrix to be copied
   */
  const MySimpleMatrix& operator = (const MySimpleMatrix&);

  /** \fn operator /=
   *  \param double, a scalar
   */
  const MySimpleMatrix& operator /= (double);

  /** \fn operator /=
   *  \param int, a scalar
   */
  const MySimpleMatrix& operator /= (int);

  /** \fn operator +=
   *  \param MySiconosMatrix : a matrix to add
   */
  const MySimpleMatrix& operator +=(const MySiconosMatrix&);

  /** \fn operator -=
   *  \param MySiconosMatrix : a matrix to subtract
   */
  const MySimpleMatrix& operator -=(const MySiconosMatrix&);

  /** \fn operator *=
   *  \param double, a scalar
   */
  const MySimpleMatrix& operator *= (double);

  /** \fn operator *=
   *  \param int, a scalar
   */
  const MySimpleMatrix& operator *= (int);

  /** \fn operator ==
   * \brief: A==B when (A-B).normInf()<tolerance
   * \param 2 MySiconosMatrix
   * \return a boolean
   */
  friend bool operator == (const MySiconosMatrix&, const MySiconosMatrix&);

  /** \fn operator + (const MySiconosMatrix& m1, const MySiconosMatrix& m2);
   *  \brief Addition of two matrices
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   *  \exception SiconosMatrixException, if the two matrices have different types, in this case use function add
   */
  friend MySimpleMatrix operator +(const MySiconosMatrix&, const MySiconosMatrix&);

  /** \fn operator - (const MySiconosMatrix& m1, const MySiconosMatrix& m2);
   *  \brief subtraction of two matrices
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   *  \exception SiconosMatrixException, if the two matrices have different types, in this case use function sub
   */
  friend MySimpleMatrix operator - (const MySiconosMatrix&, const MySiconosMatrix&);

  /** \fn operator * (const MySiconosMatrix& m1, const MySiconosMatrix& m2);
   *  \brief multiplication of two matrices
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleMatrix
   *  \exception SiconosMatrixException, if the two matrices have different types, in this case use function prod
   */
  friend MySimpleMatrix operator *(const MySiconosMatrix&, const MySiconosMatrix&);

  /** \fn operator * (const MySiconosMatrix& m1, double d);
   *  \brief multiplication of a matrix by a double
   *  \param a MySiconosMatrix
   *  \param a double
   *  \return a MySimpleMatrix
   */
  friend MySimpleMatrix operator * (const MySiconosMatrix&, double);


  /** \fn operator * (const MySiconosMatrix& m1, int d);
   *  \brief multiplication of a matrix by an int
   *  \param a MySiconosMatrix
   *  \param an int
   *  \return a MySimpleMatrix
   */
  friend MySimpleMatrix operator *(const MySiconosMatrix&, int);

  /** \fn operator * (double d, const MySiconosMatrix& m1);
   *  \brief multiplication of a matrix by a double
   *  \param a double
   *  \param a MySiconosMatrix
   *  \return a MySimpleMatrix
   */
  friend MySimpleMatrix operator * (double , const MySiconosMatrix&);


  /** \fn operator * (int d, const MySiconosMatrix& m1);
   *  \brief multiplication of a matrix by an int
   *  \param an int
   *  \param a MySiconosMatrix
   *  \return a MySimpleMatrix
   */
  friend MySimpleMatrix operator *(int, const MySiconosMatrix&);

  /** \fn operator / (const MySiconosMatrix& m1, double d);
   *  \brief division of the matrix by a double
   *  \param a MySiconosMatrix
   *  \param a double
   *  \return a MySimpleMatrix
   *  \exception SiconosMatrixException, if the double d = 0
   */
  friend MySimpleMatrix operator /(const MySiconosMatrix&, double);

  /** \fn operator / (const MySiconosMatrix& m1, int d);
   *  \brief division of the matrix by an int
   *  \param a MySiconosMatrix
   *  \param an int
   *  \return a MySimpleMatrix
   *  \exception SiconosMatrixException, if the int d = 0
   */
  friend MySimpleMatrix operator / (const MySiconosMatrix&, int);

  /** \fn trans (const MySiconosMatrix& m1);
   *  \brief transpose the matrix m1
   *  \param a MySiconosMatrix
   *  \return a MySimpleMatrix
   */
  friend MySimpleMatrix trans(const MySiconosMatrix&);

  /** \fn add (const MySiconosMatrix& m1, const MySiconosMatrix& m2);
   *  \brief Addition of two matrices
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   */
  friend MySimpleMatrix add(const MySiconosMatrix&, const MySiconosMatrix&);

  /** \fn sub (const MySiconosMatrix& m1, const MySiconosMatrix& m2);
   *  \brief subtraction of two matrices
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleMatrix
   *  \exception SiconosMatrixException, if the sizes are incompatible
   */
  friend MySimpleMatrix sub(const MySiconosMatrix&, const MySiconosMatrix&);

  /** \fn prod (const MySiconosMatrix& m1, const MySiconosMatrix& m2);
   *  \brief multiplication of two matrices
   *  \param 2 MySiconosMatrix
   *  \return a MySimpleMatrix
   */
  friend MySimpleMatrix prod(const MySiconosMatrix&, const MySiconosMatrix&);

  /** \fn MySimpleMatrix multTranspose(const MySiconosMatrix &A, const MySiconosMatrix &B)
   *  \brief compute A*Bt
   *  \param 2 MySiconosMatrix
   *  \return MySimpleMatrix : the result of the multiplication
   */
  friend MySimpleMatrix multranspose(const MySiconosMatrix&, const MySiconosMatrix&);
};
#endif
