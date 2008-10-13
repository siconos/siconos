/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>

#include "KernelConfig.h"
#ifdef HAVE_ATLAS
#include <boost/numeric/bindings/atlas/clapack.hpp>
namespace lapack = boost::numeric::bindings::atlas;
#else
#include <boost/numeric/bindings/lapack/lapack.hpp>
namespace lapack = boost::numeric::bindings::lapack;
#endif

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#include "SimpleMatrix.h"
#include "SimpleVector.h"
#include "ioMatrix.h"






// =================================================
//                CONSTRUCTORS
// =================================================

using std::cout;
using std::endl;

// parameters: dimensions and type.
SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
  SiconosMatrix(1, row, col), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  if (typ == DENSE)
  {
    mat.Dense = new DenseMat(ublas::zero_matrix<double>(row, col));
    // num = 1; default value
  }
  else if (typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat(ublas::zero_matrix<double>(row, col));
    num = 2;
  }
  else if (typ == SYMMETRIC)
  {
    mat.Sym = new SymMat(ublas::zero_matrix<double>(row, col));
    num = 3;
  }
  else if (typ == SPARSE)
  {
    mat.Sparse = new SparseMat(row, col, upper);
    num = 4;
    zero();
  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat(row, col, upper, lower);
    num = 5;
    zero();
  }
  else if (typ == ZERO)
  {
    mat.Zero = new ZeroMat(row, col);
    num = 6;
  }
  else if (typ == IDENTITY)
  {
    mat.Identity = new IdentityMat(row, col);
    num = 7;
  }
  else
    SiconosMatrixException::selfThrow("SiconosMatrix::constructor(UBLAS_TYPE type, unsigned int row, unsigned int col): invalid type.");
}

// parameters: dimensions, input value and type
SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, double inputValue, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
  SiconosMatrix(1, row, col), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  // This constructor has sense only for dense matrices ...
  if (typ == DENSE)
  {
    mat.Dense = new DenseMat(ublas::scalar_matrix<double>(row, col, inputValue));
    // num = 1; default value
  }
  else
    SiconosMatrixException::selfThrow("SiconosMatrix::constructor(UBLAS_TYPE type, unsigned int row, unsigned int col, double fillInValue): invalid type.");
}

// // parameters: a vector (stl) of double and the type.
// SimpleMatrix::SimpleMatrix(const std::vector<double>& v, unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int lower, unsigned int upper):
//   SiconosMatrix(1, row, col), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
// {
//   if( (  (v.size() != row*col) && (typ != SYMMETRIC && typ != BANDED) )
//       || (v.size() != row*row && typ == SYMMETRIC)
//       || (typ == BANDED && ( (v.size()) != (unsigned int)(std::max)(row, col)*(lower+1+upper) ) ))
//     SiconosMatrixException::selfThrow("constructor(UBLAS_TYPE, const std::vector<double>, int, int) : invalid vector size");

//   if(typ == DENSE)
//     {
//       mat.Dense = new DenseMat(row,col);
//       // num = 1; default value
//     }
//   else if(typ == TRIANGULAR)
//     {
//       mat.Triang = new TriangMat(row,col);
//       num = 2;
//     }
//   else if(typ == SYMMETRIC)
//     {
//       mat.Sym = new SymMat(row);
//       num = 3;
//     }
//   else if(typ == SPARSE)
//     {
//       SiconosMatrixException::selfThrow("SimpleMatrix::constructor(UBLAS_TYPE, const std::vector<double>, int row, int col, int lower, int upper) : warning -- use constructor(const SparseMat &m) or constructor(UBLAS_TYPE, int row, int col) with UBLAS_TYPE = SPARSE");

//     }
//   else if(typ == BANDED)
//     {
//       mat.Banded = new BandedMat(row, col, lower, upper);
//       num = 5;
//     }
//   else
//     SiconosMatrixException::selfThrow("constructor(UBLAS_TYPE, const std::vector<double>, int, int) : invalid type of matrix given");

//   std::copy(v.begin(), v.end(), (vect.Dense)->begin());


// }

// Copy constructors
SimpleMatrix::SimpleMatrix(const SimpleMatrix &smat): SiconosMatrix(smat.getNum(), smat.size(0), smat.size(1)), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  if (num == 1)
  {
    mat.Dense = new DenseMat(smat.size(0), smat.size(1));
    noalias(*mat.Dense) = (*smat.getDensePtr());
  }
  //   mat.Dense = new DenseMat(*smat.getDensePtr());

  else if (num == 2)
    mat.Triang = new TriangMat(*smat.getTriangPtr());

  else if (num == 3)

    mat.Sym = new SymMat(*smat.getSymPtr());

  else if (num == 4)
    mat.Sparse = new SparseMat(*smat.getSparsePtr());

  else if (num == 5)
    mat.Banded = new BandedMat(*smat.getBandedPtr());

  else if (num == 6)
    mat.Zero = new ZeroMat(smat.size(0), smat.size(1));

  else// if(num == 7)
    mat.Identity = new IdentityMat(smat.size(0), smat.size(1));
}

SimpleMatrix::SimpleMatrix(const SiconosMatrix &m): SiconosMatrix(m.getNum(), m.size(0), m.size(1)), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  // num is set in SiconosMatrix constructor with m.getNum() ... must be changed if m is Block
  unsigned int numM = m.getNum();
  if (numM == 0) // ie if m is Block, this matrix is set to a dense.
  {
    num = 1;
    // get number of blocks in a row/col of m.
    mat.Dense = new DenseMat(dimRow, dimCol);
    ConstBlockIterator1 it;
    ConstBlockIterator2 it2;
    unsigned int posRow = 0;
    unsigned int posCol = 0;

    for (it = m.begin(); it != m.end(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        setBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if (num == 1)
  {
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    noalias(*mat.Dense) = (*m.getDensePtr());
  }

  else if (num == 2)
    mat.Triang = new TriangMat(*m.getTriangPtr());

  else if (num == 3)
    mat.Sym = new SymMat(*m.getSymPtr());

  else if (num == 4)
    mat.Sparse = new SparseMat(*m.getSparsePtr());

  else if (num == 5)
    mat.Banded = new BandedMat(*m.getBandedPtr());

  else if (num == 6)
    mat.Zero = new ZeroMat(m.size(0), m.size(1));

  else // if(num == 7)
    mat.Identity = new IdentityMat(m.size(0), m.size(1));
}

SimpleMatrix::SimpleMatrix(const DenseMat& m): SiconosMatrix(1, m.size1(), m.size2()), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Dense = new DenseMat(m);
}

SimpleMatrix::SimpleMatrix(const TriangMat& m): SiconosMatrix(2, m.size1(), m.size2()), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Triang = new TriangMat(m);
}

SimpleMatrix::SimpleMatrix(const SymMat& m): SiconosMatrix(3, m.size1(), m.size2()), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Sym = new SymMat(m);
}

SimpleMatrix::SimpleMatrix(const SparseMat& m): SiconosMatrix(4, m.size1(), m.size2()), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Sparse = new SparseMat(m);
}

SimpleMatrix::SimpleMatrix(const BandedMat& m): SiconosMatrix(5, m.size1(), m.size2()), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Banded = new BandedMat(m);
}

SimpleMatrix::SimpleMatrix(const ZeroMat& m): SiconosMatrix(6, m.size1(), m.size2()), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Zero = new ZeroMat(m);
}

SimpleMatrix::SimpleMatrix(const IdentityMat& m): SiconosMatrix(7, m.size1(), m.size2()), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Identity = new IdentityMat(m);
}

SimpleMatrix::SimpleMatrix(const std::string &file, bool ascii): SiconosMatrix(1), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Dense = new DenseMat();
  if (ascii)
  {
    ioMatrix io(file, "ascii");
    io.read(*this);
  }
  else
  {
    ioMatrix io(file, "binary");
    io.read(*this);
  }
  dimRow = (mat.Dense)->size1();
  dimCol = (mat.Dense)->size2();

}

SimpleMatrix::~SimpleMatrix()
{
  if (num == 1)
    delete(mat.Dense);
  else if (num == 2)
    delete(mat.Triang);
  else if (num == 3)
    delete(mat.Sym);
  else if (num == 4)
    delete(mat.Sparse);
  else if (num == 5)
    delete(mat.Banded);
  else if (num == 6)
    delete(mat.Zero);
  else if (num == 7)
    delete(mat.Identity);
  if (ipiv) delete ipiv;
}

//======================================
// get Ublas component (dense, sym ...)
//======================================

const DenseMat SimpleMatrix::getDense(unsigned int, unsigned int) const
{
  if (num != 1)
    SiconosMatrixException::selfThrow("SimpleMatrix::getDense(): the current matrix is not a Dense matrix");

  return *mat.Dense;
}

const TriangMat SimpleMatrix::getTriang(unsigned int, unsigned int) const
{
  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat SimpleMatrix::getTriang(): the current matrix is not a Triangular matrix");

  return *mat.Triang;
}

const SymMat SimpleMatrix::getSym(unsigned int, unsigned int) const
{
  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat SimpleMatrix::getSym(): the current matrix is not a Symmetric matrix");

  return *mat.Sym;
}

const SparseMat SimpleMatrix::getSparse(unsigned int, unsigned int) const
{
  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat SimpleMatrix::getSparse(): the current matrix is not a Sparse matrix");

  return *mat.Sparse;
}

const BandedMat SimpleMatrix::getBanded(unsigned int, unsigned int) const
{
  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat SimpleMatrix::getBanded(): the current matrix is not a Banded matrix");

  return *mat.Banded;
}

const ZeroMat SimpleMatrix::getZero(unsigned int, unsigned int) const
{
  if (num != 6)
    SiconosMatrixException::selfThrow("ZeroMat SimpleMatrix::getZero(): the current matrix is not a Zero matrix");

  return *mat.Zero;
}

const IdentityMat SimpleMatrix::getIdentity(unsigned int, unsigned int) const
{
  if (num != 7)
    SiconosMatrixException::selfThrow("IdentityMat SimpleMatrix::getIdentity(): the current matrix is not a Identity matrix");

  return *mat.Identity;
}

DenseMat* SimpleMatrix::getDensePtr(unsigned int, unsigned int) const
{
  if (num != 1)
    SiconosMatrixException::selfThrow("DenseMat* SimpleMatrix::getDensePtr(): the current matrix is not a Dense matrix");

  return mat.Dense;
}

TriangMat* SimpleMatrix::getTriangPtr(unsigned int, unsigned int) const
{
  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat* SimpleMatrix::getTriangPtr(): the current matrix is not a Triangular matrix");

  return mat.Triang;
}

SymMat* SimpleMatrix::getSymPtr(unsigned int, unsigned int) const
{
  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat* SimpleMatrix::getSymPtr(): the current matrix is not a Symmetric matrix");

  return mat.Sym;
}

SparseMat* SimpleMatrix::getSparsePtr(unsigned int, unsigned int) const
{
  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat* SimpleMatrix::getSparsePtr(): the current matrix is not a Sparse matrix");

  return mat.Sparse;
}

BandedMat* SimpleMatrix::getBandedPtr(unsigned int, unsigned int) const
{
  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat* SimpleMatrix::getBandedPtr(): the current matrix is not a Banded matrix");

  return mat.Banded;
}

ZeroMat* SimpleMatrix::getZeroPtr(unsigned int, unsigned int) const
{
  if (num != 6)
    SiconosMatrixException::selfThrow("ZeroMat* SimpleMatrix::getZeroPtr(): the current matrix is not a Zero matrix");

  return mat.Zero;
}

IdentityMat* SimpleMatrix::getIdentityPtr(unsigned int, unsigned int) const
{
  if (num != 7)
    SiconosMatrixException::selfThrow("IdentityMat* SimpleMatrix::getIdentityPtr(): the current matrix is not a Identity matrix");

  return mat.Identity;
}

double* SimpleMatrix::getArray(unsigned int, unsigned int) const
{
  if (num == 4)
    SiconosMatrixException::selfThrow("SimpleMatrix::getArray(): not yet implemented for sparse matrix.");

  if (num == 1)
    return &(((*mat.Dense).data())[0]);
  else if (num == 2)
    return &(((*mat.Triang).data())[0]);
  else if (num == 3)
    return &(((*mat.Sym).data())[0]);
  else if (num == 6)
  {
    ZeroMat::iterator1 it = (*mat.Zero).begin1();
    return const_cast<double*>(&(*it));
  }
  else if (num == 7)
  {
    IdentityMat::iterator1 it = (*mat.Identity).begin1();
    return const_cast<double*>(&(*it));
  }
  else
    return &(((*mat.Banded).data())[0]);
}

// ===========================
//       fill matrix
// ===========================

void SimpleMatrix::zero()
{
  unsigned int size1 = dimRow;
  unsigned int size2 = dimCol;
  if (num == 1)
    *mat.Dense = ublas::zero_matrix<double>(size1, size2);
  else if (num == 2)
    *mat.Triang = ublas::zero_matrix<double>(size1, size2);

  else if (num == 3)
    *mat.Sym = ublas::zero_matrix<double>(size1, size2);

  else if (num == 4)
    *mat.Sparse = ublas::zero_matrix<double>(size1, size2);

  else if (num == 5)
    *mat.Banded = ublas::zero_matrix<double>(size1, size2);

  else if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix::zero(): you can not set to zero a matrix of type Identity!.");
  resetLU();
  // if num == 6: nothing
}

void SimpleMatrix::eye()
{
  unsigned int size1 = dimRow;
  unsigned int size2 = dimCol;
  if (num == 1)
    *mat.Dense = ublas::identity_matrix<double>(size1, size2);

  else if (num == 2)
    *mat.Triang = ublas::identity_matrix<double>(size1, size2);

  else if (num == 3)
    *mat.Sym = ublas::identity_matrix<double>(size1, size2);

  else if (num == 4)
    *mat.Sparse = ublas::identity_matrix<double>(size1, size2);

  else if (num == 5)
    *mat.Banded = ublas::identity_matrix<double>(size1, size2);

  else if (num == 6)
    SiconosMatrixException::selfThrow("SimpleMatrix::eye(): you can not set to identity a matrix of type Zero!.");
  resetLU();
}

//=======================
// set matrix dimension
//=======================

void SimpleMatrix::resize(unsigned int row, unsigned int col, unsigned int lower, unsigned int upper, bool preserve)
{

  if (num == 1)
  {
    (*mat.Dense).resize(row, col, preserve);
  }
  else if (num == 2)
  {
    (*mat.Triang).resize(row, col, preserve);
  }
  else if (num == 3)
  {
    (*mat.Sym).resize(row, col, preserve);
  }
  else if (num == 4)
  {
    (*mat.Sparse).resize(row, col, preserve);
  }
  else if (num == 5)
  {
    (*mat.Banded).resize(row, col, lower, upper, preserve);
  }
  else if (num == 6)
  {
    (*mat.Zero).resize(row, col, preserve);
  }
  else if (num == 7)
  {
    (*mat.Identity).resize(row, col, preserve);
  }
  resetLU();
  dimRow = row;
  dimCol = col;
}

//=======================
//       get norm
//=======================

const double SimpleMatrix::normInf() const
{
  if (num == 1)
    return norm_inf(*mat.Dense);
  else if (num == 2)
    return norm_inf(*mat.Triang);
  else if (num == 3)
    return norm_inf(*mat.Sym);
  else if (num == 4)
    return norm_inf(*mat.Sparse);
  else if (num == 5)
    return norm_inf(*mat.Banded);
  else if (num == 6)
    return 0;
  else // if(num==7)
    return 1;
}

//=====================
// screen display
//=====================

void SimpleMatrix::display() const
{
  if (num == 1)
    std::cout << *mat.Dense << std::endl;
  else if (num == 2)
    std::cout << *mat.Triang << std::endl;
  else if (num == 3)
    std::cout << *mat.Sym << std::endl;
  else if (num == 4)
    std::cout << *mat.Sparse << std::endl;
  else if (num == 5)
    std::cout << *mat.Banded << std::endl;
  else if (num == 6)
    std::cout << *mat.Zero << std::endl;
  else if (num == 7)
    std::cout << *mat.Identity << std::endl;
}
//=============================
// Elements access (get or set)
//=============================

double& SimpleMatrix::operator()(unsigned int row, unsigned int col)
{
  if (row >= dimRow || col >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator(): Index out of range");

  if (num == 1)
    return (*mat.Dense)(row, col);
  else if (num == 2)
    return (*mat.Triang)(row, col);
  else if (num == 3)
    return (*mat.Sym)(row, col);
  else if (num == 4)
  {
    double *d = (*mat.Sparse).find_element(row, col);
    double & ref = *d;
    return ref;
  }
  else if (num == 5)
    return (*mat.Banded)(row, col);
  else if (num == 6)
    return const_cast<double&>((*mat.Zero)(row, col));
  else // i(num==7)
    return const_cast<double&>((*mat.Identity)(row, col));
}

const double SimpleMatrix::operator()(unsigned int row, unsigned int col) const
{
  if (row >= dimRow || col >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator(): Index out of range");

  if (num == 1)
    return (*mat.Dense)(row, col);
  else if (num == 2)
    return (*mat.Triang)(row, col);
  else if (num == 3)
    return (*mat.Sym)(row, col);
  else if (num == 4)
    return (*mat.Sparse)(row, col);
  else if (num == 5)
    return (*mat.Banded)(row, col);
  else if (num == 6)
    return 0.0;
  else // if (num == 7)
    return (row == col);
}

const double SimpleMatrix::getValue(unsigned int row, unsigned int col) const
{
  if (row >= dimRow || col >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix:getValue(index): Index out of range");

  if (num == 1)
    return (*mat.Dense)(row, col);
  else if (num == 2)
    return (*mat.Triang)(row, col);
  else if (num == 3)
    return (*mat.Sym)(row, col);
  else if (num == 4)
  {
    double * d = (*mat.Sparse).find_element(row, col);
    if (d)
      return *d;
    else
      return 0.0;
  }
  else if (num == 5)
    return (*mat.Banded)(row, col);
  else if (num == 6)
    return 0;
  else //if (num==7)
    return(row == col);
}

void SimpleMatrix::setValue(unsigned int row, unsigned int col, double value)
{
  if (row >= dimRow || col >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix:setValue: Index out of range");

  if (num == 1)
    (*mat.Dense)(row, col) = value;
  else if (num == 2)
    (*mat.Triang)(row, col) = value;
  else if (num == 3)
    (*mat.Sym)(row, col) = value ;
  else if (num == 4)
  {
    double * d = (*mat.Sparse).find_element(row, col);
    if (d)
    {
      *d = value;
    }
    else
    {
      (*mat.Sparse).insert_element(row, col, value);
    }
  }
  else if (num == 5)
    (*mat.Banded)(row, col) = value;
  else if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:setValue: forbidden for Identity or Zero type matrices.");
  resetLU();

}

//============================================
// Access (get or set) to blocks of elements
//============================================

void SimpleMatrix::setBlock(unsigned int row_min, unsigned int col_min, const SiconosMatrix& m)
{
  // Set current matrix elements, starting from row row_min and column col_min, with the values of the matrix m.
  // m may be a BlockMatrix.

  // Exceptions ...
  if (&m == this)
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(pos,..., m): m = this.");

  if (row_min >= dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(row,col): row is out of range");

  if (col_min >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::setBloc(row,col)k: col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if (row_max > dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(row,col,m): m.row + row is out of range.");

  if (col_max > dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(row,col,m): m.col + col is out of range.");

  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block matrix ...
  {
    BlocksMat::const_iterator1 it;
    BlocksMat::const_iterator2 it2;
    unsigned int posRow = row_min;
    unsigned int posCol = col_min;

    for (it = m.begin(); it != m.end(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        setBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = col_min;
    }
  }
  else // if m is a SimpleMatrix
  {
    if (numM != num)
      SiconosMatrixException::selfThrow("SimpleMatrix::setBlock(i,j,m), inconsistent types.");

    if (num == 1)
      noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) = *(m.getDensePtr());
    else if (num == 2)
      noalias(ublas::subrange(*mat.Triang, row_min, row_max, col_min, col_max)) = *(m.getTriangPtr());
    else if (num == 3)
      noalias(ublas::subrange(*mat.Sym, row_min, row_max, col_min, col_max)) = *(m.getSymPtr());
    else if (num == 4)
      noalias(ublas::subrange(*mat.Sparse, row_min, row_max, col_min, col_max)) = *(m.getSparsePtr());
    else if (num == 5)
      noalias(ublas::subrange(*mat.Banded, row_min, row_max, col_min, col_max)) = *(m.getBandedPtr());
    else // if(num==6) or num == 7 nothing to do
    {}
    resetLU();
  }
}

void setBlock(SPC::SiconosMatrix  MIn, SP::SiconosMatrix MOut, const Index& dim, const Index& start)
{
  // To copy a subBlock of MIn into a subBlock of MOut.
  // dim[0], dim[1]: number of rows and columns of the sub-block
  // start[0], start[1]: position (row, column) of the first element of the subBlock in MIn
  // start[2], start[3]: position (row, column) of the first element of the subBlock in MOut

  if (MIn == MOut) // useless op => nothing to be done.
  {}// SiconosVectorException::selfThrow("");
  else
  {
    unsigned int numIn = MIn->getNum();
    unsigned int numOut = MOut->getNum();

    if (numOut == 6 || numOut == 7) // if MOut = 0 or Identity => read-only
      SiconosMatrixException::selfThrow("matrices, setBlock(MIn, MOut...): MOut is read-only (zero or identity matrix?).");

    // Check dimension
    Index MDim(4); // dim. of matrices MIn and MOut.
    MDim[0] = MIn->size(0);
    MDim[1] = MIn->size(1);
    MDim[2] = MOut->size(0);
    MDim[3] = MOut->size(1);

    for (unsigned int i = 0; i < 4 ; ++i)
      if (start[i] >= MDim[i])
        SiconosMatrixException::selfThrow("matrices, setBlock(MIn, ...): sub-block indices are out of range.");

    // index position of the last element in subBlock ...
    Index end(4);
    end[0] = dim[0] + start[0];
    end[1] = dim[1] + start[1];
    end[2] = dim[0] + start[2];
    end[3] = dim[1] + start[3];

    for (unsigned int i = 0; i < 4 ; ++i)
      if (end[i] > MDim[i])
        SiconosMatrixException::selfThrow("matrices, setBlock(MIn, ...): sub-block indices are out of range.");

    // Elements from row/col start[i] to row/col (end[i]-1) will be copied.

    // If both matrices MIn and MOut are block, exception.
    if (numIn == 0 && numOut == 0)
      SiconosMatrixException::selfThrow("matrices, setBlock(MIn, MOut ...): not yet implemented for MIn and MOut both BlockMatrix. Try to use setBlock on the sub-matrices?");

    else if (numOut == 0) // if MOut is a BlockMatrix.
    {

      // Steps:
      // A - Find the blocks of MOut that "own" indices start[2] and end[2] ie
      //     the first and last sub-block to be set in a block-column
      //         --> numbers blockStart0 and blockEnd0
      // B - Find the  Block of MOut that "owns" index start[3] and end[3] ie
      //     the first sub-block to be set in a block-row
      //         --> numbers blockStart1 and blockEnd1
      //
      //        => The blocks concerned in MOut, are those between (block) rows blockStart0 and blockEnd0
      //           and (block) columns blockStart1 and blockEnd1.
      //
      // C - Loop through the concerned blocks (name = currentBlock) of MOut and call setBlock(MIn, currentBlock, subSize, currentPos).
      //     subSize: dim. of the considered sub-block of currentBlock to be set
      //     currentPos: same as "start" vector but for currentBlock
      //

      // A - Block-Row position: we look for the block of MOut that include index start[2] and end[2].
      //
      unsigned int blockStart0 = 0;
      const Index * tab = MOut->getTabRowPtr();
      while (start[2] >= (*tab)[blockStart0] && blockStart0 < tab->size())
        blockStart0 ++;
      // Relative position in the block blockStart0 of the first element to be set.
      unsigned int posOut0 = start[2];
      if (blockStart0 != 0)
        posOut0 -= (*tab)[blockStart0 - 1];

      unsigned int blockEnd0 = blockStart0;
      while (end[2] > (*tab)[blockEnd0] && blockEnd0 < tab->size())
        blockEnd0 ++;

      // Size of the last sub-block in the column of block
      unsigned int lastBlockSize0 = end[2];
      if (blockEnd0 != 0)
        lastBlockSize0 -= (*tab)[blockEnd0 - 1];

      // B - Block-Col position: we look for the block of MOut that include index start[3] and end[3].
      unsigned int blockStart1 = 0;
      tab = MOut->getTabColPtr();
      while (start[3] >= (*tab)[blockStart1] && blockStart1 < tab->size())
        blockStart1 ++;
      // Relative position in the block blockStart1 of the first element to be set.
      unsigned int posOut1 = start[3];
      if (blockStart1 != 0)
        posOut1 -= (*tab)[blockStart1 - 1];

      unsigned int blockEnd1 = blockStart1;
      while (end[3] > (*tab)[blockEnd1] && blockEnd1 < tab->size())
        blockEnd1 ++;

      // Size of the last sub-block in the row of block
      unsigned int lastBlockSize1 = end[3];
      if (blockEnd1 != 0)
        lastBlockSize1 -= (*tab)[blockEnd1 - 1];

      //C - Next, 3 steps for each row:
      // - set first sub-block in the row (number blockStart1)
      // - set all other blocks in the row except the last one
      // - set last block (number blockEnd1)
      // Same process for other rows ...

      // The current considered block
      SP::SiconosMatrix   currentBlock = MOut->getBlockPtr(blockStart0, blockStart1);

      // dim of the subBlock of currentBlock to be set.
      Index subSize(2);
      // indices of the first element of MIn (resp. currentBlock) to be read (resp. set)  (same as start for MIn and MOut).
      Index currentPos(4);

      // currentBlock position in MOut.
      unsigned int numRow = blockStart0;
      unsigned int numCol = blockStart1;

      // Init currentPos
      // row and col position for first element to be read in MIn,
      currentPos[0] = start[0];
      currentPos[1] = start[1];
      // row and col position for first element in sub-block of Mout (namely currentBlock).
      currentPos[2] = posOut0;
      currentPos[3] = posOut1;

      while (numRow != blockEnd0 + 1)
      {

        while (numCol != blockEnd1 + 1)
        {
          // Get the block of MOut from which a sub-block will be set ...
          currentBlock = MOut->getBlockPtr(numRow, numCol);

          // Set subSize[0], dim (rows) and subSize[1], dim (columns) of the sub-block.
          // subSize[0] is only required for the first block in the row, after it remains constant.
          subSize[1] = currentBlock->size(1);

          // Warning: test "a" must be done before test "b"
          if (numCol == blockEnd1) // if last column of blocks -> test "a"
            subSize[1] = lastBlockSize1;

          if (numCol == blockStart1) // -> test "b"
          {
            subSize[1] -= posOut1;
            subSize[0] = currentBlock->size(0);
            if (numRow == blockEnd0) // if last row of blocks
              subSize[0] = lastBlockSize0;
            if (numRow == blockStart0) // if first row of blocks
              subSize[0] -= posOut0;
          }

          // Set sub-block
          setBlock(MIn, currentBlock, subSize, currentPos);

          // Update currentPos:
          // col position for first element to be read in MIn,
          currentPos[1] += subSize[1] ;
          // col position for first element to be set in sub-block.
          currentPos[3] = 0;
          numCol++;
        }

        numCol = blockStart1;
        numRow++;

        // Update currentPos:
        // row position for first element to be read in MIn,
        currentPos[0] += subSize[0] ;
        // col position for first element to be read in MIn,
        currentPos[1] = start[1] ;
        // row position for first element to be set in sub-block.
        currentPos[2] = 0;
        // col position for first element to be set in sub-block.
        currentPos[3] = posOut1;

      }

    }
    else if (numIn == 0) // If MIn is a BlockMatrix.
    {

      // Same process as for numOut == 0

      unsigned int blockStart0 = 0;
      const Index * tab = MIn->getTabRowPtr();
      while (start[0] >= (*tab)[blockStart0] && blockStart0 < tab->size())
        blockStart0 ++;
      // Relative position in the block blockStart0 of the first element to be set.
      unsigned int posOut0 = start[0];
      if (blockStart0 != 0)
        posOut0 -= (*tab)[blockStart0 - 1];

      unsigned int blockEnd0 = blockStart0;
      while (end[0] > (*tab)[blockEnd0] && blockEnd0 < tab->size())
        blockEnd0 ++;

      // Size of the last sub-block in the column of block
      unsigned int lastBlockSize0 = end[0];
      if (blockEnd0 != 0)
        lastBlockSize0 -= (*tab)[blockEnd0 - 1];

      // B - Block-Col position: we look for the block of MOut that include index start[3] and end[3].
      unsigned int blockStart1 = 0;
      tab = MIn->getTabColPtr();
      while (start[1] >= (*tab)[blockStart1] && blockStart1 < tab->size())
        blockStart1 ++;
      // Relative position in the block blockStart1 of the first element to be set.
      unsigned int posOut1 = start[1];
      if (blockStart1 != 0)
        posOut1 -= (*tab)[blockStart1 - 1];

      unsigned int blockEnd1 = blockStart1;
      while (end[1] > (*tab)[blockEnd1] && blockEnd1 < tab->size())
        blockEnd1 ++;

      // Size of the last sub-block in the row of block
      unsigned int lastBlockSize1 = end[1];
      if (blockEnd1 != 0)
        lastBlockSize1 -= (*tab)[blockEnd1 - 1];

      //C - Next, 3 steps for each row:
      // - set first sub-block in the row (number blockStart1)
      // - set all other blocks in the row except the last one
      // - set last block (number blockEnd1)
      // Same process for other rows ...

      // The current considered block
      SPC::SiconosMatrix  currentBlock = MIn->getBlockPtr(blockStart0, blockStart1);

      // dim of the subBlock of currentBlock to be set.
      Index subSize(2);
      // indices of the first element of MIn (resp. currentBlock) to be read (resp. set)  (same as start for MIn and MOut).
      Index currentPos(4);

      // currentBlock position in MOut.
      unsigned int numRow = blockStart0;
      unsigned int numCol = blockStart1;

      // Init currentPos
      // row and col position for first element to be read in MIn,
      currentPos[0] = posOut0;
      currentPos[1] = posOut1;
      // row and col position for first element in sub-block of Mout (namely currentBlock).
      currentPos[2] = start[2];
      currentPos[3] = start[3];

      while (numRow != blockEnd0 + 1)
      {

        while (numCol != blockEnd1 + 1)
        {
          // Get the block of MOut from which a sub-block will be set ...
          currentBlock = MIn->getBlockPtr(numRow, numCol);

          // Set subSize[0], dim (rows) and subSize[1], dim (columns) of the sub-block.
          // subSize[0] is only required for the first block in the row, after it remains constant.
          subSize[1] = currentBlock->size(1);
          // Warning: test "a" must be done before test "b"
          if (numCol == blockEnd1) // if last column of blocks -> test "a"
            subSize[1] = lastBlockSize1;

          if (numCol == blockStart1) // -> test "b"
          {
            subSize[1] -= posOut1;
            subSize[0] = currentBlock->size(0);
            if (numRow == blockEnd0) // if last row of blocks
              subSize[0] = lastBlockSize0;
            if (numRow == blockStart0) // if first row of blocks
              subSize[0] -= posOut0;
          }

          // Set sub-block
          setBlock(currentBlock, MOut, subSize, currentPos);

          // Update currentPos:
          // col position for first element to be read in MIn,
          currentPos[1] = 0 ;
          // col position for first element to be set in sub-block.
          currentPos[3] += subSize[1];
          numCol++;
        }

        numCol = blockStart1;
        numRow++;

        // Update currentPos:
        // row position for first element to be read in MIn,
        currentPos[0] = 0;
        // col position for first element to be read in MIn,
        currentPos[1] = posOut1;
        // row position for first element to be set in sub-block.
        currentPos[2] += subSize[0] ;
        // col position for first element to be set in sub-block.
        currentPos[3] = start[3];

      }
      MOut->resetLU();

    }
    else // neither MIn nor MOut is a BlockMatrix.
    {
      switch (numIn)
      {
      case 1:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->getDensePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getDensePtr(), start[0], end[0], start[1], end[1]);
        break;

      case 2:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->getDensePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getTriangPtr(), start[0], end[0], start[1], end[1]);
        break;

      case 3:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->getDensePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getSymPtr(), start[0], end[0], start[1], end[1]);
        break;

      case 4:
        if (numOut == 1)
          noalias(ublas::subrange(*MOut->getDensePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getSparsePtr(), start[0], end[0], start[1], end[1]);
        else if (numOut == 4)
          noalias(ublas::subrange(*MOut->getSparsePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getSparsePtr(), start[0], end[0], start[1], end[1]);
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      case 5:
        if (numOut != 1)
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        noalias(ublas::subrange(*MOut->getDensePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getBandedPtr(), start[0], end[0], start[1], end[1]);
        break;

      case 6:
        if (numOut == 1)
          ublas::subrange(*MOut->getDensePtr(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if (numOut == 2)
          ublas::subrange(*MOut->getTriangPtr(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if (numOut == 4)
          ublas::subrange(*MOut->getSparsePtr(), start[2], end[2], start[3], end[3]) *= 0.0;
        else if (numOut == 5)
          ublas::subrange(*MOut->getBandedPtr(), start[2], end[2], start[3], end[3]) *= 0.0;
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      case 7:
        if (numOut == 1)
          noalias(ublas::subrange(*MOut->getDensePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getIdentityPtr(), start[0], end[0], start[1], end[1]);
        else if (numOut == 4)
          noalias(ublas::subrange(*MOut->getSparsePtr(), start[2], end[2], start[3], end[3])) = ublas::subrange(*MIn->getIdentityPtr(), start[0], end[0], start[1], end[1]);
        else
          SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;

      default:
        SiconosMatrixException::selfThrow("matrix, setBlock(MIn, MOut, ...), unconsistent types between MIn and MOut.");
        break;
      }
      MOut->resetLU();
    }
  }
}


void SimpleMatrix::getRow(unsigned int r, SiconosVector &vOut) const
{
  // Get row number r of current matrix and copy it into vOut.
  if (r >= dimRow)
    SiconosMatrixException::selfThrow("getRow(row): row is out of range");

  if (vOut.size() != dimCol)
    SiconosMatrixException::selfThrow("getRow(row,v): inconsistent sizes between this and v.");

  if (num == 7) // identity matrix
  {
    vOut.zero();
    vOut(r) = 1.0;
  }
  else if (num == 6) // Zero matrix
    vOut.zero();
  else
  {
    unsigned int numV = vOut.getNum();
    unsigned int pos = 0;
    if (numV == 0) // vOut is Block
    {
      VectorOfVectors::iterator it;
      for (it = vOut.begin(); it != vOut.end(); ++it)
      {
        getSubRow(r, pos, *it);
        pos += (*it)->size();
      }
    }
    else if (numV == 1)
    {

      if (num == 1)
      {
        noalias(*(vOut.getDensePtr())) = ublas::row(*mat.Dense, r);
      }
      else if (num == 2)
      {
        noalias(*(vOut.getDensePtr())) = ublas::row(*mat.Triang, r);
      }
      else if (num == 3)
      {
        noalias(*(vOut.getDensePtr())) = ublas::row(*mat.Sym, r);
      }
      else if (num == 4)
      {
        noalias(*(vOut.getDensePtr())) = ublas::row(*mat.Sparse, r);
      }
      else //if(num==5){
        noalias(*(vOut.getDensePtr())) = ublas::row(*mat.Banded, r);
    }
    else // if numV == 4
    {
      if (num == 4)
      {
        noalias(*(vOut.getSparsePtr())) = ublas::row(*mat.Sparse, r);
      }
      else
        SiconosMatrixException::selfThrow("getRow(row,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void SimpleMatrix::setRow(unsigned int r, const SiconosVector& vIn)
{
  // Set row number r of current matrix with vIn.
  unsigned int numV = vIn.getNum();
  if (r >= dimRow)
    SiconosMatrixException::selfThrow("setRow(row): row is out of range");

  if (vIn.size() != dimCol)
    SiconosMatrixException::selfThrow("setRow(row,v): inconsistent sizes between this and v.");

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("setRow(row,v): current matrix is read-only (zero or identity).");

  if (numV == 0) // vIn Block
  {
    VectorOfVectors::const_iterator it;
    unsigned int pos = 0;

    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      setSubRow(r, pos, *it);
      pos += (*it)->size();
    }
  }
  else
  {
    if (num == 1)
    {
      if (numV == 1)
      {
        noalias(ublas::row(*mat.Dense, r)) = *vIn.getDensePtr();
      }
      else if (numV == 4)
      {
        noalias(ublas::row(*mat.Dense, r)) = *vIn.getSparsePtr();
      }
    }
    else if (num == 4 && numV == 4)
      noalias(ublas::row(*mat.Sparse, r)) = *vIn.getSparsePtr();
    else
      SiconosMatrixException::selfThrow("setRow(row,v): inconsistent types between current matrix and v.");
  }

  resetLU();
}

void SimpleMatrix::getCol(unsigned int r, SiconosVector &vOut)const
{
  // Get column number r of current matrix and copy it into vOut.
  if (r >= dimCol)
    SiconosMatrixException::selfThrow("getCol(col): col is out of range");

  if (vOut.size() != dimRow)
    SiconosMatrixException::selfThrow("getCol(col,v): inconsistent sizes between this and v.");

  if (num == 7) // identity matrix
  {
    vOut.zero();
    vOut(r) = 1.0;
  }
  else if (num == 6) // Zero matrix
    vOut.zero();
  else
  {
    unsigned int numV = vOut.getNum();

    if (numV == 0) // vOut is Block
    {
      VectorOfVectors::iterator it;
      unsigned int pos = 0;
      for (it = vOut.begin(); it != vOut.end(); ++it)
      {
        getSubCol(r, pos, *it);
        pos += (*it)->size();
      }
    }
    else if (numV == 1)
    {

      if (num == 1)
      {
        noalias(*(vOut.getDensePtr())) = ublas::column(*mat.Dense, r);
      }
      else if (num == 2)
      {
        noalias(*(vOut.getDensePtr())) = ublas::column(*mat.Triang, r);
      }
      else if (num == 3)
      {
        noalias(*(vOut.getDensePtr())) = ublas::column(*mat.Sym, r);
      }
      else if (num == 4)
      {
        noalias(*(vOut.getDensePtr())) = ublas::column(*mat.Sparse, r);
      }
      else //if(num==5){
        noalias(*(vOut.getDensePtr())) = ublas::column(*mat.Banded, r);
    }
    else // if numV == 4
    {
      if (num == 4)
      {
        noalias(*(vOut.getSparsePtr())) = ublas::column(*mat.Sparse, r);
      }
      else
        SiconosMatrixException::selfThrow("getCol(col,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }
}

void SimpleMatrix::setCol(unsigned int r, const SiconosVector &vIn)
{
  // Set column number r of current matrix with vIn.
  unsigned int numV = vIn.getNum();
  if (r >= dimCol)
    SiconosMatrixException::selfThrow("setCol(col): col is out of range");

  if (vIn.size() != dimRow)
    SiconosMatrixException::selfThrow("setCol(col,v): inconsistent sizes between this and v.");

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("setCol(col,v): current matrix is read-only (zero or identity).");

  if (numV == 0) // vIn Block
  {
    VectorOfVectors::const_iterator it;
    unsigned int pos = 0;
    for (it = vIn.begin(); it != vIn.end(); ++it)
    {
      setSubCol(r, pos, *it);
      pos += (*it)->size();
    }
  }
  else
  {
    if (num == 1)
    {
      if (numV == 1)
      {
        noalias(ublas::column(*mat.Dense, r)) = *vIn.getDensePtr();
      }
      else if (numV == 4)
      {
        noalias(ublas::column(*mat.Dense, r)) = *vIn.getSparsePtr();
      }
    }
    else if (num == 4 && numV == 4)
      noalias(ublas::column(*mat.Sparse, r)) = *vIn.getSparsePtr();
    else
      SiconosMatrixException::selfThrow("setCol(col,v): inconsistent types between current matrix and v.");
  }

  resetLU();
}

void SimpleMatrix::getSubRow(unsigned int r, unsigned int pos, SP::SiconosVector vOut) const
{
  // Get row number r of current matrix, starting from element at position pos, and copy it into vOut.
  if (r >= dimRow)
    SiconosMatrixException::selfThrow("getSubRow(row,pos,v): row is out of range");

  if (vOut->size() > dimCol - pos)
    SiconosMatrixException::selfThrow("getSubRow(row,pos,v): inconsistent sizes between this and v.");

  if (num == 7) // identity matrix
  {
    vOut->zero();
    if (r - pos >= 0)
      (*vOut)(r - pos) = 1.0;
  }
  else if (num == 6) // Zero matrix
    vOut->zero();
  else
  {
    unsigned int numV = vOut->getNum();
    unsigned int subPos = pos;
    unsigned int nbEl = vOut->size();
    if (numV == 0) // vOut is Block
    {
      VectorOfVectors::iterator it;
      for (it = vOut->begin(); it != vOut->end(); ++it)
      {
        getSubRow(r, subPos, *it);
        subPos += (*it)->size();
      }
    }
    else if (numV == 1)
    {
      if (num == 1)
      {
        //      noalias(*(vOut->getDensePtr())) = ublas::row(ublas::subrange(*mat.Dense, r, r+1,pos, endPos),0);
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if (num == 2)
      {
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<TriangMat >(*mat.Triang, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if (num == 3)
      {
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<SymMat >(*mat.Sym, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
      }
      else if (num == 4)
      {
#ifdef BOOST_LIMITATION
        SiconosMatrixException("SimpleMatrix::getSubRow warning - ublas::matrix_vector_slice<SparseMat> does not exist for MacOS.");
#else
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
#endif
      }
      else //if(num==5){
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<BandedMat >(*mat.Banded, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
    }
    else // if numV == 4
    {
      if (num == 4)
      {
#ifdef BOOST_LIMITATION
        SiconosMatrixException("SimpleMatrix::getSubRow warning - ublas::matrix_vector_slice<SparseMat> does not exist for MacOs.");
#else
        noalias(*(vOut->getSparsePtr())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl));
#endif
      }
      else
        SiconosMatrixException::selfThrow("getSubRow(row,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }

}

void SimpleMatrix::setSubRow(unsigned int r, unsigned int pos, SP::SiconosVector vIn)
{
  // Set row number r, starting from element at position pos, of current matrix with vIn.
  unsigned int numV = vIn->getNum();
  if (r >= dimRow)
    SiconosMatrixException::selfThrow("setSubRow(row): row is out of range");

  if (vIn->size() > dimCol - pos)
    SiconosMatrixException::selfThrow("setSubRow(row,v): inconsistent sizes between this and v.");

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("setSubRow(row,v): current matrix is read-only (zero or identity).");

  if (numV == 0) // vIn Block
  {
    VectorOfVectors::const_iterator it;
    unsigned int subPos = pos;
    for (it = vIn->begin(); it != vIn->end(); ++it)
    {
      setSubRow(r, subPos, *it);
      subPos += (*it)->size();
    }
  }
  else
  {
    unsigned int nbEl = vIn->size();
    if (num == 1)
    {
      if (numV == 1)
      {
        noalias(ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl))) = *vIn->getDensePtr();
      }
      else if (numV == 4)
      {
        ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl)) = *vIn->getSparsePtr();
      }
    }
    else if (num == 4 && numV == 4)
#ifdef BOOST_LIMITATION
      SiconosMatrixException("SimpleMatrix::setSubRow warning - ublas::matrix_vector_slice<SparseMat> does not exist for MacOs.");
#else
      ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(r, 0, nbEl), ublas::slice(pos, 1, nbEl)) = *vIn->getSparsePtr();
#endif
    else
      SiconosMatrixException::selfThrow("setSubRow(row,v): inconsistent types between current matrix and v.");
    resetLU();
  }

}

void SimpleMatrix::getSubCol(unsigned int r, unsigned int pos, SP::SiconosVector vOut) const
{
  // Get col number r of current matrix, starting from element at position pos, and copy it into vOut.
  if (r >= dimCol)
    SiconosMatrixException::selfThrow("getSubCol(col,pos,v): col is out of range");

  if (vOut->size() > dimRow - pos)
    SiconosMatrixException::selfThrow("getSubCol(col,pos,v): inconsistent sizes between this and v.");

  if (num == 7) // identity matrix
  {
    vOut->zero();
    if (r - pos >= 0)
      (*vOut)(r - pos) = 1.0;
  }
  else if (num == 6) // Zero matrix
    vOut->zero();
  else
  {
    unsigned int numV = vOut->getNum();
    unsigned int subPos = pos;
    unsigned int nbEl = vOut->size();
    if (numV == 0) // vOut is Block
    {
      VectorOfVectors::iterator it;
      for (it = vOut->begin(); it != vOut->end(); ++it)
      {
        getSubRow(r, subPos, *it);
        subPos += (*it)->size();
      }
    }

    else if (numV == 1)
    {
      if (num == 1)
      {
        //      noalias(*(vOut->getDensePtr())) = ublas::row(ublas::subrange(*mat.Dense, r, r+1,pos, endPos),0);
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if (num == 2)
      {
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<TriangMat >(*mat.Triang, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if (num == 3)
      {
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<SymMat >(*mat.Sym, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
      }
      else if (num == 4)
      {
#ifdef BOOST_LIMITATION
        SiconosMatrixException("SimpleMatrix::getSubCol warning - ublas::matrix_vector_slice<SparseMat> does not exist for MacOs.");
#else
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
#endif
      }
      else //if(num==5){
        noalias(*(vOut->getDensePtr())) = ublas::matrix_vector_slice<BandedMat >(*mat.Banded, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
    }
    else // if numV == 4
    {
      if (num == 4)
      {
#ifdef BOOST_LIMITATION
        SiconosMatrixException("SimpleMatrix::getSubCol warning - ublas::matrix_vector_slice<SparseMat> does not exist for MacOs.");
#else
        noalias(*(vOut->getSparsePtr())) = ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl));
#endif
      }
      else
        SiconosMatrixException::selfThrow("getSubCol(col,v): inconsistent types between this (not sparse) and v (sparse).");
    }
  }

}

void SimpleMatrix::setSubCol(unsigned int r, unsigned int pos, SP::SiconosVector vIn)
{
  // Set column number r, starting from element at position pos, of current matrix with vIn.
  unsigned int numV = vIn->getNum();
  if (r >= dimCol)
    SiconosMatrixException::selfThrow("setSubCol(col): col is out of range");

  if (vIn->size() > dimRow - pos)
    SiconosMatrixException::selfThrow("setSubCol(col,v): inconsistent sizes between this and v.");

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("setSubCol(col,v): current matrix is read-only (zero or identity).");

  if (numV == 0) // vIn Block
  {
    VectorOfVectors::const_iterator it;
    unsigned int subPos = pos;
    for (it = vIn->begin(); it != vIn->end(); ++it)
    {
      setSubCol(r, subPos, *it);
      subPos += (*it)->size();
    }
  }
  else
  {
    unsigned int nbEl = vIn->size();
    if (num == 1)
    {
      if (numV == 1)
      {
        noalias(ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl))) = *vIn->getDensePtr();
      }
      else if (numV == 4)
      {
        ublas::matrix_vector_slice<DenseMat >(*mat.Dense, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl)) = *vIn->getSparsePtr();
      }
    }
    else if (num == 4 && numV == 4)
#ifdef BOOST_LIMITATION
      SiconosMatrixException("SimpleMatrix::setSubCol warning - ublas::matrix_vector_slice<SparseMat> does not exist for MacOs.");
#else
      ublas::matrix_vector_slice<SparseMat >(*mat.Sparse, ublas::slice(pos, 1, nbEl), ublas::slice(r, 0, nbEl)) = *vIn->getSparsePtr();
#endif
    else
      SiconosMatrixException::selfThrow("setSubCol(row,v): inconsistent types between current matrix and v.");
    resetLU();
  }
}



void SimpleMatrix::addBlock(unsigned int row_min, unsigned int col_min, const SiconosMatrix& m)
{
  // add m to current matrix elements, starting from row row_min and column col_min, to the values of the matrix m.
  // m may be a BlockMatrix.

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(pos,..., m) forbidden for zero or identity matrix.");

  if (&m == this)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(pos,..., m): m = this.");

  if (row_min >= dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(row,col): row is out of range");

  if (col_min >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBloc(row,col)k: col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if (row_max > dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(row,col,m): m.row + row is out of range.");

  if (col_max > dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(row,col,m): m.col + col is out of range.");

  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block matrix ...
  {
    BlocksMat::const_iterator1 it;
    BlocksMat::const_iterator2 it2;
    unsigned int posRow = row_min;
    unsigned int posCol = col_min;

    for (it = m.begin(); it != m.end(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        addBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if (numM == 6) // if m = 0
  {
    // nothing to do !
  }
  else // if m is a SimpleMatrix
  {
    if (num == 1)
    {
      switch (numM)
      {
      case 1:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.getDensePtr());
        break;
      case 2:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.getTriangPtr());
        break;
      case 3:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.getSymPtr());
        break;
      case 4:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.getSparsePtr());
        break;
      case 5:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.getBandedPtr());
        break;
      case 7:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) += *(m.getIdentityPtr());
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(...,m): wrong matrix type for m.");
        break;
      }
    }
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::addBlock(...): implemeted only for dense matrices.");
    resetLU();
  }
}

void SimpleMatrix::subBlock(unsigned int row_min, unsigned int col_min, const SiconosMatrix& m)
{
  // sub m to current matrix elements, starting from row row_min and column col_min, to the values of the matrix m.
  // m may be a BlockMatrix.

  if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(pos,..., m) forbidden for zero or identity matrix.");

  if (&m == this)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(pos,..., m): m = this.");

  if (row_min >= dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col): row is out of range");

  if (col_min >= dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col): col is out of range");

  unsigned int row_max, col_max;
  row_max = m.size(0) + row_min;
  col_max = m.size(1) + col_min;

  if (row_max > dimRow)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col,m): m.row + row is out of range.");

  if (col_max > dimCol)
    SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(row,col,m): m.col + col is out of range.");

  unsigned int numM = m.getNum();

  if (numM == 0) // if m is a block matrix ...
  {
    BlocksMat::const_iterator1 it;
    BlocksMat::const_iterator2 it2;
    unsigned int posRow = row_min;
    unsigned int posCol = col_min;

    for (it = m.begin(); it != m.end(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        subBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if (numM == 6) // if m = 0
  {
    // nothing to do !
  }
  else // if m is a SimpleMatrix
  {
    if (num == 1)
    {
      switch (numM)
      {
      case 1:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.getDensePtr());
        break;
      case 2:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.getTriangPtr());
        break;
      case 3:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.getSymPtr());
        break;
      case 4:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.getSparsePtr());
        break;
      case 5:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.getBandedPtr());
        break;
      case 7:
        noalias(ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max)) -= *(m.getIdentityPtr());
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(...,m): wrong matrix type for m.");
        break;
      }
    }
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::subBlock(...): implemeted only for dense matrices.");
    resetLU();
  }
}


//=============
// Assignment
//=============

SimpleMatrix& SimpleMatrix::operator = (const SiconosMatrix& m)
{

  if (&m == this) return *this; // auto-assignment.

  unsigned int numM = m.getNum();

  if (dimRow != m.size(0) || dimCol != m.size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix::operator = failed. Inconsistent sizes.");

  if (numM == 6) // m = zero matrix
  {
    zero();
    return *this;
  }

  if (numM == 7) // m = identity matrix
  {
    eye();
    return *this;
  }

  if (numM == 0) // if m is a BlockMatrix
  {
    ConstBlockIterator1 it;
    ConstBlockIterator2 it2;
    unsigned int posRow = 0;
    unsigned int posCol = 0;

    for (it = m.begin(); it != m.end(); ++it)
    {
      for (it2 = it.begin(); it2 != it.end(); ++it2)
      {
        setBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else
  {
    switch (num)
    {
    case 1:
      switch (numM)
      {
      case 1:
        noalias(*(mat.Dense)) = *m.getDensePtr();
        break;
      case 2:
        noalias(*(mat.Dense)) = *m.getTriangPtr();
        break;
      case 3:
        noalias(*(mat.Dense)) = *m.getSymPtr();
        break;
      case 4:
        noalias(*(mat.Dense)) = *m.getSparsePtr();
        break;
      case 5:
        noalias(*(mat.Dense)) = *m.getBandedPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
        break;
      }
      break;
    case 2:
      switch (numM)
      {
      case 2:
        noalias(*(mat.Triang)) = *m.getTriangPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::assignment of a bad type of matrix into a triangular one.");
        break;
      }
      break;
    case 3:
      if (numM == 3)
        noalias(*(mat.Sym)) = *m.getSymPtr();
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::bad assignment of matrix (symetric one = dense or ...)");
      break;
    case 4:
      switch (numM)
      {
      case 2:
        noalias(*(mat.Sparse)) = *m.getTriangPtr();
        break;
      case 3:
        noalias(*(mat.Sparse)) = *m.getSymPtr();
        break;
      case 4:
        noalias(*(mat.Sparse)) = *m.getSparsePtr();
        break;
      case 5:
        noalias(*(mat.Sparse)) = *m.getBandedPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
        break;
      }
      break;
    case 5:
      switch (numM)
      {
      case 5:
        noalias(*(mat.Banded)) = *m.getBandedPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
        break;
      }
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
      break;
    }
    resetLU();
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator = (const SimpleMatrix& m)
{

  if (&m == this) return *this; // auto-assignment.

  unsigned int numM = m.getNum();

  if (dimRow != m.size(0) || dimCol != m.size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix::operator = failed. Inconsistent sizes.");

  if (numM == 6) // m = zero matrix
  {
    zero();
    return *this;
  }
  else if (numM == 7) // m = identity matrix
  {
    eye();
    return *this;
  }

  switch (num)
  {
  case 1:
    switch (numM)
    {
    case 1:
      noalias(*(mat.Dense)) = *m.getDensePtr();
      break;
    case 2:
      noalias(*(mat.Dense)) = *m.getTriangPtr();
      break;
    case 3:
      noalias(*(mat.Dense)) = *m.getSymPtr();
      break;
    case 4:
      noalias(*(mat.Dense)) = *m.getSparsePtr();
      break;
    case 5:
      noalias(*(mat.Dense)) = *m.getBandedPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (numM)
    {
    case 2:
      noalias(*(mat.Triang)) = *m.getTriangPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::assignment of a bad type of matrix into a triangular one.");
      break;
    }
    break;
  case 3:
    if (numM == 3)
      noalias(*(mat.Sym)) = *m.getSymPtr();
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::bad assignment of matrix (symetric one = dense or ...)");
    break;
  case 4:
    switch (numM)
    {
    case 2:
      noalias(*(mat.Sparse)) = *m.getTriangPtr();
      break;
    case 3:
      noalias(*(mat.Sparse)) = *m.getSymPtr();
      break;
    case 4:
      noalias(*(mat.Sparse)) = *m.getSparsePtr();
      break;
    case 5:
      noalias(*(mat.Sparse)) = *m.getBandedPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (numM)
    {
    case 5:
      noalias(*(mat.Banded)) = *m.getBandedPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix): invalid type of matrix");
    break;
    resetLU();
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator = (const DenseMat& m)
{
  if (num != 1)
    SiconosVectorException::selfThrow("SimpleMatrix::operator = DenseMat : forbidden: the current matrix is not dense.");

  if (dimRow != m.size1() || dimCol != m.size2())
    SiconosMatrixException::selfThrow("SimpleMatrix::operator = DenseMat failed. Inconsistent sizes.");

  noalias(*(mat.Dense)) = m;
  //atlas::copy(m, *(mat.Dense));

  resetLU();
  return *this;
}

//=================================
// Op. and assignment (+=, -= ... )
//=================================

SimpleMatrix& SimpleMatrix::operator +=(const SiconosMatrix& m)
{

  unsigned int numM = m.getNum();
  if (numM == 6) // m = 0
    return *this;

  if (&m == this) // auto-assignment
  {
    switch (num)
    {
    case 1:
      *mat.Dense += *mat.Dense;
      break;
    case 2:
      *mat.Triang += *mat.Triang;
      break;
    case 3:
      *mat.Sym += *mat.Sym;
      break;
    case 4:
      *mat.Sparse += *mat.Sparse;
      break;
    case 5:
      *mat.Banded += *mat.Banded;
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix op+= invalid type of matrix");
    }
    resetLU();
    return *this;
  }

  if (dimRow != m.size(0) || dimCol != m.size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix op+= inconsistent sizes.");

  if (numM == 0) // m is a BlockMatrix
  {
    ConstBlockIterator1 it1;
    ConstBlockIterator2 it2;
    unsigned int posRow = 0;
    unsigned int posCol = 0;
    // We scan all the blocks of m ...
    for (it1 = m.begin(); it1 != m.end(); ++it1)
    {
      for (it2 = it1.begin(); it2 != it1.end(); ++it2)
      {
        addBlock(posRow, posCol, **it2); // Each block of m is added into this.
        posCol += (*it2)->size(1);
      }
      posRow += (*it1)->size(0);
      posCol = 0;
    }
  }
  else // if m is a SimpleMatrix
  {
    switch (num)
    {
    case 1:
      switch (numM)
      {
      case 1:
        noalias(*(mat.Dense)) += *m.getDensePtr();
        break;
      case 2:
        noalias(*(mat.Dense)) += *m.getTriangPtr();
        break;
      case 3:
        noalias(*(mat.Dense)) += *m.getSymPtr();
        break;
      case 4:
        noalias(*(mat.Dense)) += *m.getSparsePtr();
        break;
      case 5:
        noalias(*(mat.Dense)) += *m.getBandedPtr();
        break;
      case 7:
        noalias(*(mat.Dense)) += *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op+= (const SimpleMatrix): invalid type of matrix");
        break;
      }
      break;
    case 2:
      switch (numM)
      {
      case 2:
        noalias(*(mat.Triang)) += *m.getTriangPtr();
        break;
      case 7:
        noalias(*(mat.Triang)) += *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op+= of a bad type of matrix into a triangular one.");
        break;
      }
      break;
    case 3:
      if (numM == 3)
        noalias(*(mat.Sym)) += *m.getSymPtr();
      else if (numM == 7)
        noalias(*(mat.Sym)) += *m.getIdentityPtr();
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::op+= bad assignment of matrix (symetric one = dense or ...)");
      break;
    case 4:
      switch (numM)
      {
      case 2:
        noalias(*(mat.Sparse)) += *m.getTriangPtr();
        break;
      case 3:
        noalias(*(mat.Sparse)) += *m.getSymPtr();
        break;
      case 4:
        noalias(*(mat.Sparse)) += *m.getSparsePtr();
        break;
      case 5:
        noalias(*(mat.Sparse)) += *m.getBandedPtr();
        break;
      case 7:
        noalias(*(mat.Sparse)) += *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op+=: invalid type of matrix");
        break;
      }
      break;
    case 5:
      switch (numM)
      {
      case 5:
        noalias(*(mat.Banded)) += *m.getBandedPtr();
        break;
      case 7:
        noalias(*(mat.Banded)) += *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op+= : invalid type of matrix");
        break;
      }
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op+= : invalid type of matrix");
      break;
    }
    resetLU();
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator -= (const SiconosMatrix& m)
{

  unsigned int numM = m.getNum();
  if (numM == 6) // m = 0
    return *this;

  if (&m == this) // auto-assignment
  {
    switch (num)
    {
    case 1:
      *mat.Dense -= *mat.Dense;
      break;
    case 2:
      *mat.Triang -= *mat.Triang;
      break;
    case 3:
      *mat.Sym -= *mat.Sym;
      break;
    case 4:
      *mat.Sparse -= *mat.Sparse;
      break;
    case 5:
      *mat.Banded -= *mat.Banded;
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix op-= invalid type of matrix");
    }
    resetLU();
    return *this;
  }
  if (dimRow != m.size(0) || dimCol != m.size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix op-= inconsistent sizes.");

  if (numM == 0) // m is a BlockMatrix
  {
    ConstBlockIterator1 it1;
    ConstBlockIterator2 it2;
    unsigned int posRow = 0;
    unsigned int posCol = 0;
    // We scan all the blocks of m ...
    for (it1 = m.begin(); it1 != m.end(); ++it1)
    {
      for (it2 = it1.begin(); it2 != it1.end(); ++it2)
      {
        subBlock(posRow, posCol, **it2); // Each block of m is added into this.
        posCol += (*it2)->size(1);
      }
      posRow += (*it1)->size(0);
      posCol = 0;
    }
  }
  else // if m is a SimpleMatrix
  {
    switch (num)
    {
    case 1:
      switch (numM)
      {
      case 1:
        noalias(*(mat.Dense)) -= *m.getDensePtr();
        break;
      case 2:
        noalias(*(mat.Dense)) -= *m.getTriangPtr();
        break;
      case 3:
        noalias(*(mat.Dense)) -= *m.getSymPtr();
        break;
      case 4:
        noalias(*(mat.Dense)) -= *m.getSparsePtr();
        break;
      case 5:
        noalias(*(mat.Dense)) -= *m.getBandedPtr();
        break;
      case 7:
        noalias(*(mat.Dense)) -= *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op-= (const SimpleMatrix): invalid type of matrix");
        break;
      }
      break;
    case 2:
      switch (numM)
      {
      case 2:
        noalias(*(mat.Triang)) -= *m.getTriangPtr();
        break;
      case 7:
        noalias(*(mat.Triang)) -= *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op-= of a bad type of matrix into a triangular one.");
        break;
      }
      break;
    case 3:
      if (numM == 3)
        noalias(*(mat.Sym)) -= *m.getSymPtr();
      else if (numM == 7)
        noalias(*(mat.Sym)) -= *m.getIdentityPtr();
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::op-= bad assignment of matrix (symetric one = dense or ...)");
      break;
    case 4:
      switch (numM)
      {
      case 2:
        noalias(*(mat.Sparse)) -= *m.getTriangPtr();
        break;
      case 3:
        noalias(*(mat.Sparse)) -= *m.getSymPtr();
        break;
      case 4:
        noalias(*(mat.Sparse)) -= *m.getSparsePtr();
        break;
      case 5:
        noalias(*(mat.Sparse)) -= *m.getBandedPtr();
        break;
      case 7:
        noalias(*(mat.Sparse)) -= *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op-=: invalid type of matrix");
        break;
      }
      break;
    case 5:
      switch (numM)
      {
      case 5:
        noalias(*(mat.Banded)) -= *m.getBandedPtr();
        break;
      case 7:
        noalias(*(mat.Banded)) -= *m.getIdentityPtr();
        break;
      default:
        SiconosMatrixException::selfThrow("SimpleMatrix::op-= : invalid type of matrix");
        break;
      }
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op-= : invalid type of matrix");
      break;
    }
    resetLU();
  }
  return *this;

}

void SimpleMatrix::trans()
{
  switch (num)
  {
  case 1:
    *mat.Dense = ublas::trans(*mat.Dense);
    break;
  case 2:
    SiconosMatrixException::selfThrow("SimpleMatrix::trans() failed, the matrix is triangular matrix and can not be transposed in place.");
    break;
  case 3:
    break;
  case 4:
    *mat.Sparse = ublas::trans(*mat.Sparse);
  case 5:
    *mat.Banded = ublas::trans(*mat.Banded);
    break;
  case 6:
    break;
  case 7:
    break;
  }
  unsigned int tmp = dimRow;
  dimRow = dimCol;
  dimCol = tmp;
  resetLU();
}

void SimpleMatrix::trans(const SiconosMatrix &m)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, not yet implemented for m being a BlockMatrix.");


  if (&m == this)
    trans();//SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, m = this, use this->trans().");
  else
  {
    unsigned int numM = m.getNum();
    switch (numM)
    {
    case 1:
      if (num != 1)
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a dense matrix into another type.");
      noalias(*mat.Dense) = ublas::trans(*m.getDensePtr());
      break;
    case 2:
      if (num != 1)
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a triangular matrix into a non-dense one.");
      noalias(*mat.Dense) = ublas::trans(*m.getTriangPtr());
      break;
    case 3:
      *this = m;
      break;
    case 4:
      if (num == 1)
        noalias(*mat.Dense) = ublas::trans(*m.getSparsePtr());
      else if (num == 4)
        noalias(*mat.Sparse) = ublas::trans(*m.getSparsePtr());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a sparse matrix into a forbidden type (not dense nor sparse).");
      break;
    case 5:
      if (num == 1)
        noalias(*mat.Dense) = ublas::trans(*m.getBandedPtr());
      else if (num == 5)
        noalias(*mat.Banded) = ublas::trans(*m.getBandedPtr());
      else
        SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a banded matrix into a forbidden type (not dense nor banded).");
      break;
    case 6:
      *this = m;
      break;
    case 7:
      *this = m;
    }
    unsigned int tmp = dimRow;
    dimRow = dimCol;
    dimCol = tmp;
    resetLU();
  }
}

void SimpleMatrix::PLUFactorizationInPlace()
{
  if (num != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUFactorizationInPlace: only implemented for dense matrices.");

  if (isPLUFactorized)
  {
    std::cout << "SimpleMatrix::PLUFactorizationInPlace warning: this matrix is already PLUFactorized. " << std::endl;
    return;
  }
  if (!ipiv)
    ipiv = new std::vector<int>(dimRow);
  else
    ipiv->resize(dimRow);
  int info = lapack::getrf(*mat.Dense, *ipiv);
  if (info != 0)
  {
    isPLUFactorized = false;
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace failed: the matrix is singular.");
  }
  //std::cout<<"SimpleMatrix::PLUFactorizationInPlace warning: the matrix is singular." << std::endl;
  else isPLUFactorized = true;
}

void SimpleMatrix::PLUInverseInPlace()
{
  if (!isPLUFactorized)
    PLUFactorizationInPlace();

#ifdef HAVE_ATLAS
  int info = lapack::getri(*mat.Dense, *ipiv);   // solve from factorization

  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUInverseInPlace failed, the matrix is singular.");

  isPLUInversed = true;
#else
  SiconosMatrixException::selfThrow("SimpleMatrix::PLUInverseInPlace not implemented with lapack.");
#endif
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(M) failed. Not yet implemented for M being a BlockMatrix.");
  int info;
  if (!isPLUFactorized) // call gesv => LU-factorize+solve
  {
    // solve system:
    if (!ipiv)
      ipiv = new std::vector<int>(dimRow);
    else
      ipiv->resize(dimRow);
    info = lapack::gesv(*mat.Dense, *ipiv, *(B.getDensePtr()));
    isPLUFactorized = true;
    // B now contains solution:
  }
  else // call getrs: only solve using previous lu-factorization
    info = lapack::getrs(*mat.Dense, *ipiv, *(B.getDensePtr()));

  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(V) failed. Not yet implemented for V being a BlockVector.");
  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.getDensePtr()); // Conversion of vector to matrix. Temporary solution.
  int info;
  if (!isPLUFactorized) // call gesv => LU-factorize+solve
  {
    // solve system:
    if (!ipiv)
      ipiv = new std::vector<int>(dimRow);
    else
      ipiv->resize(dimRow);
    info = lapack::gesv(*mat.Dense, *ipiv, tmpB);
    isPLUFactorized = true;
    // B now contains solution:
  }
  else // call getrs: only solve using previous lu-factorization
    info = lapack::getrs(*mat.Dense, *ipiv, tmpB);

  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
  *(B.getDensePtr()) = ublas::column(tmpB, 0);
}

void SimpleMatrix::resetLU()
{
  if (ipiv) ipiv->clear();
  isPLUFactorized = false;
  isPLUInversed = false;
}

const SimpleMatrix operator * (const SiconosMatrix & A, double a)
{
  // To compute B = a * A

  unsigned int numA = A.getNum();

  if (numA == 6) // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if (numA == 7)
  {
    return (DenseMat)(a**A.getIdentityPtr());
  }
  else if (numA == 0) // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp *= a;
    return tmp;
  }
  else if (numA == 1) // dense)
    return (DenseMat)(a** A.getDensePtr());
  else if (numA == 2)
    return (TriangMat)(a ** A.getTriangPtr());
  else if (numA == 3)
    return (SymMat)(a ** A.getSymPtr());
  else if (numA == 4)
    return (SparseMat)(a ** A.getSparsePtr());
  else //if(numA==5)
    return (BandedMat)(a ** A.getBandedPtr());
}

const SimpleMatrix operator * (double a, const SiconosMatrix & A)
{
  // To compute B = a * A

  unsigned int numA = A.getNum();

  if (numA == 6) // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if (numA == 7)
  {
    return (DenseMat)(a**A.getIdentityPtr());
  }
  else if (numA == 0) // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp *= a;
    return tmp;
  }
  else if (numA == 1) // dense)
    return (DenseMat)(a** A.getDensePtr());
  else if (numA == 2)
    return (TriangMat)(a ** A.getTriangPtr());
  else if (numA == 3)
    return (SymMat)(a ** A.getSymPtr());
  else if (numA == 4)
    return (SparseMat)(a ** A.getSparsePtr());
  else //if(numA==5)
    return (BandedMat)(a ** A.getBandedPtr());
}

const SimpleMatrix operator / (const SiconosMatrix & A, double a)
{
  // To compute B = A/a

  if (a == 0.0)
    SiconosMatrixException::selfThrow(" Matrix, operator / , division by zero.");

  unsigned int numA = A.getNum();

  if (numA == 6) // if A = 0
  {
    //DenseMat p(zero_matrix(A.size(0),A.size(1)));
    //return p;
    return A;
  }
  else if (numA == 7)
  {
    return (DenseMat)(*A.getIdentityPtr() / a);
  }
  else if (numA == 0) // A block
  {
    SimpleMatrix tmp(A); // ... copy ...
    tmp /= a;
    return tmp;
  }
  else if (numA == 1) // dense)
    return (DenseMat)(*A.getDensePtr() / a);
  else if (numA == 2)
    return (TriangMat)(*A.getTriangPtr() / a);
  else if (numA == 3)
    return (SymMat)(*A.getSymPtr() / a);
  else if (numA == 4)
    return (SparseMat)(*A.getSparsePtr() / a);
  else //if(numA==5)
    return (BandedMat)(*A.getBandedPtr() / a);
}

const SimpleMatrix operator + (const  SiconosMatrix& A, const  SiconosMatrix& B)
{
  // To compute C = A + B

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix operator +: inconsistent sizes");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();

  // == A or B equal to null ==
  if (numA == 6) // A = 0
  {
    if (numB == 6) // B = 0
      return SimpleMatrix(A.size(0), A.size(1));
    else
      return SimpleMatrix(B);
  }

  if (numB == 6)
    return SimpleMatrix(A);

  // == A and B different from 0 ==

  if (numA == numB && numA != 0) // all matrices are of the same type and NOT block
  {
    if (numA == 1)
      return (DenseMat)(*A.getDensePtr() + *B.getDensePtr());
    else if (numA == 2)
      return (TriangMat)(*A.getTriangPtr() + *B.getTriangPtr());
    else if (numA == 3)
      return (SymMat)(*A.getSymPtr() + *B.getSymPtr());
    else if (numA == 4)
    {
      SparseMat tmp(*A.getSparsePtr());
      tmp += *B.getSparsePtr();
      return tmp;
      // return (SparseMat)(*A.getSparsePtr() + *B.getSparsePtr());
    }
    else //if(numA==5)
    {
      BandedMat tmp(*A.getBandedPtr());
      tmp += *B.getBandedPtr();
      return tmp;
    }
  }
  else if (numA != 0 && numB != 0 && numA != numB) // A and B of different types and none is block
  {
    if (numA == 1)
    {
      if (numB == 2)
        return (DenseMat)(*A.getDensePtr() + *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getDensePtr() + *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getDensePtr() + *B.getSparsePtr());
      else if (numB == 5)
        return (DenseMat)(*A.getDensePtr() + *B.getBandedPtr());
      else // if(numB ==7)
        return (DenseMat)(*A.getDensePtr() + *B.getIdentityPtr());
      SiconosMatrixException::selfThrow("Matrix operator +: invalid type of matrix");
    }
    else if (numA == 2)
    {
      if (numB == 1)
        return (DenseMat)(*A.getTriangPtr() + *B.getDensePtr());
      else if (numB == 3)
        return (DenseMat)(*A.getTriangPtr() + *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getTriangPtr() + *B.getSparsePtr());
      else if (numB == 5)
        return (DenseMat)(*A.getTriangPtr() + *B.getBandedPtr());
      else // if(numB ==7:
        return (DenseMat)(*A.getTriangPtr() + *B.getIdentityPtr());
    }
    else if (numA == 3)
    {
      if (numB == 1)
        return (DenseMat)(*A.getSymPtr() + *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getSymPtr() + *B.getTriangPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getSymPtr() + *B.getSparsePtr());
      else if (numB == 5)
        return (DenseMat)(*A.getSymPtr() + *B.getBandedPtr());
      else // if(numB ==7)
        return (DenseMat)(*A.getSymPtr() + *B.getIdentityPtr());
    }
    else if (numA == 4)
    {
      if (numB == 1)
        return (DenseMat)(*A.getSparsePtr() + *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getSparsePtr() + *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getSparsePtr() + *B.getSymPtr());
      else if (numB == 5)
        return (DenseMat)(*A.getSparsePtr() + *B.getBandedPtr());
      else // if(numB ==7)
        return (DenseMat)(*A.getSparsePtr() + *B.getIdentityPtr());
    }

    else if (numA == 5)
    {
      if (numB == 1)
        return (DenseMat)(*A.getBandedPtr() + *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getBandedPtr() + *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getBandedPtr() + *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getBandedPtr() + *B.getSparsePtr());
      else //if(numB ==7)
        return (DenseMat)(*A.getBandedPtr() + *B.getIdentityPtr());
    }

    else //if(numA==7)
    {
      if (numB == 1)
        return (DenseMat)(*A.getIdentityPtr() + *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getIdentityPtr() + *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getIdentityPtr() + *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getIdentityPtr() + *B.getSparsePtr());
      else //if(numB ==5)
        return (DenseMat)(*A.getIdentityPtr() + *B.getBandedPtr());
    }
  }
  else if (numB != 0) // B Simple, whatever is A
  {
    SimpleMatrix tmp(B);
    tmp += A;
    return tmp;
  }
  else // B Block, A simple or block
  {
    SimpleMatrix tmp(A);
    tmp += B;
    return tmp;
  }
}

void add(const SiconosMatrix & A, const SiconosMatrix& B, SiconosMatrix& C)
{
  // To compute C = A + B in an "optimized" way (in comparison with operator +)

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");
  if ((A.size(0) != C.size(0)) || (A.size(1) != C.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();

  // === if C is zero or identity => read-only ===
  if (numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");

  // === common memory between A, B, C ===
  if (&A == &C) // A and C have common memory
  {
    C += B;
  }
  else if (&B == &C)  // B and C have common memory
  {
    C += A;
  }
  else // No common memory between C and A or B.
  {
    if (numA == 6) // A = 0
      C = B ;
    else if (numB == 6) // B = 0
      C = A;
    else // A and B different from 0
    {
      if (numC == 0) // if C is Block
      {
        if (numA != 0) // A simple, whatever is B
        {
          C = A;
          C += B;
        }
        else  // A Block
        {
          C = B;
          C += A;
        }
      }
      else // if C is a SimpleMatrix
      {
        if (numA == numB && numA != 0) // A and B are of the same type and NOT block
        {
          if (numC == numA)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) = *A.getDensePtr() + *B.getDensePtr();
            else if (numA == 2)
              noalias(*C.getTriangPtr()) = *A.getTriangPtr() + *B.getTriangPtr();
            else if (numA == 3)
              noalias(*C.getSymPtr()) = *A.getSymPtr() + *B.getSymPtr();
            else if (numA == 4)
              noalias(*C.getSparsePtr()) = *A.getSparsePtr() + *B.getSparsePtr();
            else //if(numA==5)
              noalias(*C.getBandedPtr()) = *A.getBandedPtr() + *B.getBandedPtr();
          }
          else // C and A of different types.
          {
            if (numC != 1)
              SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
            // Only dense matrices are allowed for output.

            if (numA == 1)
              noalias(*C.getDensePtr()) = *A.getDensePtr() + *B.getDensePtr();
            else if (numA == 2)
              noalias(*C.getDensePtr()) = *A.getTriangPtr() + *B.getTriangPtr();
            else if (numA == 3)
              noalias(*C.getDensePtr()) = *A.getSymPtr() + *B.getSymPtr();
            else if (numA == 4)
              noalias(*C.getDensePtr()) = *A.getSparsePtr() + *B.getSparsePtr();
            else //if(numA==5)
              noalias(*C.getDensePtr()) = *A.getBandedPtr() + *B.getBandedPtr();
          }
          C.resetLU();
        }
        else if (numA != 0 && numB != 0 && numA != numB) // A and B of different types and none is block
        {
          if (numC != 1)
            SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if (numA == 1)
            switch (numB)
            {
            case 2:
              noalias(*C.getDensePtr()) = *A.getDensePtr() + *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getDensePtr() + *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getDensePtr() + *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getDensePtr() + *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getDensePtr() + *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 2)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() + *B.getDensePtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() + *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() + *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() + *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() + *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 3)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getSymPtr() + *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getSymPtr() + *B.getTriangPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getSymPtr() + *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getSymPtr() + *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getSymPtr() + *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 4)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() + *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() + *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() + *B.getSymPtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() + *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() + *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 5)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() + *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() + *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() + *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() + *B.getSparsePtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() + *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 7)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() + *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() + *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() + *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() + *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() + *B.getBandedPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          C.resetLU();
        }
        else // A and/or B is Block
        {
          if (numA != 0) // A Simple, whatever is B
          {
            C = A;
            C += B;
          }
          else // A Block
          {
            C = B;
            C += A;
          }
        }
      }
    }
  }

}

const SimpleMatrix operator - (const  SiconosMatrix& A, const  SiconosMatrix& B)
{
  // To compute C = A - B

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix operator +: inconsistent sizes");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();

  // == B equal to null ==
  if (numB == 6)
    return SimpleMatrix(A);

  // == B different from 0 ==

  if (numA == numB && numA != 0) // all matrices are of the same type and NOT block
  {
    if (numA == 1)
      return (DenseMat)(*A.getDensePtr() - *B.getDensePtr());
    else if (numA == 2)
      return (TriangMat)(*A.getTriangPtr() - *B.getTriangPtr());
    else if (numA == 3)
      return (SymMat)(*A.getSymPtr() - *B.getSymPtr());
    else if (numA == 4)
    {
      SparseMat tmp(*A.getSparsePtr());
      tmp -= *B.getSparsePtr();
      return tmp;
      //return (SparseMat)(*A.getSparsePtr() - *B.getSparsePtr());
    }
    else //if(numA==5)
    {
      BandedMat tmp(*A.getBandedPtr());
      tmp -= *B.getBandedPtr();
      return tmp;
      //return (BandedMat)(*A.getBandedPtr() - *B.getBandedPtr());
    }
  }
  else if (numA != 0 && numB != 0 && numA != numB) // A and B of different types and none is block
  {
    if (numA == 1)
    {
      if (numB == 2)
        return (DenseMat)(*A.getDensePtr() - *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getDensePtr() - *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getDensePtr() - *B.getSparsePtr());
      else if (numB == 5)
        return (DenseMat)(*A.getDensePtr() - *B.getBandedPtr());
      else // if(numB ==7)
        return (DenseMat)(*A.getDensePtr() - *B.getIdentityPtr());
      SiconosMatrixException::selfThrow("Matrix operator +: invalid type of matrix");
    }
    else if (numA == 2)
    {
      if (numB == 1)
        return (DenseMat)(*A.getTriangPtr() - *B.getDensePtr());
      else if (numB == 3)
        return (DenseMat)(*A.getTriangPtr() - *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getTriangPtr() - *B.getSparsePtr());
      else if (numB == 5)
        return (DenseMat)(*A.getTriangPtr() - *B.getBandedPtr());
      else // if(numB ==7:
        return (DenseMat)(*A.getTriangPtr() - *B.getIdentityPtr());
    }
    else if (numA == 3)
    {
      if (numB == 1)
        return (DenseMat)(*A.getSymPtr() - *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getSymPtr() - *B.getTriangPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getSymPtr() - *B.getSparsePtr());
      else if (numB == 5)
        return (DenseMat)(*A.getSymPtr() - *B.getBandedPtr());
      else // if(numB ==7)
        return (DenseMat)(*A.getSymPtr() - *B.getIdentityPtr());
    }
    else if (numA == 4)
    {
      if (numB == 1)
        return (DenseMat)(*A.getSparsePtr() - *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getSparsePtr() - *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getSparsePtr() - *B.getSymPtr());
      else if (numB == 5)
        return (DenseMat)(*A.getSparsePtr() - *B.getBandedPtr());
      else // if(numB ==7)
        return (DenseMat)(*A.getSparsePtr() - *B.getIdentityPtr());
    }

    else if (numA == 5)
    {
      if (numB == 1)
        return (DenseMat)(*A.getBandedPtr() - *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getBandedPtr() - *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getBandedPtr() - *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getBandedPtr() - *B.getSparsePtr());
      else //if(numB ==7)
        return (DenseMat)(*A.getBandedPtr() - *B.getIdentityPtr());
    }

    else if (numA == 6)
    {
      if (numB == 1)
        return (DenseMat)(*A.getZeroPtr() - *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getZeroPtr() - *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getZeroPtr() - *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getZeroPtr() - *B.getSparsePtr());
      else //if(numB ==7)
        return (DenseMat)(*A.getZeroPtr() - *B.getIdentityPtr());
    }
    else //if(numA==7)
    {
      if (numB == 1)
        return (DenseMat)(*A.getIdentityPtr() - *B.getDensePtr());
      else if (numB == 2)
        return (DenseMat)(*A.getIdentityPtr() - *B.getTriangPtr());
      else if (numB == 3)
        return (DenseMat)(*A.getIdentityPtr() - *B.getSymPtr());
      else if (numB == 4)
        return (DenseMat)(*A.getIdentityPtr() - *B.getSparsePtr());
      else //if(numB ==5)
        return (DenseMat)(*A.getIdentityPtr() - *B.getBandedPtr());
    }
  }
  else // A and/or B are/is Block
  {
    SimpleMatrix tmp(A);
    tmp -= B;
    return tmp;
  }
}

void sub(const SiconosMatrix & A, const SiconosMatrix& B, SiconosMatrix& C)
{
  // To compute C = A - B in an "optimized" way (in comparison with operator +)

  if ((A.size(0) != B.size(0)) || (A.size(1) != B.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");
  if ((A.size(0) != C.size(0)) || (A.size(1) != C.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();

  // === if C is zero or identity => read-only ===
  if (numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");

  // === common memory between A, B, C ===
  if (&A == &C) // A and C have common memory
  {
    C -= B;
  }
  else if (&B == &C)  // B and C have common memory
  {
    if (numB == 0 || numA == 0) // if A or B(C) is Block
    {
      C *= -1.0;
      C += A;
    }
    else
    {
      if (numC == 0) // if C is Block
      {
        C = A;
        C -= B;
      }
      else // if C is a SimpleMatrix
      {
        if (numA == numB && numA != 0) // A and B are of the same type and NOT block
        {
          if (numA == 1)
            *C.getDensePtr() = *A.getDensePtr() - *B.getDensePtr();
          else if (numA == 2)
            *C.getTriangPtr() = *A.getTriangPtr() - *B.getTriangPtr();
          else if (numA == 3)
            *C.getSymPtr() = *A.getSymPtr() - *B.getSymPtr();
          else if (numA == 4)
            *C.getSparsePtr() = *A.getSparsePtr() - *B.getSparsePtr();
          else //if(numA==5)
            *C.getBandedPtr() = *A.getBandedPtr() - *B.getBandedPtr();
        }
        else if (numA != 0 && numB != 0 && numA != numB) // A and B of different types and none is block
        {
          if (numC != 1)  // => numB == 1
            SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if (numA == 1)
            *C.getDensePtr() = *A.getDensePtr() - *B.getDensePtr();
          else if (numA == 2)
            *C.getDensePtr() = *A.getTriangPtr() - *B.getDensePtr();
          else if (numA == 3)
            *C.getDensePtr() = *A.getSymPtr() - *B.getDensePtr();
          else if (numA == 4)
            *C.getDensePtr() = *A.getSparsePtr() - *B.getDensePtr();
          else if (numA == 5)
            *C.getDensePtr() = *A.getBandedPtr() - *B.getDensePtr();
          else if (numA == 6)
            *C.getDensePtr() = *A.getZeroPtr() - *B.getDensePtr();
          else //if(numA==7)
            *C.getDensePtr() = *A.getIdentityPtr() - *B.getDensePtr();
        }
        else // A and/or B is Block
        {
          C = A;
          C -= B;
        }
        C.resetLU();
      }
    }
  }
  else // No common memory between C and A or B.
  {
    if (numB == 6) // B = 0
      C = A;
    else // B different from 0
    {
      if (numC == 0) // if C is Block
      {
        C = A;
        C -= B;
      }
      else // if C is a SimpleMatrix
      {
        if (numA == numB && numA != 0) // A and B are of the same type and NOT block
        {
          if (numC == numA)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) = *A.getDensePtr() - *B.getDensePtr();
            else if (numA == 2)
              noalias(*C.getTriangPtr()) = *A.getTriangPtr() - *B.getTriangPtr();
            else if (numA == 3)
              noalias(*C.getSymPtr()) = *A.getSymPtr() - *B.getSymPtr();
            else if (numA == 4)
              noalias(*C.getSparsePtr()) = *A.getSparsePtr() - *B.getSparsePtr();
            else //if(numA==5)
              noalias(*C.getBandedPtr()) = *A.getBandedPtr() - *B.getBandedPtr();
          }
          else // C and A of different types.
          {
            if (numC != 1)
              SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
            // Only dense matrices are allowed for output.

            if (numA == 1)
              noalias(*C.getDensePtr()) = *A.getDensePtr() - *B.getDensePtr();
            else if (numA == 2)
              noalias(*C.getDensePtr()) = *A.getTriangPtr() - *B.getTriangPtr();
            else if (numA == 3)
              noalias(*C.getDensePtr()) = *A.getSymPtr() - *B.getSymPtr();
            else if (numA == 4)
              noalias(*C.getDensePtr()) = *A.getSparsePtr() - *B.getSparsePtr();
            else //if(numA==5)
              noalias(*C.getDensePtr()) = *A.getBandedPtr() - *B.getBandedPtr();
          }
          C.resetLU();
        }
        else if (numA != 0 && numB != 0 && numA != numB) // A and B of different types and none is block
        {
          if (numC != 1)
            SiconosMatrixException::selfThrow("Matrix addition ( add(A,B,C) ): wrong type for resulting matrix C.");
          // Only dense matrices are allowed for output.

          if (numA == 1)
            switch (numB)
            {
            case 2:
              noalias(*C.getDensePtr()) = *A.getDensePtr() - *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getDensePtr() - *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getDensePtr() - *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getDensePtr() - *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getDensePtr() - *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 2)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() - *B.getDensePtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() - *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() - *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() - *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getTriangPtr() - *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 3)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getSymPtr() - *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getSymPtr() - *B.getTriangPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getSymPtr() - *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getSymPtr() - *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getSymPtr() - *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 4)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() - *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() - *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() - *B.getSymPtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() - *B.getBandedPtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getSparsePtr() - *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 5)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() - *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() - *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() - *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() - *B.getSparsePtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getBandedPtr() - *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 6)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getZeroPtr() - *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getZeroPtr() - *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getZeroPtr() - *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getZeroPtr() - *B.getSparsePtr();
              break;
            case 7:
              noalias(*C.getDensePtr()) = *A.getZeroPtr() - *B.getIdentityPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else if (numA == 7)
            switch (numB)
            {
            case 1:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() - *B.getDensePtr();
              break;
            case 2:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() - *B.getTriangPtr();
              break;
            case 3:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() - *B.getSymPtr();
              break;
            case 4:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() - *B.getSparsePtr();
              break;
            case 5:
              noalias(*C.getDensePtr()) = *A.getIdentityPtr() - *B.getBandedPtr();
              break;
            default:
              SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
            }
          else
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          C.resetLU();
        }
        else // A and/or B is Block
        {
          C = A;
          C -= B;
        }
      }
    }
  }
}

//========================
// Matrices comparison
//========================

bool operator == (const SiconosMatrix &m, const SiconosMatrix &x)
{
  //  if( ! isComparableTo( &m, &x))
  //    return false;
  // Warning : two block matrices may be "equal" but have blocks of different sizes.
  double norm = (m - x).normInf();
  return (norm < tolerance);
}

//======================
// Product of matrices
//======================

const SimpleMatrix prod(const SiconosMatrix &A, const SiconosMatrix& B)
{
  // To compute C = A * B

  if ((A.size(1) != B.size(0)))
    SiconosMatrixException::selfThrow("Matrix function C=prod(A,B): inconsistent sizes");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();

  // == TODO: implement block product ==
  if (numA == 0 || numB == 0)
    SiconosMatrixException::selfThrow("Matrix product ( C=prod(A,B) ): not yet implemented for BlockMatrix objects.");

  if (numA == 7 || numB == 6) // A = identity or B = 0
    return SimpleMatrix(B);

  else if (numB == 7 || numA == 6) // B = identity or A = 0
    return SimpleMatrix(A);

  else // neither A or B is equal to identity or zero.
  {
    if (numB == 1)
    {
      if (numA == 1)
      {
        DenseMat p(A.size(0), B.size(1));
        atlas::gemm(*A.getDensePtr(), *B.getDensePtr(), p);
        //      return (DenseMat)(prod(*A.getDensePtr(),*B.getDensePtr()));
        return p;
      }
      else if (numA == 2)
        return (DenseMat)(prod(*A.getTriangPtr(), *B.getDensePtr()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.getSymPtr(), *B.getDensePtr()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.getSparsePtr(), *B.getDensePtr()));
      else// if(numA==5)
        return (DenseMat)(prod(*A.getBandedPtr(), *B.getDensePtr()));
    }
    else if (numB == 2)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.getDensePtr(), *B.getTriangPtr()));
      else if (numA == 2)
        return (TriangMat)(prod(*A.getTriangPtr(), *B.getTriangPtr()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.getSymPtr(), *B.getTriangPtr()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.getSparsePtr(), *B.getTriangPtr()));
      else //if(numA==5)
        return (DenseMat)(prod(*A.getBandedPtr(), *B.getTriangPtr()));
    }
    else if (numB == 3)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.getDensePtr(), *B.getSymPtr()));
      else if (numA == 2)
        return (DenseMat)(prod(*A.getTriangPtr(), *B.getSymPtr()));
      else if (numA == 3)
        return (SymMat)(prod(*A.getSymPtr(), *B.getSymPtr()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.getSparsePtr(), *B.getSymPtr()));
      else // if (numA == 5)
        return (DenseMat)(prod(*A.getBandedPtr(), *B.getSymPtr()));
    }
    else if (numB == 4)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.getDensePtr(), *B.getSparsePtr()));
      else if (numA == 2)
        return (DenseMat)(prod(*A.getTriangPtr(), *B.getSparsePtr()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.getSymPtr(), *B.getSparsePtr()));
      else if (numA == 4)
        return (SparseMat)(prod(*A.getSparsePtr(), *B.getSparsePtr()));
      else //if(numA==5){
        return (DenseMat)(prod(*A.getBandedPtr(), *B.getSparsePtr()));
    }
    else //if(numB==5)
    {
      if (numA == 1)
        return (DenseMat)(prod(*A.getDensePtr(), *B.getBandedPtr()));
      else if (numA == 2)
        return (DenseMat)(prod(*A.getTriangPtr(), *B.getBandedPtr()));
      else if (numA == 3)
        return (DenseMat)(prod(*A.getSymPtr(), *B.getBandedPtr()));
      else if (numA == 4)
        return (DenseMat)(prod(*A.getSparsePtr(), *B.getBandedPtr()));
      else //if(numA==5)
        return (DenseMat)(prod(*A.getBandedPtr(), *B.getBandedPtr()));
    }
  }
}

void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init)
{
  // To compute C = A * B

  if ((A.size(1) != B.size(0)))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): inconsistent sizes");

  if (A.size(0) != C.size(0) || B.size(1) != C.size(1))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): inconsistent sizes");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();

  // == TODO: implement block product ==
  if (numA == 0 || numB == 0)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): not yet implemented for BlockMatrix objects.");

  // === if C is zero or identity => read-only ===
  if (numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");


  if (numA == 7) // A = identity ...
  {
    if (init)
    {
      if (&C != &B) C = B; // if C and B are two different objects.
      // else nothing
    }
    else
      C += B;
  }

  else if (numB == 7) // B = identity
  {
    if (init)
    {
      if (&C != &A) C = A; // if C and A are two different objects.
      // else nothing
    }
    else
      C += A;
  }

  else if (numA == 6 || numB == 6) // if A or B = 0
  {
    if (init)
      C.zero();
    //else nothing
  }
  else if (numC == 0) // if C is Block - Temp. solution
  {
    SimpleMatrix tmp(C);
    prod(A, B, tmp, init);
    C = tmp;
  }
  else // neither A or B is equal to identity or zero.
  {
    if (init)
    {
      if (&C == &A) // if common memory between A and C
      {
        switch (numA)
        {
        case 1:
          if (numB == 1)
            //*C.getDensePtr() = prod(*A.getDensePtr(),*B.getDensePtr());
            atlas::gemm(*A.getDensePtr(), *B.getDensePtr(), *C.getDensePtr());
          else if (numB == 2)
            *C.getDensePtr()  = prod(*A.getDensePtr(), *B.getTriangPtr());
          else if (numB == 3)
            *C.getDensePtr()  = prod(*A.getDensePtr(), *B.getSymPtr());
          else if (numB == 4)
            *C.getDensePtr()  = prod(*A.getDensePtr(), *B.getSparsePtr());
          else //if(numB==5)
            *C.getDensePtr() = prod(*A.getDensePtr(), *B.getBandedPtr());
          break;
        case 2:
          if (numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getTriangPtr() = prod(*A.getTriangPtr(), *B.getTriangPtr());
          break;
        case 3:
          if (numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSymPtr() = prod(*A.getSymPtr(), *B.getSymPtr());
          break;
        case 4:
          if (numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSparsePtr() = prod(*A.getSparsePtr(), *B.getSparsePtr());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else if (&C == &B)
      {
        switch (numB)
        {
        case 1:
          if (numA == 1)
            *C.getDensePtr() = prod(*A.getDensePtr(), *B.getDensePtr());
          else if (numA == 2)
            *C.getDensePtr()  = prod(*A.getTriangPtr(), *B.getDensePtr());
          else if (numA == 3)
            *C.getDensePtr()  = prod(*A.getSymPtr(), *B.getDensePtr());
          else if (numA == 4)
            *C.getDensePtr()  = prod(*A.getSparsePtr(), *B.getDensePtr());
          else //if(numB==5)
            *C.getDensePtr() = prod(*A.getBandedPtr(), *B.getDensePtr());
          break;
        case 2:
          if (numA != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getTriangPtr() = prod(*A.getTriangPtr(), *B.getTriangPtr());
          break;
        case 3:
          if (numA != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSymPtr() = prod(*A.getSymPtr(), *B.getSymPtr());
          break;
        case 4:
          if (numA != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSparsePtr() = prod(*A.getSparsePtr(), *B.getSparsePtr());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else // if no alias between C and A or B.
      {
        switch (numC)
        {
        case 1:
          if (numB == 1)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) = prod(*A.getDensePtr(), *B.getDensePtr());
            else if (numA == 2)
              noalias(*C.getDensePtr()) = prod(*A.getTriangPtr(), *B.getDensePtr());
            else if (numA == 3)
              noalias(*C.getDensePtr())  = prod(*A.getSymPtr(), *B.getDensePtr());
            else if (numA == 4)
              noalias(*C.getDensePtr()) = prod(*A.getSparsePtr(), *B.getDensePtr());
            else// if(numA==5)
              noalias(*C.getDensePtr())  = prod(*A.getBandedPtr(), *B.getDensePtr());
          }
          else if (numB == 2)
          {
            if (numA == 1)
              noalias(*C.getDensePtr())  = prod(*A.getDensePtr(), *B.getTriangPtr());
            else if (numA == 2)
              noalias(*C.getDensePtr())  = prod(*A.getTriangPtr(), *B.getTriangPtr());
            else if (numA == 3)
              noalias(*C.getDensePtr())  = prod(*A.getSymPtr(), *B.getTriangPtr());
            else if (numA == 4)
              noalias(*C.getDensePtr())  = prod(*A.getSparsePtr(), *B.getTriangPtr());
            else //if(numA==5)
              noalias(*C.getDensePtr())  = prod(*A.getBandedPtr(), *B.getTriangPtr());
          }
          else if (numB == 3)
          {
            if (numA == 1)
              noalias(*C.getDensePtr())  = prod(*A.getDensePtr(), *B.getSymPtr());
            else if (numA == 2)
              noalias(*C.getDensePtr())  = prod(*A.getTriangPtr(), *B.getSymPtr());
            else if (numA == 3)
              noalias(*C.getDensePtr())  = prod(*A.getSymPtr(), *B.getSymPtr());
            else if (numA == 4)
              noalias(*C.getDensePtr())  = prod(*A.getSparsePtr(), *B.getSymPtr());
            else // if (numA == 5)
              noalias(*C.getDensePtr())  = prod(*A.getBandedPtr(), *B.getSymPtr());
          }
          else if (numB == 4)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) = prod(*A.getDensePtr(), *B.getSparsePtr());
            else if (numA == 2)
              noalias(*C.getDensePtr()) = prod(*A.getTriangPtr(), *B.getSparsePtr());
            else if (numA == 3)
              noalias(*C.getDensePtr()) = prod(*A.getSymPtr(), *B.getSparsePtr());
            else if (numA == 4)
              noalias(*C.getDensePtr()) = prod(*A.getSparsePtr(), *B.getSparsePtr());
            else //if(numA==5){
              noalias(*C.getDensePtr()) = prod(*A.getBandedPtr(), *B.getSparsePtr());
          }
          else //if(numB==5)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) = prod(*A.getDensePtr(), *B.getBandedPtr());
            else if (numA == 2)
              noalias(*C.getDensePtr()) = prod(*A.getTriangPtr(), *B.getBandedPtr());
            else if (numA == 3)
              noalias(*C.getDensePtr()) = prod(*A.getSymPtr(), *B.getBandedPtr());
            else if (numA == 4)
              noalias(*C.getDensePtr()) = prod(*A.getSparsePtr(), *B.getBandedPtr());
            else //if(numA==5)
              noalias(*C.getDensePtr()) = prod(*A.getBandedPtr(), *B.getBandedPtr());
          }
          break;
        case 2:
          if (numA != 2 || numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.getTriangPtr()) = prod(*A.getTriangPtr(), *B.getTriangPtr());
          break;
        case 3:
          if (numA != 3 || numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.getSymPtr()) = prod(*A.getSymPtr(), *B.getSymPtr());
          break;
        case 4:
          if (numA != 4 || numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.getSparsePtr()) = prod(*A.getSparsePtr(), *B.getSparsePtr());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
    }
    else // += case
    {
      if (&C == &A) // if common memory between A and C
      {
        switch (numA)
        {
        case 1:
          if (numB == 1)
            *C.getDensePtr() += prod(*A.getDensePtr(), *B.getDensePtr());
          else if (numB == 2)
            *C.getDensePtr()  += prod(*A.getDensePtr(), *B.getTriangPtr());
          else if (numB == 3)
            *C.getDensePtr()  += prod(*A.getDensePtr(), *B.getSymPtr());
          else if (numB == 4)
            *C.getDensePtr()  += prod(*A.getDensePtr(), *B.getSparsePtr());
          else //if(numB==5)
            *C.getDensePtr() += prod(*A.getDensePtr(), *B.getBandedPtr());
          break;
        case 2:
          if (numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getTriangPtr() += prod(*A.getTriangPtr(), *B.getTriangPtr());
          break;
        case 3:
          if (numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSymPtr() += prod(*A.getSymPtr(), *B.getSymPtr());
          break;
        case 4:
          if (numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSparsePtr() += prod(*A.getSparsePtr(), *B.getSparsePtr());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else if (&C == &B)
      {
        switch (numB)
        {
        case 1:
          if (numA == 1)
            *C.getDensePtr() += prod(*A.getDensePtr(), *B.getDensePtr());
          else if (numA == 2)
            *C.getDensePtr()  += prod(*A.getTriangPtr(), *B.getDensePtr());
          else if (numA == 3)
            *C.getDensePtr()  += prod(*A.getSymPtr(), *B.getDensePtr());
          else if (numA == 4)
            *C.getDensePtr()  += prod(*A.getSparsePtr(), *B.getDensePtr());
          else //if(numB==5)
            *C.getDensePtr() += prod(*A.getBandedPtr(), *B.getDensePtr());
          break;
        case 2:
          if (numA != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getTriangPtr() += prod(*A.getTriangPtr(), *B.getTriangPtr());
          break;
        case 3:
          if (numA != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSymPtr() += prod(*A.getSymPtr(), *B.getSymPtr());
          break;
        case 4:
          if (numA != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          *C.getSparsePtr() += prod(*A.getSparsePtr(), *B.getSparsePtr());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
      else // if no alias between C and A or B.
      {
        switch (numC)
        {
        case 1:
          if (numB == 1)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) += prod(*A.getDensePtr(), *B.getDensePtr());
            else if (numA == 2)
              noalias(*C.getDensePtr()) += prod(*A.getTriangPtr(), *B.getDensePtr());
            else if (numA == 3)
              noalias(*C.getDensePtr())  += prod(*A.getSymPtr(), *B.getDensePtr());
            else if (numA == 4)
              noalias(*C.getDensePtr()) += prod(*A.getSparsePtr(), *B.getDensePtr());
            else// if(numA==5)
              noalias(*C.getDensePtr())  += prod(*A.getBandedPtr(), *B.getDensePtr());
          }
          else if (numB == 2)
          {
            if (numA == 1)
              noalias(*C.getDensePtr())  += prod(*A.getDensePtr(), *B.getTriangPtr());
            else if (numA == 2)
              noalias(*C.getDensePtr())  += prod(*A.getTriangPtr(), *B.getTriangPtr());
            else if (numA == 3)
              noalias(*C.getDensePtr())  += prod(*A.getSymPtr(), *B.getTriangPtr());
            else if (numA == 4)
              noalias(*C.getDensePtr())  += prod(*A.getSparsePtr(), *B.getTriangPtr());
            else //if(numA==5)
              noalias(*C.getDensePtr())  += prod(*A.getBandedPtr(), *B.getTriangPtr());
          }
          else if (numB == 3)
          {
            if (numA == 1)
              noalias(*C.getDensePtr())  += prod(*A.getDensePtr(), *B.getSymPtr());
            else if (numA == 2)
              noalias(*C.getDensePtr())  += prod(*A.getTriangPtr(), *B.getSymPtr());
            else if (numA == 3)
              noalias(*C.getDensePtr())  += prod(*A.getSymPtr(), *B.getSymPtr());
            else if (numA == 4)
              noalias(*C.getDensePtr())  += prod(*A.getSparsePtr(), *B.getSymPtr());
            else // if (numA == 5)
              noalias(*C.getDensePtr())  += prod(*A.getBandedPtr(), *B.getSymPtr());
          }
          else if (numB == 4)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) += prod(*A.getDensePtr(), *B.getSparsePtr());
            else if (numA == 2)
              noalias(*C.getDensePtr()) += prod(*A.getTriangPtr(), *B.getSparsePtr());
            else if (numA == 3)
              noalias(*C.getDensePtr()) += prod(*A.getSymPtr(), *B.getSparsePtr());
            else if (numA == 4)
              noalias(*C.getDensePtr()) += prod(*A.getSparsePtr(), *B.getSparsePtr());
            else //if(numA==5){
              noalias(*C.getDensePtr()) += prod(*A.getBandedPtr(), *B.getSparsePtr());
          }
          else //if(numB==5)
          {
            if (numA == 1)
              noalias(*C.getDensePtr()) += prod(*A.getDensePtr(), *B.getBandedPtr());
            else if (numA == 2)
              noalias(*C.getDensePtr()) += prod(*A.getTriangPtr(), *B.getBandedPtr());
            else if (numA == 3)
              noalias(*C.getDensePtr()) += prod(*A.getSymPtr(), *B.getBandedPtr());
            else if (numA == 4)
              noalias(*C.getDensePtr()) += prod(*A.getSparsePtr(), *B.getBandedPtr());
            else //if(numA==5)
              noalias(*C.getDensePtr()) += prod(*A.getBandedPtr(), *B.getBandedPtr());
          }
          break;
        case 2:
          if (numA != 2 || numB != 2)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.getTriangPtr()) += prod(*A.getTriangPtr(), *B.getTriangPtr());
          break;
        case 3:
          if (numA != 3 || numB != 3)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.getSymPtr()) += prod(*A.getSymPtr(), *B.getSymPtr());
          break;
        case 4:
          if (numA != 4 || numB != 4)
            SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
          noalias(*C.getSparsePtr()) += prod(*A.getSparsePtr(), *B.getSparsePtr());
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): wrong type for C (according to A and B types).");
        }
      }
    }
    C.resetLU();
  }
}


void axpy_prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init)
{
  // To compute C = A * B (init = true) or C += A * B (init = false) using ublas axpy_prod.
  if ((A.size(1) != B.size(0)))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): inconsistent sizes");

  if (A.size(0) != C.size(0) || B.size(1) != C.size(1))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C): inconsistent sizes");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  // == TODO: implement block product ==
  if (numA == 0 || numB == 0)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): not yet implemented for BlockMatrix objects.");

  // === if C is zero or identity => read-only ===
  if (numC == 6 || numC == 7)
    SiconosMatrixException::selfThrow("Matrix product ( prod(A,B,C) ): wrong type for resulting matrix C (read-only: zero or identity).");


  if (numA == 7) // A = identity ...
  {
    if (!init) C += B;
    else
    {
      if (&C != &B)
        C = B; // if C and B are two different objects.
      // else nothing
    }
  }

  else if (numB == 7) // B = identity
  {
    if (!init) C += A;
    else
    {
      if (&C != &A) C = A; // if C and A are two different objects.
      // else nothing
    }
  }

  else if (numA == 6 || numB == 6) // if A or B = 0
  {
    if (init) C.zero(); // else nothing
  }
  else if (numC == 0) // if C is Block - Temp. solution
  {
    SimpleMatrix tmp(C);
    axpy_prod(A, B, tmp, init);
    C = tmp;
  }
  else // neither A or B is equal to identity or zero.
  {
    if (&C == &A) // if common memory between A and C
    {
      switch (numA)
      {
      case 1:
        if (numB == 1)
          ublas::axpy_prod(*A.getDensePtr(), *B.getDensePtr(), *A.getDensePtr(), init);
        else if (numB == 2)
          ublas::axpy_prod(*A.getDensePtr(), *B.getTriangPtr(), *A.getDensePtr(), init);
        else if (numB == 3)
          ublas::axpy_prod(*A.getDensePtr(), *B.getSymPtr(), *A.getDensePtr(), init);
        else if (numB == 4)
          ublas::axpy_prod(*A.getDensePtr(), *B.getSparsePtr(), *A.getDensePtr(), init);
        else //if(numB==5)
          ublas::axpy_prod(*A.getDensePtr(), *B.getBandedPtr(), *A.getDensePtr(), init);
        break;
      case 2:
        //        if(numB != 2)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //        ublas::axpy_prod(*A.getTriangPtr(), *B.getTriangPtr(), *A.getTriangPtr(), init);
        break;
      case 3:
        //if(numB != 3)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //ublas::axpy_prod(*A.getSymPtr(), *B.getSymPtr(), *A.getSymPtr(), init);
        break;
      case 4:
        //        if(numB != 4)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //ublas::axpy_prod(*A.getSparsePtr(), *B.getSparsePtr(), *A.getSparsePtr(),init);
        break;
      default:
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      }
    }
    else if (&C == &B)
    {
      switch (numB)
      {
      case 1:
        if (numA == 1)
          ublas::axpy_prod(*A.getDensePtr(), *B.getDensePtr(), *B.getDensePtr(), init);
        else if (numA == 2)
          ublas::axpy_prod(*A.getTriangPtr(), *B.getDensePtr(), *B.getDensePtr(), init);
        else if (numA == 3)
          ublas::axpy_prod(*A.getSymPtr(), *B.getDensePtr(), *B.getDensePtr(), init);
        else if (numA == 4)
          ublas::axpy_prod(*A.getSparsePtr(), *B.getDensePtr(), *B.getDensePtr(), init);
        else //if(numB==5)
          ublas::axpy_prod(*A.getBandedPtr(), *B.getDensePtr(), *B.getDensePtr(), init);
        break;
      case 2:
        //if(numA != 2)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //        ublas::axpy_prod(*A.getTriangPtr(), *B.getTriangPtr(),*B.getTriangPtr(), init);
        break;
      case 3:
        //        if(numA != 3)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //ublas::axpy_prod(*A.getSymPtr(), *B.getSymPtr(), *B.getSymPtr(), init);
        break;
      case 4:
        //        if(numA != 4)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //ublas::axpy_prod(*A.getSparsePtr(), *B.getSparsePtr(), *B.getSparsePtr(), init);
        break;
      default:
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      }
    }
    else // if no alias between C and A or B.
    {
      switch (numC)
      {
      case 1:
        if (numB == 1)
        {
          if (numA == 1)
            ublas::axpy_prod(*A.getDensePtr(), *B.getDensePtr(), *C.getDensePtr(), init);
          else if (numA == 2)
            ublas::axpy_prod(*A.getTriangPtr(), *B.getDensePtr(), *C.getDensePtr(), init);
          else if (numA == 3)
            ublas::axpy_prod(*A.getSymPtr(), *B.getDensePtr(), *C.getDensePtr(), init);
          else if (numA == 4)
            ublas::axpy_prod(*A.getSparsePtr(), *B.getDensePtr(), *C.getDensePtr(), init);
          else// if(numA==5)
            ublas::axpy_prod(*A.getBandedPtr(), *B.getDensePtr(), *C.getDensePtr(), init);
        }
        else if (numB == 2)
        {
          if (numA == 1)
            ublas::axpy_prod(*A.getDensePtr(), *B.getTriangPtr(), *C.getDensePtr(), init);
          else if (numA == 2)
            ublas::axpy_prod(*A.getTriangPtr(), *B.getTriangPtr(), *C.getDensePtr(), init);
          else if (numA == 3)
            ublas::axpy_prod(*A.getSymPtr(), *B.getTriangPtr(), *C.getDensePtr(), init);
          else if (numA == 4)
            ublas::axpy_prod(*A.getSparsePtr(), *B.getTriangPtr(), *C.getDensePtr(), init);
          else //if(numA==5)
            ublas::axpy_prod(*A.getBandedPtr(), *B.getTriangPtr(), *C.getDensePtr(), init);
        }
        else if (numB == 3)
        {
          if (numA == 1)
            ublas::axpy_prod(*A.getDensePtr(), *B.getSymPtr(), *C.getDensePtr(), init);
          else if (numA == 2)
            ublas::axpy_prod(*A.getTriangPtr(), *B.getSymPtr(), *C.getDensePtr(), init);
          else if (numA == 3)
            ublas::axpy_prod(*A.getSymPtr(), *B.getSymPtr(), *C.getDensePtr(), init);
          else if (numA == 4)
            ublas::axpy_prod(*A.getSparsePtr(), *B.getSymPtr(), *C.getDensePtr(), init);
          else // if (numA == 5)
            ublas::axpy_prod(*A.getBandedPtr(), *B.getSymPtr(), *C.getDensePtr(), init);
        }
        else if (numB == 4)
        {
          if (numA == 1)
            ublas::axpy_prod(*A.getDensePtr(), *B.getSparsePtr(), *C.getDensePtr(), init);
          else if (numA == 2)
            ublas::axpy_prod(*A.getTriangPtr(), *B.getSparsePtr(), *C.getDensePtr(), init);
          else if (numA == 3)
            ublas::axpy_prod(*A.getSymPtr(), *B.getSparsePtr(), *C.getDensePtr(), init);
          else if (numA == 4)
            ublas::axpy_prod(*A.getSparsePtr(), *B.getSparsePtr(), *C.getDensePtr(), init);
          else //if(numA==5){
            ublas::axpy_prod(*A.getBandedPtr(), *B.getSparsePtr(), *C.getDensePtr(), init);
        }
        else //if(numB==5)
        {
          if (numA == 1)
            ublas::axpy_prod(*A.getDensePtr(), *B.getBandedPtr(), *C.getDensePtr(), init);
          else if (numA == 2)
            ublas::axpy_prod(*A.getTriangPtr(), *B.getBandedPtr(), *C.getDensePtr(), init);
          else if (numA == 3)
            ublas::axpy_prod(*A.getSymPtr(), *B.getBandedPtr(), *C.getDensePtr(), init);
          else if (numA == 4)
            ublas::axpy_prod(*A.getSparsePtr(), *B.getBandedPtr(), *C.getDensePtr(), init);
          else //if(numA==5)
            ublas::axpy_prod(*A.getBandedPtr(), *B.getBandedPtr(), *C.getDensePtr(), init);
        }
        break;
      case 2:
        // if(numA!= 2 || numB != 2)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //ublas::axpy_prod(*A.getTriangPtr(), *B.getTriangPtr(),*C.getTriangPtr(), init);
        break;
      case 3:
        //        if(numA!= 3 || numB != 3)
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        //ublas::axpy_prod(*A.getSymPtr(), *B.getSymPtr(),*C.getSymPtr(),init);
        break;
      case 4:
        if (numA != 4 || numB != 4)
          SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
        ublas::sparse_prod(*A.getSparsePtr(), *B.getSparsePtr(), *C.getSparsePtr(), init);
        break;
      default:
        SiconosMatrixException::selfThrow("Matrix function axpy_prod(A,B,C): wrong type for C (according to A and B types).");
      }
    }
    C.resetLU();
  }
}


const SimpleMatrix pow(const SimpleMatrix& m, unsigned int power)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, pow function: not yet implemented for BlockMatrix.");
  if (!m.isSquare())
    SiconosMatrixException::selfThrow("pow(SimpleMatrix), matrix is not square.");

  if (power < 0)
    SiconosMatrixException::selfThrow("pow(SimpleMatrix,n) with negative value is not supported");

  if (power > 0)
  {
    unsigned int num = m.getNum();
    if (num == 1)
    {
      DenseMat p = *m.getDensePtr();
      for (unsigned int i = 1; i < power; i++)
        p = prod(p, *m.getDensePtr());
      return p;
    }
    else if (num == 2)
    {
      TriangMat t = *m.getTriangPtr();
      for (unsigned int i = 1; i < power; i++)
        t = prod(t, *m.getTriangPtr());
      return t;
    }
    else if (num == 3)
    {
      SymMat s = *m.getSymPtr();
      for (unsigned int i = 1; i < power; i++)
        s = prod(s, *m.getSymPtr());
      return s;
    }
    else if (num == 4)
    {
      SparseMat sp = *m.getSparsePtr();
      for (unsigned int i = 1; i < power; i++)
        sp = prod(sp, *m.getSparsePtr());
      return sp;
    }
    else if (num == 5)
    {
      DenseMat b = *m.getBandedPtr();
      for (unsigned int i = 1; i < power; i++)
        b = prod(b, *m.getBandedPtr());
      return b;
    }
    else if (num == 6)
    {
      ZeroMat z(m.size(0), m.size(1));
      return z;
    }
    else // if (num==7)
    {
      IdentityMat I(m.size(0), m.size(1));;
      return I;
    }
  }
  else// if(power == 0)
  {
    IdentityMat I = ublas::identity_matrix<double>(m.size(0), m.size(1));
    return I;
  }
}

// ========== Products matrix - vector //

// Computation of y = A*x
//
// Two specific functions are used to handle all the cases where x or y are blocks.
// All of their blocks can also be blocks. Then we use:
// - private_prod to "slice" A when y is block, ie according to its rows.
// - private_addprod to "slice" A when x is block, and to sum over the columns of blocks to compute y = sum subA x[i].

// The following function is private and used inside prod(...) public functions.
// It is required to deal with block vectors of blocks ( of blocks ...).
// It computes res = subA*x +res, subA being a submatrix of A (rows from startRow to startRow+sizeY and columns between startCol and startCol+sizeX).
// If x is a block vector, it call the present function for all blocks.
void private_addprod(SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SPC::SiconosVector x, SP::SiconosVector y)
{
  if (A->isBlock())
    SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for block matrix.");

  if (!x->isBlock()) // if input vector is not block
  {
    // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
    // Then computation of y = subA*x + y.
    unsigned int numA = A->getNum();
    unsigned int numY = y->getNum();
    unsigned int numX = x->getNum();
    unsigned int sizeX = x->size();
    unsigned int sizeY = y->size();

    if (numX != numY)
      SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x and y of different types.");

    if (numY == 1 && numX == 1)
    {
      if (numA == 1)
        *y->getDensePtr() += prod(ublas::subrange(*A->getDensePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 2)
        *y->getDensePtr() += prod(ublas::subrange(*A->getTriangPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 3)
        *y->getDensePtr() += prod(ublas::subrange(*A->getSymPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 4)
        *y->getDensePtr() += prod(ublas::subrange(*A->getSparsePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else //if(numA==5)
        *y->getDensePtr() += prod(ublas::subrange(*A->getBandedPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
    }
    else // x and y sparse
    {
      if (numA == 4)
        *y->getSparsePtr() += prod(ublas::subrange(*A->getSparsePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getSparsePtr());
      else
        SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x, y  sparse and A not sparse.");
    }

  }
  else // if block
  {
    VectorOfVectors::const_iterator it;
    unsigned int startColBis = startCol;
    for (it = x->begin(); it != x->end(); ++it)
    {
      private_addprod(A, startRow, startColBis, (*it), y);
      startColBis += (*it)->size();
    }
  }
}

// x and y blocks
void private_prod(SPC::SiconosMatrix A, unsigned int startRow, SPC::SiconosVector x, SP::SiconosVector y, bool init)
{

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

  if (!y->isBlock()) // if y is not a block vector, private_addProd, to consider cases where x is block.
  {
    if (init) // y = subA * x , else y += subA * x
      y->zero();
    private_addprod(A, startRow, 0, x, y);
  }
  else // if y is a block: call private_prod on each block and so on until the considered block is a simple vector.
  {
    unsigned int row = startRow;
    VectorOfVectors::const_iterator it;
    for (it = y->begin(); it != y->end(); ++it)
    {
      private_prod(A, row, x, *it, init);
      row += (*it)->size();
    }
  }
}


// With trans(A) ...
void private_addprod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SP::SiconosVector y)
{
  if (A->isBlock())
    SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for block matrix.");

  if (!x->isBlock()) // if input vector is not block
  {
    // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
    // Then computation of y = subA*x + y.
    unsigned int numA = A->getNum();
    unsigned int numY = y->getNum();
    unsigned int numX = x->getNum();
    unsigned int sizeX = x->size();
    unsigned int sizeY = y->size();

    if (numX != numY)
      SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for x and y of different types.");

    if (numY == 1 && numX == 1)
    {
      if (numA == 1)
        *y->getDensePtr() += prod(ublas::subrange(trans(*A->getDensePtr()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 2)
        *y->getDensePtr() += prod(ublas::subrange(trans(*A->getTriangPtr()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 3)
        *y->getDensePtr() += prod(ublas::subrange(trans(*A->getSymPtr()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 4)
        *y->getDensePtr() += prod(ublas::subrange(trans(*A->getSparsePtr()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else //if(numA==5)
        *y->getDensePtr() += prod(ublas::subrange(trans(*A->getBandedPtr()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
    }
    else // x and y sparse
    {
      if (numA == 4)
        *y->getSparsePtr() += prod(ublas::subrange(trans(*A->getSparsePtr()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getSparsePtr());
      else
        SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for x, y  sparse and A not sparse.");
    }

  }
  else // if block
  {
    VectorOfVectors::const_iterator it;
    unsigned int startColBis = startCol;
    for (it = x->begin(); it != x->end(); ++it)
    {
      private_addprod((*it), A, startRow, startColBis, y);
      startColBis += (*it)->size();
    }
  }
}

void private_prod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::SiconosVector  y, bool init)
{

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between el. of A of index (col) startCol and startCol + sizeY

  if (!y->isBlock()) // if y is not a block vector, private_addProd, to consider cases where x is block.
  {
    if (init) // y = subA * x , else y += subA * x
      y->zero();
    private_addprod(x, A, startCol, 0 , y);
  }
  else // if y is a block: call private_prod on each block and so on until the considered block is a simple vector.
  {
    unsigned int col = startCol;
    VectorOfVectors::const_iterator it;
    for (it = y->begin(); it != y->end(); ++it)
    {
      private_prod(x, A, col, *it, init);
      col += (*it)->size();
    }
  }
}

void private_addprod(double a, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SPC::SiconosVector x, SP::SiconosVector y)
{
  if (A->isBlock())
    SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for block matrix.");

  if (!x->isBlock()) // if input vector is not block
  {
    // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
    // Then computation of y = subA*x + y.
    unsigned int numA = A->getNum();
    unsigned int numY = y->getNum();
    unsigned int numX = x->getNum();
    unsigned int sizeX = x->size();
    unsigned int sizeY = y->size();

    if (numX != numY)
      SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x and y of different types.");

    if (numY == 1 && numX == 1)
    {
      if (numA == 1)
        *y->getDensePtr() += a * prod(ublas::subrange(*A->getDensePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 2)
        *y->getDensePtr() += a * prod(ublas::subrange(*A->getTriangPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 3)
        *y->getDensePtr() += a * prod(ublas::subrange(*A->getSymPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else if (numA == 4)
        *y->getDensePtr() += a * prod(ublas::subrange(*A->getSparsePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
      else //if(numA==5)
        *y->getDensePtr() += a * prod(ublas::subrange(*A->getBandedPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
    }
    else // x and y sparse
    {
      if (numA == 4)
        *y->getSparsePtr() += a * prod(ublas::subrange(*A->getSparsePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getSparsePtr());
      else
        SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x, y  sparse and A not sparse.");
    }

  }
  else // if block
  {
    VectorOfVectors::const_iterator it;
    unsigned int startColBis = startCol;
    for (it = x->begin(); it != x->end(); ++it)
    {
      private_addprod(a, A, startRow, startColBis, (*it), y);
      startColBis += (*it)->size();
    }
  }
}

void private_prod(double a, SPC::SiconosMatrix A, unsigned int startRow, SPC::SiconosVector x, SP::SiconosVector  y, bool init)
{

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

  if (!y->isBlock()) // if y is not a block vector, private_addProd, to consider cases where x is block.
  {
    if (init) // y = subA * x , else y += subA * x
      y->zero();
    private_addprod(a, A, startRow, 0, x, y);
  }
  else // if y is a block: call private_prod on each block and so on until the considered block is a simple vector.
  {
    unsigned int row = startRow;
    VectorOfVectors::const_iterator it;
    for (it = y->begin(); it != y->end(); ++it)
    {
      private_prod(a, A, row, x, *it, init);
      row += (*it)->size();
    }
  }
}

const SimpleVector prod(const SiconosMatrix& A, const SiconosVector& x)
{
  // To compute y = A * x

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: inconsistent sizes.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();

  if (numA == 0) // if A is block ...
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: not yet implemented for block matrix.");

  if (numA == 6) // A = 0
    return (DenseVect)(ublas::zero_vector<double>(x.size()));

  else if (numA == 7) // A = Identity
    return x;

  else
  {
    if (numX != 0) // if x is not a block vector.
    {
      if (numX == 1)
      {
        if (numA == 1)
          return (DenseVect)(prod(*A.getDensePtr(), *x.getDensePtr()));
        else if (numA == 2)
          return (DenseVect)(prod(*A.getTriangPtr(), *x.getDensePtr()));
        else if (numA == 3)
          return (DenseVect)(prod(*A.getSymPtr(), *x.getDensePtr()));
        else if (numA == 4)
          return (DenseVect)(prod(*A.getSparsePtr(), *x.getDensePtr()));
        else // if(numA==5)
          return (DenseVect)(prod(*A.getBandedPtr(), *x.getDensePtr()));
      }
      else //if(numX == 4)
      {
        if (numA == 1)
          return (DenseVect)(prod(*A.getDensePtr(), *x.getSparsePtr()));
        else if (numA == 2)
          return (DenseVect)(prod(*A.getTriangPtr(), *x.getSparsePtr()));
        else if (numA == 3)
          return (DenseVect)(prod(*A.getSymPtr(), *x.getSparsePtr()));
        else if (numA == 4)
          return (DenseVect)(prod(*A.getSparsePtr(), *x.getSparsePtr()));
        else // if(numA==5)
          return (DenseVect)(prod(*A.getBandedPtr(), *x.getSparsePtr()));
      }
    }
    else // if (x.isBlock())
    {
      VectorOfVectors::const_iterator it;
      SimpleVector res(A.size(0));
      unsigned int start = 0;

      for (it = x.begin(); it != x.end(); ++it)
      {
        private_addprod(createSiconosMatrixSPtrConst(A), 0, start, *it, createSiconosVectorSPtr(res));
        start += (*it)->size();
      }
      return res;
    }
  }
}

void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += A*x if init = false.

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

  if (A.size(0) != y.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("prod(A,x,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init)
      y.zero();
    //else nothing
  }

  else if (numA == 7) // A = identity
  {
    if (!init)
      y += x;
    else
    {
      if (&x != &y) y = x ; // if x and y do not share memory (ie are different objects)
      // else nothing
    }
  }

  else // A is not 0 or identity
  {

    // === First case: y is not a block vector ===
    if (numY != 0)
    {
      // if x is block: call of a specific function to treat each block
      if (numX == 0)
      {
        if (init)
          y.zero();
        unsigned int startRow = 0;
        unsigned int startCol = 0;
        // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA x[i]), with subA a submatrix of A,
        // starting from position startRow in rows and startCol in columns.
        // private_prod takes also into account the fact that each block of x can also be a block.
        VectorOfVectors::const_iterator it;
        for (it = x.begin(); it != x.end(); ++it)
        {
          private_addprod(createSiconosMatrixSPtrConst(A), startRow, startCol, *it, createSiconosVectorSPtr(y));
          startCol += (*it)->size();
        }
      }
      else // If neither x nor y are block: direct call to ublas::prod.
      {
        if (init)
        {
          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.getDensePtr()) = ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
                else
                  noalias(*y.getSparsePtr()) = ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              }
              else //if(numA==5)
                noalias(*y.getDensePtr()) = ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.getDensePtr() = ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                *y.getDensePtr() = ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                *y.getDensePtr() = ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                *y.getDensePtr() = ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                *y.getDensePtr() = ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.getSparsePtr() = ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                *y.getSparsePtr() = ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                *y.getSparsePtr() = ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
                *y.getSparsePtr() = ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              else //if(numA==5)
                *y.getSparsePtr() = ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
        }
        else // += case
        {
          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.getDensePtr()) += ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
                else
                  noalias(*y.getSparsePtr()) += ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              }
              else //if(numA==5)
                noalias(*y.getDensePtr()) += ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.getDensePtr() += ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                *y.getDensePtr() += ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                *y.getDensePtr() += ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                *y.getDensePtr() += ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                *y.getDensePtr() += ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.getSparsePtr() += ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                *y.getSparsePtr() += ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                *y.getSparsePtr() += ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
                *y.getSparsePtr() += ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              else //if(numA==5)
                *y.getSparsePtr() += ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
        }
      }
    }
    else // === Second case: y is a block vector ===
    {
      unsigned int startRow = 0;
      VectorOfVectors::const_iterator it;
      // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
      // private_prod takes into account the fact that x and y[i] may be block vectors.
      for (it = y.begin(); it != y.end(); ++it)
      {
        private_prod(createSiconosMatrixSPtrConst(A), startRow, createSiconosVectorSPtrConst(x), *it, init);
        startRow += (*it)->size();
      }
    }
  }
}

void subprod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, const std::vector<unsigned int>& coord, bool init)
{
  // To compute subY = subA * subX in an "optimized" way (in comparison with y = prod(A,x) )
  // or subY += subA*subX if init = false.

  // coord is [r0A r1A c0A c1A r0x r1x r0y r1y]
  //
  // subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and columns between c0A and c1A-1.
  // The same for x and y with rix and riy.

  // Check dims
  unsigned int rowA = coord[1] - coord[0];
  unsigned int colA = coord[3] - coord[2];
  unsigned int dimX = coord[5] - coord[4];
  unsigned int dimY = coord[7] - coord[6];

  if (colA != dimX)
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: inconsistent sizes between A and x.");

  if (rowA != dimY)
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: inconsistent sizes between A and y.");

  if (dimX > x.size() || dimY > y.size() || rowA > A.size(0) || colA > A.size(1))
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: input index too large.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0 || numY == 0) // If A,x or y is Block
    SiconosMatrixException::selfThrow("subprod(A,x,y) error: not yet implemented for A/y block matrices/vector.");

  if (numA == 6) // A = 0
  {
    if (init)
    {
      if (numY == 1)
        ublas::subrange(*y.getDensePtr(), coord[6], coord[7]) *= 0.0;
      else //if(numY==4)
        ublas::subrange(*y.getSparsePtr(), coord[6], coord[7]) *= 0.0;
    }
    //else nothing
  }

  else if (numA == 7) // A = identity
  {
    if (!init)
      ublas::subrange(*y.getDensePtr(), coord[6], coord[7]) += ublas::subrange(*x.getDensePtr(), coord[4], coord[5]);
    else
    {
      // if x and y do not share memory (ie are different objects)
      if (&x != &y)
        noalias(ublas::subrange(*y.getDensePtr(), coord[6], coord[7])) = ublas::subrange(*x.getDensePtr(), coord[4], coord[5]);

      // else nothing
    }
  }

  else // A is not 0 or identity
  {
    if (numX == 0) // ie if x is a block vector
    {
      VectorOfVectors::const_iterator it;
      // Number of the subvector of x that handles element at position coord[4]
      unsigned int firstBlockNum = x.getNumVectorAtPos(coord[4]);
      // Number of the subvector of x that handles element at position coord[5]
      unsigned int lastBlockNum = x.getNumVectorAtPos(coord[5]);
      std::vector<unsigned int> subCoord = coord;
      SPC::SiconosVector  tmp = x[firstBlockNum];
      unsigned int subSize =  x[firstBlockNum]->size(); // Size of the sub-vector
      const Index * xTab = x.getTabIndexPtr();
      if (firstBlockNum != 0)
      {
        subCoord[4] -= (*xTab)[firstBlockNum - 1];
        subCoord[5] =  std::min(coord[5] - (*xTab)[firstBlockNum - 1], subSize);
      }
      else
        subCoord[5] =  std::min(coord[5], subSize);

      if (firstBlockNum == lastBlockNum)
      {
        subprod(A, *tmp, y, subCoord, init);
      }
      else
      {
        unsigned int xPos = 0 ; // Position in x of the current sub-vector of x
        bool firstLoop = true;
        subCoord[3] = coord[2] + subCoord[5] - subCoord[4];
        for (it = x.begin(); it != x.end(); ++it)
        {
          if ((*it)->getNum() == 0)
            SiconosMatrixException::selfThrow("subprod(A,x,y) error: not yet implemented for x block of blocks ...");
          if (xPos >= firstBlockNum && xPos <= lastBlockNum)
          {
            tmp = x[xPos];
            if (firstLoop)
            {
              subprod(A, *tmp, y, subCoord, init);
              firstLoop = false;
            }
            else
            {
              subCoord[2] += subCoord[5] - subCoord[4]; // !! old values for 4 and 5
              subSize = tmp->size();
              subCoord[4] = 0;
              subCoord[5] = std::min(coord[5] - (*xTab)[xPos - 1], subSize);
              subCoord[3] = subCoord[2] + subCoord[5] - subCoord[4];
              subprod(A, *tmp, y, subCoord, false);
            }
          }
          xPos++;
        }
      }
    }
    else
    {
      if (init)
      {
        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subX(*x.getDensePtr(), ublas::range(coord[4], coord[5]));

            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");
            ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));

            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subX(*x.getSparsePtr(), ublas::range(coord[4], coord[5]));
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));

              if (numY == 1)
              {
                ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
                noalias(subY) = ublas::prod(subA, subX);
              }
              else
              {
                ublas::vector_range<SparseVect> subY(*y.getSparsePtr(), ublas::range(coord[6], coord[7]));
                noalias(subY) = ublas::prod(subA, subX);
              }
#endif
            }
            else //if(numA==5)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) = ublas::prod(subA, subX);
            }
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> and vector_range<SparseVect> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subY(*y.getSparsePtr(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY = ublas::prod(subA, subY);
            }
          }
        }
      }
      else // += case
      {
        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subX(*x.getDensePtr(), ublas::range(coord[4], coord[5]));

            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");
            ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));

            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subX(*x.getSparsePtr(), ublas::range(coord[4], coord[5]));
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 2)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 3)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              if (numY == 1)
              {
                ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
                noalias(subY) += ublas::prod(subA, subX);
              }
              else
              {
                ublas::vector_range<SparseVect> subY(*y.getSparsePtr(), ublas::range(coord[6], coord[7]));
                noalias(subY) += ublas::prod(subA, subX);
              }
#endif
            }
            else //if(numA==5)
            {
              ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[6], coord[7]));
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              noalias(subY) += ublas::prod(subA, subX);
            }
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            ublas::vector_range<DenseVect> subY(*y.getDensePtr(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
          }
          else //if(numX == 4)
          {
            ublas::vector_range<SparseVect> subY(*y.getSparsePtr(), ublas::range(coord[4], coord[5]));
            if (numA == 1)
            {
              ublas::matrix_range<DenseMat> subA(*A.getDensePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 2)
            {
              ublas::matrix_range<TriangMat> subA(*A.getTriangPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 3)
            {
              ublas::matrix_range<SymMat> subA(*A.getSymPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
            else if (numA == 4)
            {
#ifdef BOOST_LIMITATION
              SiconosMatrixException("SimpleMatrix::subprod warning - ublas::matrix_range<SparseMat> does not exist for MacOs.");
#else
              ublas::matrix_range<SparseMat> subA(*A.getSparsePtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
#endif
            }
            else //if(numA==5)
            {
              ublas::matrix_range<BandedMat> subA(*A.getBandedPtr(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
              subY += ublas::prod(subA, subY);
            }
          }
        }
      }
    }
  }
}

void prod(double a, const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = a*A * x in an "optimized" way (in comparison with y = prod(A,x) )
  // or y += a*A*x if init = false.

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

  if (A.size(0) != y.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("prod(A,x,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init)
      y.zero();
    //else nothing
  }

  else if (numA == 7) // A = identity
  {
    scal(a, x, y, init);
  }

  else // A is not 0 or identity
  {

    // === First case: y is not a block vector ===
    if (numY != 0)
    {
      // if x is block: call of a specific function to treat each block
      if (numX == 0)
      {
        if (init)
          y.zero();
        unsigned int startRow = 0;
        unsigned int startCol = 0;
        // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA x[i]), with subA a submatrix of A,
        // starting from position startRow in rows and startCol in columns.
        // private_prod takes also into account the fact that each block of x can also be a block.
        VectorOfVectors::const_iterator it;
        for (it = x.begin(); it != x.end(); ++it)
        {
          private_addprod(a, createSiconosMatrixSPtrConst(A), startRow, startCol, *it, createSiconosVectorSPtr(y));
          startCol += (*it)->size();
        }
      }
      else // If neither x nor y are block: direct call to ublas::prod.
      {
        if (init)
        {
          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.getDensePtr()) = a * ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
                else
                  noalias(*y.getSparsePtr()) = a * ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              }
              else //if(numA==5)
                noalias(*y.getDensePtr()) = a * ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.getDensePtr() = a * ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                *y.getDensePtr() = a * ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                *y.getDensePtr() = a * ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                *y.getDensePtr() = a * ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                *y.getDensePtr() = a * ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.getSparsePtr() = a * ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                *y.getSparsePtr() = a * ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                *y.getSparsePtr() = a * ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
                *y.getSparsePtr() = a * ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              else //if(numA==5)
                *y.getSparsePtr() = a * ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
        }
        else // += case
        {
          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.getDensePtr()) += a * ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
                else
                  noalias(*y.getSparsePtr()) += a * ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              }
              else //if(numA==5)
                noalias(*y.getDensePtr()) += a * ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.getDensePtr() += a * ublas::prod(*A.getDensePtr(), *x.getDensePtr());
              else if (numA == 2)
                *y.getDensePtr() += a * ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
              else if (numA == 3)
                *y.getDensePtr() += a * ublas::prod(*A.getSymPtr(), *x.getDensePtr());
              else if (numA == 4)
                *y.getDensePtr() += a * ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
              else //if(numA==5)
                *y.getDensePtr() += a * ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.getSparsePtr() += a * ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
              else if (numA == 2)
                *y.getSparsePtr() += a * ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
              else if (numA == 3)
                *y.getSparsePtr() += a * ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
              else if (numA == 4)
                *y.getSparsePtr() += a * ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
              else //if(numA==5)
                *y.getSparsePtr() += a * ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
            }
          }
        }
      }
    }
    else // === Second case: y is a block vector ===
    {
      unsigned int startRow = 0;
      VectorOfVectors::const_iterator it;
      // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
      // private_prod takes into account the fact that x and y[i] may be block vectors.
      for (it = y.begin(); it != y.end(); ++it)
      {
        private_prod(a, createSiconosMatrixSPtrConst(A), startRow, createSiconosVectorSPtrConst(x), *it, init);
        startRow += (*it)->size();
      }
    }
  }
}

void prod(const SiconosVector& x, const SiconosMatrix& A, SiconosVector& y, bool init)
{
  // To compute y = trans(A) * x in an "optimized" way, if init = true
  // (or y = trans(A) * x + y if init = false

  if (A.size(0) != x.size())
    SiconosMatrixException::selfThrow("prod(x,A,y) error: inconsistent sizes between A and x.");

  if (A.size(1) != y.size())
    SiconosMatrixException::selfThrow("prod(x,A,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("prod(x,A,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init)
      y.zero();
    // else nothing
  }

  else if (numA == 7) // A = identity
  {
    if (!init)
      y += x;
    else
    {
      if (&x != &y) y = x ; // if x and y do not share memory (ie are different objects)
      // else nothing
    }
  }

  else // A is not 0 or identity
  {
    // === First case: y is not a block vector ===
    if (numY != 0)
    {
      // if x is block: call of a specific function to treat each block
      if (numX == 0)
      {
        if (init)
          y.zero();
        unsigned int startRow = 0;
        unsigned int startCol = 0;
        // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA x[i]), with subA a submatrix of A,
        // starting from position startRow in rows and startCol in columns.
        // private_prod takes also into account the fact that each block of x can also be a block.
        VectorOfVectors::const_iterator it;
        for (it = x.begin(); it != x.end(); ++it)
        {
          private_addprod(*it, createSiconosMatrixSPtrConst(A), startRow, startCol, createSiconosVectorSPtr(y));
          startCol += (*it)->size();
        }
      }
      else // If neither x nor y are block: direct call to ublas::prod.
      {
        if (init)
        {

          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getDensePtr()), *x.getDensePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getTriangPtr()), *x.getDensePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getSymPtr()), *x.getDensePtr());
              else if (numA == 4)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getSparsePtr()), *x.getDensePtr());
              else //if(numA==5)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getBandedPtr()), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getDensePtr()), *x.getSparsePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getTriangPtr()), *x.getSparsePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getSymPtr()), *x.getSparsePtr());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getSparsePtr()), *x.getSparsePtr());
                else
                  noalias(*y.getSparsePtr()) = ublas::prod(trans(*A.getSparsePtr()), *x.getSparsePtr());
              }
              else //if(numA==5)
                noalias(*y.getDensePtr()) = ublas::prod(trans(*A.getBandedPtr()), *x.getSparsePtr());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.getDensePtr() = ublas::prod(trans(*A.getDensePtr()), *x.getDensePtr());
              else if (numA == 2)
                *y.getDensePtr() = ublas::prod(trans(*A.getTriangPtr()), *x.getDensePtr());
              else if (numA == 3)
                *y.getDensePtr() = ublas::prod(trans(*A.getSymPtr()), *x.getDensePtr());
              else if (numA == 4)
                *y.getDensePtr() = ublas::prod(trans(*A.getSparsePtr()), *x.getDensePtr());
              else //if(numA==5)
                *y.getDensePtr() = ublas::prod(trans(*A.getBandedPtr()), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.getSparsePtr() = ublas::prod(trans(*A.getDensePtr()), *x.getSparsePtr());
              else if (numA == 2)
                *y.getSparsePtr() = ublas::prod(trans(*A.getTriangPtr()), *x.getSparsePtr());
              else if (numA == 3)
                *y.getSparsePtr() = ublas::prod(trans(*A.getSymPtr()), *x.getSparsePtr());
              else if (numA == 4)
                *y.getSparsePtr() = ublas::prod(trans(*A.getSparsePtr()), *x.getSparsePtr());
              else //if(numA==5)
                *y.getSparsePtr() = ublas::prod(trans(*A.getBandedPtr()), *x.getSparsePtr());
            }
          }
        }
        else // += case
        {

          if (&x != &y) // if no common memory between x and y.
          {
            if (numX == 1)
            {
              if (numY != 1)
                SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getDensePtr()), *x.getDensePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getTriangPtr()), *x.getDensePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getSymPtr()), *x.getDensePtr());
              else if (numA == 4)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getSparsePtr()), *x.getDensePtr());
              else //if(numA==5)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getBandedPtr()), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numY != 1 && numA != 4)
                SiconosMatrixException::selfThrow("prod(x,A,y) error: y (output) must be a dense vector.");

              if (numA == 1)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getDensePtr()), *x.getSparsePtr());
              else if (numA == 2)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getTriangPtr()), *x.getSparsePtr());
              else if (numA == 3)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getSymPtr()), *x.getSparsePtr());
              else if (numA == 4)
              {
                if (numY == 1)
                  noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getSparsePtr()), *x.getSparsePtr());
                else
                  noalias(*y.getSparsePtr()) += ublas::prod(trans(*A.getSparsePtr()), *x.getSparsePtr());
              }
              else //if(numA==5)
                noalias(*y.getDensePtr()) += ublas::prod(trans(*A.getBandedPtr()), *x.getSparsePtr());
            }
          }
          else // if x and y are the same object => alias
          {
            if (numX == 1)
            {
              if (numA == 1)
                *y.getDensePtr() += ublas::prod(trans(*A.getDensePtr()), *x.getDensePtr());
              else if (numA == 2)
                *y.getDensePtr() += ublas::prod(trans(*A.getTriangPtr()), *x.getDensePtr());
              else if (numA == 3)
                *y.getDensePtr() += ublas::prod(trans(*A.getSymPtr()), *x.getDensePtr());
              else if (numA == 4)
                *y.getDensePtr() += ublas::prod(trans(*A.getSparsePtr()), *x.getDensePtr());
              else //if(numA==5)
                *y.getDensePtr() += ublas::prod(trans(*A.getBandedPtr()), *x.getDensePtr());
            }
            else //if(numX == 4)
            {
              if (numA == 1)
                *y.getSparsePtr() += ublas::prod(trans(*A.getDensePtr()), *x.getSparsePtr());
              else if (numA == 2)
                *y.getSparsePtr() += ublas::prod(trans(*A.getTriangPtr()), *x.getSparsePtr());
              else if (numA == 3)
                *y.getSparsePtr() += ublas::prod(trans(*A.getSymPtr()), *x.getSparsePtr());
              else if (numA == 4)
                *y.getSparsePtr() += ublas::prod(trans(*A.getSparsePtr()), *x.getSparsePtr());
              else //if(numA==5)
                *y.getSparsePtr() += ublas::prod(trans(*A.getBandedPtr()), *x.getSparsePtr());
            }
          }
        }
      }
    }
    else // === Second case: y is a block vector ===
    {
      unsigned int startRow = 0;
      VectorOfVectors::const_iterator it;
      // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
      // private_prod takes into account the fact that x and y[i] may be block vectors.
      for (it = y.begin(); it != y.end(); ++it)
      {
        private_prod(createSiconosVectorSPtrConst(x), createSiconosMatrixSPtrConst(A), startRow, *it, init);
        startRow += (*it)->size();
      }
    }
  }
}

void axpy_prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init)
{
  // To compute y = A * x ( init = true) or y += A * x (init = false) using ublas::axpy_prod

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

  if (A.size(0) != y.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numA == 0) // If A is Block
    SiconosMatrixException::selfThrow("axpy_prod(A,x,y) error: not yet implemented for block matrices.");

  if (numA == 6) // A = 0
  {
    if (init) y.zero(); // else nothing ...
  }

  else if (numA == 7) // A = identity
  {
    if (!init) y += x;
    else
    {
      if (&x != &y)
        y = x ; // if x and y do not share memory (ie are different objects)
    }
    // else nothing
  }

  else // A is not 0 or identity
  {
    // === First case: y is not a block vector ===
    if (numY != 0)
    {
      // if x is block: call of a specific function to treat each block
      if (numX == 0)
      {
        if (init) y.zero();
        unsigned int startRow = 0;
        unsigned int startCol = 0;
        // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA x[i]), with subA a submatrix of A,
        // starting from position startRow in rows and startCol in columns.
        // private_prod takes also into account the fact that each block of x can also be a block.
        VectorOfVectors::const_iterator it;
        for (it = x.begin(); it != x.end(); ++it)
        {
          private_addprod(createSiconosMatrixSPtrConst(A), startRow, startCol, *it, createSiconosVectorSPtr(y));
          startCol += (*it)->size();
        }
      }
      else // If neither x nor y are block: direct call to ublas::prod.
      {
        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            if (numY != 1)
              SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              ublas::axpy_prod(*A.getDensePtr(), *x.getDensePtr(), *y.getDensePtr(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.getTriangPtr(), *x.getDensePtr(), *y.getDensePtr(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.getSymPtr(), *x.getDensePtr(), *y.getDensePtr(), init);
            else if (numA == 4)
              ublas::axpy_prod(*A.getSparsePtr(), *x.getDensePtr(), *y.getDensePtr(), init);
            else //if(numA==5)
              ublas::axpy_prod(*A.getBandedPtr(), *x.getDensePtr(), *y.getDensePtr(), init);
          }
          else //if(numX == 4)
          {
            if (numY != 1 && numA != 4)
              SiconosMatrixException::selfThrow("axpy_prod(A,x,y) error: y (output) must be a dense vector.");

            if (numA == 1)
              ublas::axpy_prod(*A.getDensePtr(), *x.getSparsePtr(), *y.getDensePtr(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.getTriangPtr(), *x.getSparsePtr(), *y.getDensePtr(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.getSymPtr(), *x.getSparsePtr(), *y.getDensePtr(), init);
            else if (numA == 4)
            {
              if (numY == 1)
                ublas::axpy_prod(*A.getSparsePtr(), *x.getSparsePtr(), *y.getDensePtr(), init);
              else
                ublas::axpy_prod(*A.getSparsePtr(), *x.getSparsePtr(), *y.getSparsePtr(), init);
            }
            else //if(numA==5)
              ublas::axpy_prod(*A.getBandedPtr(), *x.getSparsePtr(), *y.getDensePtr(), init);
          }
        }
        else // if x and y are the same object => alias
        {
          if (numX == 1)
          {
            if (numA == 1)
              ublas::axpy_prod(*A.getDensePtr(), *x.getDensePtr(), *x.getDensePtr(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.getTriangPtr(), *x.getDensePtr(), *x.getDensePtr(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.getSymPtr(), *x.getDensePtr(), *x.getDensePtr(), init);
            else if (numA == 4)
              ublas::axpy_prod(*A.getSparsePtr(), *x.getDensePtr(), *x.getDensePtr(), init);
            else //if(numA==5)
              ublas::axpy_prod(*A.getBandedPtr(), *x.getDensePtr(), *x.getDensePtr(), init);
          }
          else //if(numX == 4)
          {
            if (numA == 1)
              ublas::axpy_prod(*A.getDensePtr(), *x.getSparsePtr(), *x.getSparsePtr(), init);
            else if (numA == 2)
              ublas::axpy_prod(*A.getTriangPtr(), *x.getSparsePtr(), *x.getSparsePtr(), init);
            else if (numA == 3)
              ublas::axpy_prod(*A.getSymPtr(), *x.getSparsePtr(), *x.getSparsePtr(), init);
            else if (numA == 4)
              ublas::axpy_prod(*A.getSparsePtr(), *x.getSparsePtr(), *x.getSparsePtr(), init);
            else //if(numA==5)
              ublas::axpy_prod(*A.getBandedPtr(), *x.getSparsePtr(), *x.getSparsePtr(), init);
          }
        }
      }
    }
    else // === Second case: y is a block vector ===
    {
      unsigned int startRow = 0;
      VectorOfVectors::const_iterator it;
      // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
      // private_prod takes into account the fact that x and y[i] may be block vectors.
      for (it = y.begin(); it != y.end(); ++it)
      {
        private_prod(createSiconosMatrixSPtrConst(A), startRow, createSiconosVectorSPtrConst(x), *it, init);
        startRow += (*it)->size();
      }
    }
  }
}

void gemv(const CBLAS_TRANSPOSE transA, double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock() || A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  atlas::gemv(transA, a, *A.getDensePtr(), *x.getDensePtr(), b, *y.getDensePtr());
}

void gemv(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock() || A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  atlas::gemv(a, *A.getDensePtr(), *x.getDensePtr(), b, *y.getDensePtr());
}

void gemv(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock() || A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  atlas::gemv(*A.getDensePtr(), *x.getDensePtr(), *y.getDensePtr());
}

void gemm(const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  if (A.isBlock() || B.isBlock() || C.isBlock())
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA != 1 || numB != 1 || numC != 1)
    SiconosMatrixException::selfThrow("gemm(...) failed: reserved to dense matrices.");

  atlas::gemm(transA, transB, a, *A.getDensePtr(), *B.getDensePtr(), b, *C.getDensePtr());
  C.resetLU();
}

void gemm(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();

  if (numA == 0 || numB == 0 || numC == 0)
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");

  if (numA == 1 && numB == 1 && numC == 1)
    atlas::gemm(a, *A.getDensePtr(), *B.getDensePtr(), b, *C.getDensePtr());
  else if (numA == 1 && numB == 1 && numC != 1)
  {
    // To be improved ...

    DenseMat * tmpA = NULL;
    DenseMat * tmpB = NULL;
    DenseMat * tmpC = NULL;

    if (numA != 1)
      tmpA = new DenseMat(*A.getDensePtr());
    else
      tmpA = A.getDensePtr();

    if (numB != 1)
      tmpB = new DenseMat(*B.getDensePtr());
    else
      tmpB = B.getDensePtr();

    if (numC != 1)
      tmpC = new DenseMat(*C.getDensePtr());
    else
      tmpC = C.getDensePtr();

    atlas::gemm(a, *tmpA, *tmpB, b, *tmpC);
    if (numC != 1)
    {
      noalias(*C.getDensePtr()) = *tmpC;
      delete tmpC;
    }

    if (numA != 1)
      delete tmpA;
    if (numB != 1)
      delete tmpB;
  }
  C.resetLU();
}

void gemm(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C)
{
  if (A.isBlock() || B.isBlock() || C.isBlock())
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA != 1 || numB != 1 || numC != 1)
    SiconosMatrixException::selfThrow("gemm(...) failed: reserved to dense matrices.");

  atlas::gemm(*A.getDensePtr(), *B.getDensePtr(), *C.getDensePtr());
  C.resetLU();
}

void scal(double a, const SiconosMatrix& A, SiconosMatrix& B, bool init)
{
  // To compute B = a * A (init = true) or B += a*A (init = false).

  if (&A == &B)
  {
    if (init) B *= a;
    else B *= (1.0 + a);
  }
  else
  {
    unsigned int numA = A.getNum();
    unsigned int numB = B.getNum();

    if (numB == 6 || numB == 7) // B = 0 or identity.
      SiconosMatrixException::selfThrow("scal(a,A,B) : forbidden for B being a zero or identity matrix.");

    if (numA == 6)
    {
      if (init) B.zero(); // else nothing
    }
    else if (numA == 7)
    {
      if (init)
      {
        B.eye();
        B *= a;
      }
      else
      {
        // Assuming B is square ...
        for (unsigned int i = 0; i < B.size(0); ++i)
          B(i, i) += a;
      }
    }
    else
    {
      if (numA == numB) // if A and B are of the same type ...
      {
        switch (numA)
        {

        case 0: // A and B are block
          if (isComparableTo(&A, &B))
          {
            BlockIterator1 itB1;
            BlockIterator2 itB2;
            ConstBlockIterator1 itA1 = A.begin();
            ConstBlockIterator2 itA2;
            for (itB1 = B.begin(); itB1 != B.end(); ++itB1)
            {
              itA2 = itA1.begin();
              for (itB2 = itB1.begin(); itB2 != itB1.end(); ++itB2)
              {
                scal(a, **itA2++, **itB2, init);
              }
              itA1++;
            }
          }
          else // if A and B are not "block-consistent"
          {
            if (init)
            {
              for (unsigned int i = 0; i < A.size(0); ++i)
                for (unsigned int j = 0; j < A.size(1); ++j)
                  B(i, j) = a * A(i, j);
            }
            else
            {
              for (unsigned int i = 0; i < A.size(0); ++i)
                for (unsigned int j = 0; j < A.size(1); ++j)
                  B(i, j) += a * A(i, j);
            }
          }
          break;

        case 1: // if both are dense
          if (init)
            noalias(*B.getDensePtr()) = a ** A.getDensePtr();
          else
            noalias(*B.getDensePtr()) += a ** A.getDensePtr();
          break;
        case 2:
          if (init)
            noalias(*B.getTriangPtr()) = a ** A.getTriangPtr();
          else
            noalias(*B.getTriangPtr()) += a ** A.getTriangPtr();
          break;
        case 3:
          if (init)
            noalias(*B.getSymPtr()) = a ** A.getSymPtr();
          else
            noalias(*B.getSymPtr()) += a ** A.getSymPtr();
          break;
        case 4:
          if (init)
            noalias(*B.getSparsePtr()) = a ** A.getSparsePtr();
          else
            noalias(*B.getSparsePtr()) += a ** A.getSparsePtr();
          break;
        case 5:
          if (init)
            noalias(*B.getBandedPtr()) = a ** A.getBandedPtr();
          else
            noalias(*B.getBandedPtr()) += a ** A.getBandedPtr();
          break;
        }
      }
      else // if A and B are of different types.
      {
        if (numA == 0 || numB == 0) // if A or B is block
        {
          if (init)
          {
            B = A;
            B *= a;
          }
          else
          {
            SimpleMatrix tmp(A);
            tmp *= a;
            B += tmp; // bof bof ...
          }
        }
        else
        {
          if (numB != 1)
            SiconosMatrixException::selfThrow("scal(a,A,B) failed. A and B types do not fit together.");

          if (init)
          {
            switch (numB)
            {
            case 1:
              noalias(*B.getDensePtr()) = a ** A.getDensePtr();
              break;
            case 2:
              noalias(*B.getDensePtr()) = a ** A.getTriangPtr();
              break;
            case 3:
              noalias(*B.getDensePtr()) = a ** A.getSymPtr();
              break;
            case 4:
              noalias(*B.getDensePtr()) = a ** A.getSparsePtr();
              break;
            case 5:
              noalias(*B.getDensePtr()) = a ** A.getBandedPtr();
              break;
            }
          }
          else

          {
            switch (numB)
            {
            case 1:
              noalias(*B.getDensePtr()) += a ** A.getDensePtr();
              break;
            case 2:
              noalias(*B.getDensePtr()) += a ** A.getTriangPtr();
              break;
            case 3:
              noalias(*B.getDensePtr()) += a ** A.getSymPtr();
              break;
            case 4:
              noalias(*B.getDensePtr()) += a ** A.getSparsePtr();
              break;
            case 5:
              noalias(*B.getDensePtr()) += a ** A.getBandedPtr();
              break;
            }
          }
        }
      }
    }
  }
}
