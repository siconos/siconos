/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "SiconosConfig.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>

#include "ioMatrix.hpp"
#include "BlockVector.hpp"
#include "BlockMatrixIterators.hpp"
#include "BlockMatrix.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

#include "SiconosAlgebra.hpp"
// Useful function (print ...) from boost bindings examples.
#include "bindings_utils.hpp"

#include "Tools.hpp"

using namespace Siconos;

#include <boost/numeric/bindings/blas.hpp>
namespace siconosBindings = boost::numeric::bindings::blas;


// =================================================
//                CONSTRUCTORS
// =================================================

using std::cout;
using std::endl;


// Default (protected, used only for derived classes)
SimpleMatrix::SimpleMatrix(int i): SiconosMatrix(1), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Dense = new DenseMat(ublas::zero_matrix<double>());
}

SimpleMatrix::SimpleMatrix(): SiconosMatrix(1), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Dense = new DenseMat(ublas::zero_matrix<double>());
}

// parameters: dimensions and type.
SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
  SiconosMatrix(1), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  if (typ == DENSE)
  {
    mat.Dense = new DenseMat(ublas::zero_matrix<double>(row, col));
    // _num = 1; default value
  }
  else if (typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat(ublas::zero_matrix<double>(row, col));
    _num = 2;
  }
  else if (typ == SYMMETRIC)
  {
    mat.Sym = new SymMat(ublas::zero_matrix<double>(row, col));
    _num = 3;
  }
  else if (typ == SPARSE)
  {
    mat.Sparse = new SparseMat(row, col, upper);
    _num = 4;
    zero();
  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat(row, col, upper, lower);
    _num = 5;
    zero();
  }
  else if (typ == ZERO)
  {
    mat.Zero = new ZeroMat(row, col);
    _num = 6;
  }
  else if (typ == IDENTITY)
  {
    mat.Identity = new IdentityMat(row, col);
    _num = 7;
  }
  else
    SiconosMatrixException::selfThrow("SiconosMatrix::constructor(UBLAS_TYPE type, unsigned int row, unsigned int col): invalid type.");
}

// parameters: dimensions, input value and type
SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, double inputValue, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
  SiconosMatrix(1), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  // This constructor has sense only for dense matrices ...
  if (typ == DENSE)
  {
    mat.Dense = new DenseMat(ublas::scalar_matrix<double>(row, col, inputValue));
    // _num = 1; default value
  }
  else
    SiconosMatrixException::selfThrow("SiconosMatrix::constructor(UBLAS_TYPE type, unsigned int row, unsigned int col, double fillInValue): invalid type.");
}

// // parameters: a vector (stl) of double and the type.
// SimpleMatrix::SimpleMatrix(const std::vector<double>& v, unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int lower, unsigned int upper):
//   SiconosMatrix(1, row, col), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
// {
//   if( (  (v.size() != row*col) && (typ != SYMMETRIC && typ != BANDED) )
//       || (v.size() != row*row && typ == SYMMETRIC)
//       || (typ == BANDED && ( (v.size()) != (unsigned int)(std::max)(row, col)*(lower+1+upper) ) ))
//     SiconosMatrixException::selfThrow("constructor(UBLAS_TYPE, const std::vector<double>, int, int) : invalid vector size");

//   if(typ == DENSE)
//     {
//       mat.Dense = new DenseMat(row,col);
//       // _num = 1; default value
//     }
//   else if(typ == TRIANGULAR)
//     {
//       mat.Triang = new TriangMat(row,col);
//       _num = 2;
//     }
//   else if(typ == SYMMETRIC)
//     {
//       mat.Sym = new SymMat(row);
//       _num = 3;
//     }
//   else if(typ == SPARSE)
//     {
//       SiconosMatrixException::selfThrow("SimpleMatrix::constructor(UBLAS_TYPE, const std::vector<double>, int row, int col, int lower, int upper) : warning -- use constructor(const SparseMat &m) or constructor(UBLAS_TYPE, int row, int col) with UBLAS_TYPE = SPARSE");

//     }
//   else if(typ == BANDED)
//     {
//       mat.Banded = new BandedMat(row, col, lower, upper);
//       _num = 5;
//     }
//   else
//     SiconosMatrixException::selfThrow("constructor(UBLAS_TYPE, const std::vector<double>, int, int) : invalid type of matrix given");

//   std::copy(v.begin(), v.end(), (vect.Dense)->begin());


// }

// Copy constructors
SimpleMatrix::SimpleMatrix(const SimpleMatrix &smat): SiconosMatrix(smat.num()), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  if (_num == 1)
  {
    mat.Dense = new DenseMat(smat.size(0), smat.size(1));
    noalias(*mat.Dense) = (*smat.dense());
  }
  //   mat.Dense = new DenseMat(*smat.dense());

  else if (_num == 2)
    mat.Triang = new TriangMat(*smat.triang());

  else if (_num == 3)

    mat.Sym = new SymMat(*smat.sym());

  else if (_num == 4)
    mat.Sparse = new SparseMat(*smat.sparse());

  else if (_num == 5)
    mat.Banded = new BandedMat(*smat.banded());

  else if (_num == 6)
    mat.Zero = new ZeroMat(smat.size(0), smat.size(1));

  else// if(_num == 7)
    mat.Identity = new IdentityMat(smat.size(0), smat.size(1));
}

/** copy constructor of a block given by the coord = [r0A r1A c0A c1A]
 *  \param A the matrix for extracting the block
 */
SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ):  SiconosMatrix(A.num()), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  if (coord[0]>=coord[1])
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). Empty row range coord[0]>= coord[1]");
  if (coord[2]>=coord[3])
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). Empty column range coord[2]>= coord[3]");
  if (coord[1] > A.size(0) )
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). row index too large.");
  if (coord[3] > A.size(1) )
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). column index too large.");

  if (_num== 1)
  {
    ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Dense=new DenseMat(subA);
  }
  else if (_num == 2)
  {
    ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Triang=new TriangMat(subA);
  }
  else if (_num == 3)
  {
    ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Sym=new SymMat(subA);
  }
  else if (_num == 4)
  {
    ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Sparse=new SparseMat(subA);
  }
  else if (_num == 5)
  {
    ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Banded=new BandedMat(subA);
  }
  else if (_num == 6)
  {
    mat.Zero = new ZeroMat(coord[1]-coord[0], coord[3]-coord[2]);
  }
  else// if(_num == 7)
    mat.Identity = new IdentityMat(coord[1]-coord[0], coord[3]-coord[2] );
}




SimpleMatrix::SimpleMatrix(const SiconosMatrix &m): SiconosMatrix(m.num()), _isPLUFactorized(), _isQRFactorized(false), _isPLUInversed(false)
{
  // _num is set in SiconosMatrix constructor with m.num() ... must be changed if m is Block
  unsigned int numM = m.num();


  _isPLUFactorized= m.isPLUFactorized();
  _isPLUInversed= m.isPLUInversed();

  if (m.ipiv())
    _ipiv.reset(new VInt(*(m.ipiv())));

  if (numM == 0) // ie if m is Block, this matrix is set to a dense.
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    _num = 1;
    // get number of blocks in a row/col of m.
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    ConstBlocksIterator1 it;
    ConstBlocksIterator2 it2;
    unsigned int posRow = 0;
    unsigned int posCol = 0;

    for (it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
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
  else if (_num == 1)
  {
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    noalias(*mat.Dense) = (*m.dense());
  }

  else if (_num == 2)
    mat.Triang = new TriangMat(*m.triang());

  else if (_num == 3)
    mat.Sym = new SymMat(*m.sym());

  else if (_num == 4)
    mat.Sparse = new SparseMat(*m.sparse());

  else if (_num == 5)
    mat.Banded = new BandedMat(*m.banded());

  else if (_num == 6)
    mat.Zero = new ZeroMat(m.size(0), m.size(1));

  else // if(_num == 7)
    mat.Identity = new IdentityMat(m.size(0), m.size(1));
}

SimpleMatrix::SimpleMatrix(const DenseMat& m): SiconosMatrix(1), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Dense = new DenseMat(m);
}

SimpleMatrix::SimpleMatrix(const TriangMat& m): SiconosMatrix(2), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Triang = new TriangMat(m);
}

SimpleMatrix::SimpleMatrix(const SymMat& m): SiconosMatrix(3), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Sym = new SymMat(m);
}

SimpleMatrix::SimpleMatrix(const SparseMat& m): SiconosMatrix(4), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Sparse = new SparseMat(m);
}

SimpleMatrix::SimpleMatrix(const BandedMat& m): SiconosMatrix(5), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Banded = new BandedMat(m);
}

SimpleMatrix::SimpleMatrix(const ZeroMat& m): SiconosMatrix(6), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Zero = new ZeroMat(m);
}

SimpleMatrix::SimpleMatrix(const IdentityMat& m): SiconosMatrix(7), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Identity = new IdentityMat(m);
}

SimpleMatrix::SimpleMatrix(const std::string &file, bool ascii): SiconosMatrix(1), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false)
{
  mat.Dense = new DenseMat();
  if (ascii)
  {
    ioMatrix::read(file, "ascii", *this);
  }
  else
  {
    ioMatrix::read(file, "binary", *this);
  }
}

SimpleMatrix::~SimpleMatrix()
{
  if (_num == 1)
    delete(mat.Dense);
  else if (_num == 2)
    delete(mat.Triang);
  else if (_num == 3)
    delete(mat.Sym);
  else if (_num == 4)
    delete(mat.Sparse);
  else if (_num == 5)
    delete(mat.Banded);
  else if (_num == 6)
    delete(mat.Zero);
  else if (_num == 7)
    delete(mat.Identity);
}

bool SimpleMatrix::isSymmetric(double tol) const
{
  SP::SimpleMatrix  m_trans (new SimpleMatrix(*this));
  m_trans->trans();
  double err = (*this-*m_trans).normInf();
  if ((*m_trans).normInf() > 0.0 )
  {
    err /= (*m_trans).normInf();
  }
  // std::cout << "err_rel  ="<< err <<std::endl;
  return (err < tol);
}
//======================================
// get Ublas component (dense, sym ...)
//======================================

const DenseMat SimpleMatrix::getDense(unsigned int, unsigned int) const
{
  if (_num != 1)
    SiconosMatrixException::selfThrow("SimpleMatrix::getDense(): the current matrix is not a Dense matrix");

  return *mat.Dense;
}

const TriangMat SimpleMatrix::getTriang(unsigned int, unsigned int) const
{
  if (_num != 2)
    SiconosMatrixException::selfThrow("TriangMat SimpleMatrix::getTriang(): the current matrix is not a Triangular matrix");

  return *mat.Triang;
}

const SymMat SimpleMatrix::getSym(unsigned int, unsigned int) const
{
  if (_num != 3)
    SiconosMatrixException::selfThrow("SymMat SimpleMatrix::getSym(): the current matrix is not a Symmetric matrix");

  return *mat.Sym;
}

const SparseMat SimpleMatrix::getSparse(unsigned int, unsigned int) const
{
  if (_num != 4)
    SiconosMatrixException::selfThrow("SparseMat SimpleMatrix::getSparse(): the current matrix is not a Sparse matrix");

  return *mat.Sparse;
}

const BandedMat SimpleMatrix::getBanded(unsigned int, unsigned int) const
{
  if (_num != 5)
    SiconosMatrixException::selfThrow("BandedMat SimpleMatrix::getBanded(): the current matrix is not a Banded matrix");

  return *mat.Banded;
}

const ZeroMat SimpleMatrix::getZero(unsigned int, unsigned int) const
{
  if (_num != 6)
    SiconosMatrixException::selfThrow("ZeroMat SimpleMatrix::getZero(): the current matrix is not a Zero matrix");

  return *mat.Zero;
}

const IdentityMat SimpleMatrix::getIdentity(unsigned int, unsigned int) const
{
  if (_num != 7)
    SiconosMatrixException::selfThrow("IdentityMat SimpleMatrix::getIdentity(): the current matrix is not a Identity matrix");

  return *mat.Identity;
}

DenseMat* SimpleMatrix::dense(unsigned int, unsigned int) const
{
  if (_num != 1)
    SiconosMatrixException::selfThrow("DenseMat* SimpleMatrix::dense(): the current matrix is not a Dense matrix");

  return mat.Dense;
}

TriangMat* SimpleMatrix::triang(unsigned int, unsigned int) const
{
  if (_num != 2)
    SiconosMatrixException::selfThrow("TriangMat* SimpleMatrix::triang(): the current matrix is not a Triangular matrix");

  return mat.Triang;
}

SymMat* SimpleMatrix::sym(unsigned int, unsigned int) const
{
  if (_num != 3)
    SiconosMatrixException::selfThrow("SymMat* SimpleMatrix::sym(): the current matrix is not a Symmetric matrix");

  return mat.Sym;
}

SparseMat* SimpleMatrix::sparse(unsigned int, unsigned int) const
{
  if (_num != 4)
    SiconosMatrixException::selfThrow("SparseMat* SimpleMatrix::sparse(): the current matrix is not a Sparse matrix");

  return mat.Sparse;
}

BandedMat* SimpleMatrix::banded(unsigned int, unsigned int) const
{
  if (_num != 5)
    SiconosMatrixException::selfThrow("BandedMat* SimpleMatrix::banded(): the current matrix is not a Banded matrix");

  return mat.Banded;
}

ZeroMat* SimpleMatrix::zero_mat(unsigned int, unsigned int) const
{
  if (_num != 6)
    SiconosMatrixException::selfThrow("ZeroMat* SimpleMatrix::zero_mat(): the current matrix is not a Zero matrix");

  return mat.Zero;
}

IdentityMat* SimpleMatrix::identity(unsigned int, unsigned int) const
{
  if (_num != 7)
    SiconosMatrixException::selfThrow("IdentityMat* SimpleMatrix::identity(): the current matrix is not a Identity matrix");

  return mat.Identity;
}

double* SimpleMatrix::getArray(unsigned int, unsigned int) const
{
  if (_num == 4)
    SiconosMatrixException::selfThrow("SimpleMatrix::getArray(): not yet implemented for sparse matrix.");

  if (_num == 1)
    return (((*mat.Dense).data()).data());
  else if (_num == 2)
    return &(((*mat.Triang).data())[0]);
  else if (_num == 3)
    return &(((*mat.Sym).data())[0]);
  else if (_num == 6)
  {
    ZeroMat::iterator1 it = (*mat.Zero).begin1();
    return const_cast<double*>(&(*it));
  }
  else if (_num == 7)
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
  unsigned int size1 = size(0);
  unsigned int size2 = size(1);
  if (_num == 1)
    *mat.Dense = ublas::zero_matrix<double>(size1, size2);
  else if (_num == 2)
    *mat.Triang = ublas::zero_matrix<double>(size1, size2);

  else if (_num == 3)
    *mat.Sym = ublas::zero_matrix<double>(size1, size2);

  else if (_num == 4)
    *mat.Sparse = ublas::zero_matrix<double>(size1, size2);

  else if (_num == 5)
    *mat.Banded = ublas::zero_matrix<double>(size1, size2);

  else if (_num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix::zero(): you can not set to zero a matrix of type Identity!.");
  resetLU();
  // if _num == 6: nothing
}

void SimpleMatrix::randomize()
{
  if (_num == 1)
    Siconos::algebra::fill(*mat.Dense);
  else
    SiconosMatrixException::selfThrow("SimpleMatrix::randomize(): only implemented for dense matrices.");
  resetLU();
}

void SimpleMatrix::randomize_sym()
{
  if (_num == 1)
    Siconos::algebra::fill_sym(*mat.Dense);
  else
    SiconosMatrixException::selfThrow("SimpleMatrix::randomize_sym(): only implemented for dense matrices.");
  resetLU();
}

void SimpleMatrix::eye()
{
  unsigned int size1 = size(0);
  unsigned int size2 = size(1);
  if (_num == 1)
    *mat.Dense = ublas::identity_matrix<double>(size1, size2);

  else if (_num == 2)
    *mat.Triang = ublas::identity_matrix<double>(size1, size2);

  else if (_num == 3)
    *mat.Sym = ublas::identity_matrix<double>(size1, size2);

  else if (_num == 4)
    *mat.Sparse = ublas::identity_matrix<double>(size1, size2);

  else if (_num == 5)
    *mat.Banded = ublas::identity_matrix<double>(size1, size2);

  else if (_num == 6)
    SiconosMatrixException::selfThrow("SimpleMatrix::eye(): you can not set to identity a matrix of type Zero!.");
  resetLU();
}



unsigned int SimpleMatrix::size(unsigned int index) const
{
  if (_num == 1)
  {
    if (index == 0) return (*mat.Dense).size1();
    else  return (*mat.Dense).size2();
  }
  else if (_num == 2)
  {
   if (index == 0) return (*mat.Triang).size1();
   else return (*mat.Triang).size2();
  }
  else if (_num == 3)
  {
   if (index == 0) return (*mat.Sym).size1();
   else  return (*mat.Sym).size2();
  }
  else if (_num == 4)
  {
   if (index == 0) return (*mat.Sparse).size1();
   else return (*mat.Sparse).size2();
  }
  else if (_num == 5)
  {
  if (index == 0) return  (*mat.Banded).size1();
  else  return  (*mat.Banded).size2();
  }
  else if (_num == 6)
  {
  if (index == 0) return (*mat.Zero).size1();
  else  return (*mat.Zero).size2();
  }
  else if (_num == 7)
  {
   if (index == 0) return (*mat.Identity).size1();
   else  return (*mat.Identity).size2();
  }
  else return 0;


};


//=======================
// set matrix dimension
//=======================

void SimpleMatrix::resize(unsigned int row, unsigned int col, unsigned int lower, unsigned int upper, bool preserve)
{

  if (_num == 1)
  {
    (*mat.Dense).resize(row, col, preserve);
  }
  else if (_num == 2)
  {
    (*mat.Triang).resize(row, col, preserve);
  }
  else if (_num == 3)
  {
    (*mat.Sym).resize(row, col, preserve);
  }
  else if (_num == 4)
  {
    (*mat.Sparse).resize(row, col, preserve);
  }
  else if (_num == 5)
  {
    (*mat.Banded).resize(row, col, lower, upper, preserve);
  }
  else if (_num == 6)
  {
    (*mat.Zero).resize(row, col, preserve);
  }
  else if (_num == 7)
  {
    (*mat.Identity).resize(row, col, preserve);
  }
  resetLU();
 }


//=====================
// screen display
//=====================

void SimpleMatrix::display() const
{
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  if (_num == 1)
    Siconos::algebra::print_m(*mat.Dense);
    //std::cout << *mat.Dense << std::endl;
  else if (_num == 2)
    std::cout << *mat.Triang << std::endl;
  else if (_num == 3)
    std::cout << *mat.Sym << std::endl;
  else if (_num == 4)
    std::cout << *mat.Sparse << std::endl;
  else if (_num == 5)
    std::cout << *mat.Banded << std::endl;
  else if (_num == 6)
    std::cout << *mat.Zero << std::endl;
  else if (_num == 7)
    std::cout << *mat.Identity << std::endl;
}

//=====================
// convert to a string
//=====================

std::string SimpleMatrix::toString() const
{
  return ::toString(*this);
}

//=====================
// convert to an ostream
//=====================

std::ostream& operator<<(std::ostream& os, const SimpleMatrix& sm)
{
  if (sm._num == 1)
    os << *sm.mat.Dense;
  else if (sm._num == 2)
    os << *sm.mat.Triang;
  else if (sm._num == 3)
    os << *sm.mat.Sym;
  else if (sm._num == 4)
    os << *sm.mat.Sparse;
  else if (sm._num == 5)
    os << *sm.mat.Banded;
  else if (sm._num == 6)
    os << *sm.mat.Zero;
  else if (sm._num == 7)
    os << *sm.mat.Identity;
  return os;
}

void prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, bool init)
{

  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");


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
    private_addprod(A, startRow, startCol, **it, y);
    startCol += (*it)->size();
  }
}

void prod(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y, bool init)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int startRow = 0;
  VectorOfVectors::const_iterator it;
  // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
  //       // private_prod takes into account the fact that x and y[i] may be block vectors.
  for (it = y.begin(); it != y.end(); ++it)
  {
    private_prod(A, startRow, x, **it, init);
    startRow += (*it)->size();
  }

}

void subprod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, const Index& coord, bool init)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Number of the subvector of x that handles element at position coord[4]
  std::size_t firstBlockNum = x.getNumVectorAtPos(coord[4]);
  // Number of the subvector of x that handles element at position coord[5]
  unsigned int lastBlockNum = x.getNumVectorAtPos(coord[5]);
  Index subCoord = coord;
  SPC::SiconosVector  tmp = x[firstBlockNum];
  std::size_t subSize =  tmp->size(); // Size of the sub-vector
  const SP::Index xTab = x.tabIndex();
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
    for (VectorOfVectors::const_iterator it = x.begin(); it != x.end(); ++it)
    {
      if ((*it)->num() == 0)
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

void prod(const SiconosVector& x, const SiconosMatrix& A, BlockVector& y, bool init)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");
  unsigned int startRow = 0;
  VectorOfVectors::const_iterator it;
  // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
  // private_prod takes into account the fact that x and y[i] may be block vectors.
  for (it = y.begin(); it != y.end(); ++it)
  {
    private_prod(createSPtrConstSiconosVector(x), createSPtrConstSiconosMatrix(A), startRow, *it, init);
    startRow += (*it)->size();
  }
}

void private_addprod(const SiconosMatrix& A, unsigned startRow, unsigned int startCol, const SiconosVector& x, SiconosVector& y)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");
  assert(!A.isBlock() && "private_addprod(A,start,x,y) error: not yet implemented for block matrix.");


  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
  // Then computation of y = subA*x + y.
  unsigned int numA = A.num();
  unsigned int numY = y.num();
  unsigned int numX = x.num();
  unsigned int sizeX = x.size();
  unsigned int sizeY = y.size();

  assert(numX == numY && "private_addprod(A,start,x,y) error: not yet implemented for x and y of different types.");

  if (numY == 1 && numX == 1)
  {

    assert(y.dense() != x.dense());

    if (numA == 1)
      noalias(*y.dense()) += prod(ublas::subrange(*A.dense(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else if (numA == 2)
      noalias(*y.dense()) += prod(ublas::subrange(*A.triang(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else if (numA == 3)
      noalias(*y.dense()) += prod(ublas::subrange(*A.sym(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else if (numA == 4)
      noalias(*y.dense()) += prod(ublas::subrange(*A.sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
    else //if(numA==5)
      noalias(*y.dense()) += prod(ublas::subrange(*A.banded(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.dense());
  }
  else // x and y sparse
  {
    if (numA == 4)
      *y.sparse() += prod(ublas::subrange(*A.sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x.sparse());
    else
      SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x, y  sparse and A not sparse.");
  }
}

void private_addprod(const SiconosMatrix& A, unsigned int startRow, unsigned int startCol, const BlockVector& x, SiconosVector& y)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  assert(!A.isBlock() && "private_addprod(A,start,x,y) error: not yet implemented for block matrix.");

  VectorOfVectors::const_iterator it;
  unsigned int startColBis = startCol;
  for (it = x.begin(); it != x.end(); ++it)
  {
    private_addprod(A, startRow, startColBis, **it, y);
    startColBis += (*it)->size();
  }

}

// x block, y siconos
void private_prod(const SiconosMatrix& A, unsigned int startRow, const BlockVector& x, SiconosVector& y, bool init)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

  if (init) // y = subA * x , else y += subA * x
    y.zero();
  private_addprod(A, startRow, 0, x, y);
}

// x and y blocks
void private_prod(SPC::SiconosMatrix A, const unsigned int startRow, SPC::BlockVector x, SP::BlockVector y, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int row = startRow;
  VectorOfVectors::const_iterator it;
  for (it = y->begin(); it != y->end(); ++it)
  {
    private_prod(*A, row, *x, **it, init);
    row += (*it)->size();
  }
}

// x block, y siconos
void private_prod(const SiconosMatrix& A, unsigned int startRow, const SiconosVector& x, SiconosVector& y, bool init)
{
  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

  if (init) // y = subA * x , else y += subA * x
    y.zero();
  private_addprod(A, startRow, 0, x, y);
}

// x and y blocks
void private_prod(SPC::SiconosMatrix A, const unsigned int startRow, SPC::SiconosVector x, SP::BlockVector y, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int row = startRow;
  VectorOfVectors::const_iterator it;
  for (it = y->begin(); it != y->end(); ++it)
  {
    private_prod(*A, row, *x, **it, init);
    row += (*it)->size();
  }
}
// With trans(A) ...
void private_addprod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SP::SiconosVector y)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A->isBlock())
    SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for block matrix.");

  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
  // Then computation of y = subA*x + y.
  unsigned int numA = A->num();
  unsigned int numY = y->num();
  unsigned int numX = x->num();
  unsigned int sizeX = x->size();
  unsigned int sizeY = y->size();

  if (numX != numY)
    SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for x and y of different types.");

  if (numY == 1 && numX == 1)
  {

    assert(y->dense() != x->dense());

    if (numA == 1)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->dense()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 2)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->triang()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 3)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->sym()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 4)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->sparse()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else //if(numA==5)
      noalias(*y->dense()) += prod(ublas::subrange(trans(*A->banded()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
  }
  else // x and y sparse
  {
    if (numA == 4)
      *y->sparse() += prod(ublas::subrange(trans(*A->sparse()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->sparse());
    else
      SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for x, y  sparse and A not sparse.");
  }
}

void private_addprod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SP::SiconosVector y)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  VectorOfVectors::const_iterator it;
  unsigned int startColBis = startCol;
  for (it = x->begin(); it != x->end(); ++it)
  {
    private_addprod((*it), A, startRow, startColBis, y);
    startColBis += (*it)->size();
  }

}

void private_prod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::SiconosVector  y, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between el. of A of index (col) startCol and startCol + sizeY
  if (init) // y = subA * x , else y += subA * x
    y->zero();
  private_addprod(x, A, startCol, 0 , y);

}

void private_prod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::BlockVector  y, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int col = startCol;
  VectorOfVectors::const_iterator it;
  for (it = y->begin(); it != y->end(); ++it)
  {
    private_prod(x, A, col, *it, init);
    col += (*it)->size();
  }
}

void private_prod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::SiconosVector  y, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between el. of A of index (col) startCol and startCol + sizeY
  if (init) // y = subA * x , else y += subA * x
    y->zero();
  private_addprod(x, A, startCol, 0 , y);

}

void private_prod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::BlockVector  y, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  unsigned int col = startCol;
  VectorOfVectors::const_iterator it;
  for (it = y->begin(); it != y->end(); ++it)
  {
    private_prod(x, A, col, *it, init);
    col += (*it)->size();
  }
}

void private_addprod(double a, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SPC::SiconosVector x, SP::SiconosVector y)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A->isBlock())
    SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for block matrix.");

  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
  // Then computation of y = subA*x + y.
  unsigned int numA = A->num();
  unsigned int numY = y->num();
  unsigned int numX = x->num();
  unsigned int sizeX = x->size();
  unsigned int sizeY = y->size();

  if (numX != numY)
    SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x and y of different types.");

  if (numY == 1 && numX == 1)
  {

    assert(y->dense() != x->dense());

    if (numA == 1)
      noalias(*y->dense()) += a * prod(ublas::subrange(*A->dense(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 2)
      noalias(*y->dense()) += a * prod(ublas::subrange(*A->triang(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 3)
      noalias(*y->dense()) += a * prod(ublas::subrange(*A->sym(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 4)
      noalias(*y->dense()) += a * prod(ublas::subrange(*A->sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else //if(numA==5)
      noalias(*y->dense()) += a * prod(ublas::subrange(*A->banded(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
  }
  else // x and y sparse
  {
    if (numA == 4)
      *y->sparse() += a * prod(ublas::subrange(*A->sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->sparse());
    else
      SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x, y  sparse and A not sparse.");
  }

}

void private_prod(double a, SPC::SiconosMatrix A, unsigned int startRow, SPC::SiconosVector x, SP::SiconosVector  y, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

  if (init) // y = subA * x , else y += subA * x
    y->zero();
  private_addprod(a, A, startRow, 0, x, y);

}

unsigned SimpleMatrix::copyData(double* data) const
{
  assert((_num == 1) && "SiconosMatrix::copyData : forbidden: the current matrix is not dense.");

  unsigned size = mat.Dense->size1() * mat.Dense->size2();
  siconosBindings::detail::copy(size, getArray(), 1, data, 1);
  return size;
}
