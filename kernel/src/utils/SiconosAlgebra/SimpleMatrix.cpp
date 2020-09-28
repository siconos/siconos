/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "SimpleMatrix.hpp"

#include <assert.h>                                             // for assert
#include <memory>                                               // for __sha...
#include <iostream>                                             // for ostream
#include <boost/numeric/ublas/io.hpp>                           // for opera...
#include <boost/numeric/ublas/matrix_proxy.hpp>                 // for matri...
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include "SimpleMatrixFriends.hpp"                              // for subprod
#include "SiconosAlgebra.hpp"    // for symmetric, triangular ...
#include "BlockMatrix.hpp"                                      // for Block...
#include "BlockMatrixIterators.hpp"                             // for Const...
#include "SiconosMatrixException.hpp"                           // for Sicon...
#include "ioMatrix.hpp"                                         // for read
#include "Tools.hpp"                                            // for toString
#include "bindings_utils.hpp"                                   // for fill
#include "NumericsMatrix.h"

#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"

//#define DEBUG_MESSAGES
#include "debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#include <cs.h>
#endif





using namespace Siconos;
namespace siconosBindings = boost::numeric::bindings::blas;
using std::cout;
using std::endl;


// =================================================
//                CONSTRUCTORS
// =================================================


// Default (protected, used only for derived classes)
SimpleMatrix::SimpleMatrix(int i):
  SiconosMatrix(Siconos::DENSE),
  _isPLUFactorized(false),
  _isPLUFactorizedInPlace(false),
  _isQRFactorized(false),
  _isPLUInversed(false),
  _isCholeskyFactorized(false),
  _isCholeskyFactorizedInPlace(false)
{
  mat.Dense = new DenseMat(ublas::zero_matrix<double>());
}

SimpleMatrix::SimpleMatrix():
  SiconosMatrix(Siconos::DENSE),
  _isPLUFactorized(false),
  _isPLUFactorizedInPlace(false),
  _isQRFactorized(false),
  _isPLUInversed(false),
  _isCholeskyFactorized(false),
  _isCholeskyFactorizedInPlace(false)
{
  mat.Dense = new DenseMat(ublas::zero_matrix<double>());
}

// parameters: dimensions and type.
SimpleMatrix::SimpleMatrix(unsigned int row,
                           unsigned int col,
                           UBLAS_TYPE typ,
                           unsigned int upper,
                           unsigned int lower):
  SiconosMatrix(Siconos::DENSE),
  _isPLUFactorized(false),
  _isPLUFactorizedInPlace(false),
  _isQRFactorized(false),
  _isPLUInversed(false),
  _isCholeskyFactorized(false),
  _isCholeskyFactorizedInPlace(false)
{
  if(typ == DENSE)
  {
    mat.Dense = new DenseMat(ublas::zero_matrix<double>(row, col));
    // _num = 1; default value
  }
  else if(typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat(ublas::zero_matrix<double>(row, col));
    _num = TRIANGULAR;
  }
  else if(typ == SYMMETRIC)
  {
    mat.Sym = new SymMat(ublas::zero_matrix<double>(row, col));
    _num = SYMMETRIC;
  }
  else if(typ == SPARSE)
  {
    mat.Sparse = new SparseMat(row, col, upper);
    _num = SPARSE;
    zero();
  }
  else if(typ == SPARSE_COORDINATE)
  {
    mat.SparseCoordinate = new SparseCoordinateMat(row, col, upper);
    _num = SPARSE_COORDINATE;
    zero();
  }
  else if(typ == BANDED)
  {
    mat.Banded = new BandedMat(row, col, upper, lower);
    _num = BANDED;
    zero();
  }
  else if(typ == ZERO)
  {
    mat.Zero = new ZeroMat(row, col);
    _num = ZERO;
  }
  else if(typ == IDENTITY)
  {
    mat.Identity = new IdentityMat(row, col);
    _num = IDENTITY;
  }
  else
    SiconosMatrixException::selfThrow("SiconosMatrix::constructor(UBLAS_TYPE type, unsigned int row, unsigned int col): invalid type.");
}

// parameters: dimensions, input value and type
SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, double inputValue, UBLAS_TYPE typ, unsigned int upper, unsigned int lower):
  SiconosMatrix(typ),
  _isPLUFactorized(false),
  _isPLUFactorizedInPlace(false),
  _isQRFactorized(false),
  _isPLUInversed(false),
  _isCholeskyFactorized(false),
  _isCholeskyFactorizedInPlace(false)
{
  // This constructor has sense only for dense matrices ...
  if(typ == DENSE)
  {
    mat.Dense = new DenseMat(ublas::scalar_matrix<double>(row, col, inputValue));
    // _num = Siconos::DENSE; default value
  }
  else
    SiconosMatrixException::selfThrow("SiconosMatrix::constructor(UBLAS_TYPE type, unsigned int row, unsigned int col, double fillInValue): invalid type.");
}

// // parameters: a vector (stl) of double and the type.
// SimpleMatrix::SimpleMatrix(const std::vector<double>& v, unsigned int row, unsigned int col, UBLAS_TYPE typ, unsigned int lower, unsigned int upper):
//   SiconosMatrix(1, row, col), _isPLUFactorized(false), _isQRFactorized(false), _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
// {
//   if( (  (v.size() != row*col) && (typ != SYMMETRIC && typ != BANDED) )
//       || (v.size() != row*row && typ == SYMMETRIC)
//       || (typ == BANDED && ( (v.size()) != (unsigned int)(std::max)(row, col)*(lower+1+upper) ) ))
//     SiconosMatrixException::selfThrow("constructor(UBLAS_TYPE, const std::vector<double>, int, int) : invalid vector size");

//   if(typ == DENSE)
//     {
//       mat.Dense = new DenseMat(row,col);
//       // _num = Siconos::DENSE; default value
//     }
//   else if(typ == TRIANGULAR)
//     {
//       mat.Triang = new TriangMat(row,col);
//       _num = Siconos::TRIANGULAR;
//     }
//   else if(typ == SYMMETRIC)
//     {
//       mat.Sym = new SymMat(row);
//       _num = Siconos::SYMMETRIC;
//     }
//   else if(typ == SPARSE)
//     {
//       SiconosMatrixException::selfThrow("SimpleMatrix::constructor(UBLAS_TYPE, const std::vector<double>, int row, int col, int lower, int upper) : warning -- use constructor(const SparseMat &m) or constructor(UBLAS_TYPE, int row, int col) with UBLAS_TYPE = SPARSE");

//     }
//   else if(typ == BANDED)
//     {
//       mat.Banded = new BandedMat(row, col, lower, upper);
//       _num = Siconos::BANDED;
//     }
//   else
//     SiconosMatrixException::selfThrow("constructor(UBLAS_TYPE, const std::vector<double>, int, int) : invalid type of matrix given");

//   std::copy(v.begin(), v.end(), (vect.Dense)->begin());


// }

// Copy constructors
SimpleMatrix::SimpleMatrix(const SimpleMatrix &m):
  SiconosMatrix(m.num()),
  _isPLUFactorized(false),
  _isPLUFactorizedInPlace(false),
  _isQRFactorized(false),
  _isPLUInversed(false),
  _isCholeskyFactorized(false),
  _isCholeskyFactorizedInPlace(false)
{

  _isSymmetric = m.isSymmetric();
  _isPositiveDefinite = m.isPositiveDefinite();
  

  _isPLUFactorized= m.isPLUFactorized();
  _isPLUFactorizedInPlace= m.isPLUFactorizedInPlace();
  _isPLUInversed= m.isPLUInversed();
  

  
  if(_num == Siconos::DENSE)
  {
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    noalias(*mat.Dense) = (*m.dense());
  }
  //   mat.Dense = new DenseMat(*m.dense());

  else if(_num == Siconos::TRIANGULAR)
    mat.Triang = new TriangMat(*m.triang());

  else if(_num == Siconos::SYMMETRIC)

    mat.Sym = new SymMat(*m.sym());

  else if(_num == Siconos::SPARSE)
    mat.Sparse = new SparseMat(*m.sparse());

  else if(_num == Siconos::SPARSE_COORDINATE)
    mat.SparseCoordinate = new SparseCoordinateMat(*m.sparseCoordinate());

  else if(_num == Siconos::BANDED)
    mat.Banded = new BandedMat(*m.banded());

  else if(_num == Siconos::ZERO)
    mat.Zero = new ZeroMat(m.size(0), m.size(1));

  else// if(_num == Siconos::IDENTITY)
    mat.Identity = new IdentityMat(m.size(0), m.size(1));
}

/** copy constructor of a block given by the coord = [r0A r1A c0A c1A]
 *  \param A the matrix for extracting the block
 */
SimpleMatrix::SimpleMatrix(const SimpleMatrix& A, const Index& coord):
  SiconosMatrix(A.num()),
  _isPLUFactorized(false),
  _isPLUFactorizedInPlace(false),
  _isQRFactorized(false),
  _isPLUInversed(false),
  _isCholeskyFactorized(false),
  _isCholeskyFactorizedInPlace(false)
{
  if(coord[0]>=coord[1])
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). Empty row range coord[0]>= coord[1]");
  if(coord[2]>=coord[3])
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). Empty column range coord[2]>= coord[3]");
  if(coord[1] > A.size(0))
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). row index too large.");
  if(coord[3] > A.size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix::SimpleMatrix(const SimpleMatrix& A , const Index& coord ). column index too large.");

  if(_num == Siconos::DENSE)
  {
    ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Dense=new DenseMat(subA);
  }
  else if(_num == Siconos::TRIANGULAR)
  {
    ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Triang=new TriangMat(subA);
  }
  else if(_num == Siconos::SYMMETRIC)
  {
    ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Sym=new SymMat(subA);
  }
  else if(_num == Siconos::SPARSE)
  {
    ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Sparse=new SparseMat(subA);
  }
  else if(_num == Siconos::SPARSE_COORDINATE)
  {
    ublas::matrix_range<SparseCoordinateMat> subA(*A.sparseCoordinate(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.SparseCoordinate=new SparseCoordinateMat(subA);
  }
  else if(_num == Siconos::BANDED)
  {
    ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]), ublas::range(coord[2], coord[3]));
    mat.Banded=new BandedMat(subA);
  }
  else if(_num == Siconos::ZERO)
  {
    mat.Zero = new ZeroMat(coord[1]-coord[0], coord[3]-coord[2]);
  }
  else// if(_num == Siconos::IDENTITY)
    mat.Identity = new IdentityMat(coord[1]-coord[0], coord[3]-coord[2]);
}




SimpleMatrix::SimpleMatrix(const SiconosMatrix &m):
  SiconosMatrix(m.num()),
  _isPLUFactorized(),
  _isPLUFactorizedInPlace(false),
  _isQRFactorized(false),
  _isPLUInversed(false),
  _isCholeskyFactorized(false),
  _isCholeskyFactorizedInPlace(false)
{
  // _num is set in SiconosMatrix constructor with m.num() ... must be changed if m is Block
  unsigned int numM = m.num();

  _isSymmetric = m.isSymmetric();
  _isPositiveDefinite = m.isPositiveDefinite();

  _isPLUFactorized= m.isPLUFactorized();
  _isPLUFactorizedInPlace= m.isPLUFactorizedInPlace();
  _isPLUInversed= m.isPLUInversed();


  if(m.ipiv())
    _ipiv.reset(new VInt(*(m.ipiv())));

  if(numM == 0)  // ie if m is Block, this matrix is set to a dense.
  {
    const BlockMatrix& mB = static_cast<const BlockMatrix&>(m);
    _num = Siconos::DENSE;
    // get number of blocks in a row/col of m.
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    ConstBlocksIterator1 it;
    ConstBlocksIterator2 it2;
    unsigned int posRow = 0;
    unsigned int posCol = 0;

    for(it = mB._mat->begin1(); it != mB._mat->end1(); ++it)
    {
      for(it2 = it.begin(); it2 != it.end(); ++it2)
      {
        setBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  }
  else if(_num == Siconos::DENSE)
  {
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    noalias(*mat.Dense) = (*m.dense());
  }

  else if(_num == Siconos::TRIANGULAR)
    mat.Triang = new TriangMat(*m.triang());

  else if(_num == Siconos::SYMMETRIC)
    mat.Sym = new SymMat(*m.sym());

  else if(_num == Siconos::SPARSE)
    mat.Sparse = new SparseMat(*m.sparse());

  else if(_num == Siconos::SPARSE_COORDINATE)
    mat.SparseCoordinate = new SparseCoordinateMat(*m.sparseCoordinate());

  else if(_num == Siconos::BANDED)
    mat.Banded = new BandedMat(*m.banded());

  else if(_num == Siconos::ZERO)
    mat.Zero = new ZeroMat(m.size(0), m.size(1));

  else // if(_num == Siconos::IDENTITY)
    mat.Identity = new IdentityMat(m.size(0), m.size(1));
}

SimpleMatrix::SimpleMatrix(const DenseMat& m):
  SiconosMatrix(Siconos::DENSE),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Dense = new DenseMat(m);
}

SimpleMatrix::SimpleMatrix(const TriangMat& m):
  SiconosMatrix(Siconos::TRIANGULAR),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Triang = new TriangMat(m);
}

SimpleMatrix::SimpleMatrix(const SymMat& m):
  SiconosMatrix(Siconos::SYMMETRIC),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Sym = new SymMat(m);
}

SimpleMatrix::SimpleMatrix(const SparseMat& m):
  SiconosMatrix(Siconos::SPARSE),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Sparse = new SparseMat(m);
}

SimpleMatrix::SimpleMatrix(const SparseCoordinateMat& m):
  SiconosMatrix(SPARSE_COORDINATE),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.SparseCoordinate = new SparseCoordinateMat(m);
}

SimpleMatrix::SimpleMatrix(const BandedMat& m):
  SiconosMatrix(Siconos::BANDED),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Banded = new BandedMat(m);
}

SimpleMatrix::SimpleMatrix(const ZeroMat& m):
  SiconosMatrix(Siconos::ZERO),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Zero = new ZeroMat(m);
}

SimpleMatrix::SimpleMatrix(const IdentityMat& m):
  SiconosMatrix(Siconos::IDENTITY),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Identity = new IdentityMat(m);
}

SimpleMatrix::SimpleMatrix(const std::string &file, bool ascii):
  SiconosMatrix(Siconos::DENSE),
  _isPLUFactorized(false), _isPLUFactorizedInPlace(false), _isQRFactorized(false),
  _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
{
  mat.Dense = new DenseMat();
  if(ascii)
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
  if(_num == Siconos::DENSE)
  {
    delete(mat.Dense);
    if (_numericsMatrix)
    {
      // _numericsMatrix->matrix0 points to the array contained in the ublas matrix
      // To avoid double free on this pointer, we set it to NULL before deletion
      if (_numericsMatrix->matrix0)
        _numericsMatrix->matrix0 =nullptr;
    }
  }
  else if(_num == Siconos::TRIANGULAR)
    delete(mat.Triang);
  else if(_num == Siconos::SYMMETRIC)
    delete(mat.Sym);
  else if(_num == Siconos::SPARSE)
    delete(mat.Sparse);
  else if(_num == Siconos::BANDED)
    delete(mat.Banded);
  else if(_num == Siconos::ZERO)
    delete(mat.Zero);
  else if(_num == Siconos::IDENTITY)
    delete(mat.Identity);
}


void SimpleMatrix::updateNumericsMatrix()
{
  /* set the numericsMatrix */
  NumericsMatrix * NM;
  if(_num == DENSE)
  {
    _numericsMatrix.reset(NM_new(),NM_clear_not_dense); // When we reset, we do not free the matrix0
                                                        //that is linked to the array of the boost container
    NM = _numericsMatrix.get();
    double * data = (double*)(getArray());
    DEBUG_EXPR(NV_display(data,size(0)*size(1)););
    NM_fill(NM, NM_DENSE, size(0), size(1), data ); // Pointer link
  }
  else
  {
    // For all the other cases, we build a sparse matrix and we call numerics for the factorization of a sparse matrix.
    _numericsMatrix.reset(NM_create(NM_SPARSE, size(0), size(1)),NM_clear);
    NM = _numericsMatrix.get();
    _numericsMatrix->matrix2->origin = NSM_CSC;
    NM_csc_alloc(NM, nnz());
    fillCSC(numericsSparseMatrix(NM)->csc, std::numeric_limits<double>::epsilon());
    DEBUG_EXPR(cs_print(numericsSparseMatrix(NM)->csc, 0););
  }
}



bool SimpleMatrix::checkSymmetry(double tol) const
{
  SP::SimpleMatrix  m_trans(new SimpleMatrix(*this));
  m_trans->trans();
  double err = (*this-*m_trans).normInf();
  if((*m_trans).normInf() > 0.0)
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
  if(_num != Siconos::DENSE)
    SiconosMatrixException::selfThrow("SimpleMatrix::getDense(): the current matrix is not a Dense matrix");

  return *mat.Dense;
}

const TriangMat SimpleMatrix::getTriang(unsigned int, unsigned int) const
{
  if(_num != Siconos::TRIANGULAR)
    SiconosMatrixException::selfThrow("TriangMat SimpleMatrix::getTriang(): the current matrix is not a Triangular matrix");

  return *mat.Triang;
}

const SymMat SimpleMatrix::getSym(unsigned int, unsigned int) const
{
  if(_num != Siconos::SYMMETRIC)
    SiconosMatrixException::selfThrow("SymMat SimpleMatrix::getSym(): the current matrix is not a Symmetric matrix");

  return *mat.Sym;
}

const SparseMat SimpleMatrix::getSparse(unsigned int, unsigned int) const
{
  if(_num != Siconos::SPARSE)
    SiconosMatrixException::selfThrow("SparseMat SimpleMatrix::getSparse(): the current matrix is not a Sparse matrix");

  return *mat.Sparse;
}

const SparseCoordinateMat SimpleMatrix::getSparseCoordinate(unsigned int, unsigned int) const
{
  if(_num != Siconos::SPARSE_COORDINATE)
    SiconosMatrixException::selfThrow("SparseCoordinateMat SimpleMatrix::getSparseCoordinate(): the current matrix is not a Sparse Coordinate matrix");

  return *mat.SparseCoordinate;
}
const BandedMat SimpleMatrix::getBanded(unsigned int, unsigned int) const
{
  if(_num != Siconos::BANDED)
    SiconosMatrixException::selfThrow("BandedMat SimpleMatrix::getBanded(): the current matrix is not a Banded matrix");

  return *mat.Banded;
}

const ZeroMat SimpleMatrix::getZero(unsigned int, unsigned int) const
{
  if(_num != Siconos::ZERO)
    SiconosMatrixException::selfThrow("ZeroMat SimpleMatrix::getZero(): the current matrix is not a Zero matrix");

  return *mat.Zero;
}

const IdentityMat SimpleMatrix::getIdentity(unsigned int, unsigned int) const
{
  if(_num != Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("IdentityMat SimpleMatrix::getIdentity(): the current matrix is not a Identity matrix");

  return *mat.Identity;
}

DenseMat* SimpleMatrix::dense(unsigned int, unsigned int) const
{
  if(_num != Siconos::DENSE)
    SiconosMatrixException::selfThrow("DenseMat* SimpleMatrix::dense(): the current matrix is not a Dense matrix");

  return mat.Dense;
}

TriangMat* SimpleMatrix::triang(unsigned int, unsigned int) const
{
  if(_num != Siconos::TRIANGULAR)
    SiconosMatrixException::selfThrow("TriangMat* SimpleMatrix::triang(): the current matrix is not a Triangular matrix");

  return mat.Triang;
}

SymMat* SimpleMatrix::sym(unsigned int, unsigned int) const
{
  if(_num != Siconos::SYMMETRIC)
    SiconosMatrixException::selfThrow("SymMat* SimpleMatrix::sym(): the current matrix is not a Symmetric matrix");

  return mat.Sym;
}

SparseMat* SimpleMatrix::sparse(unsigned int, unsigned int) const
{
  if(_num != Siconos::SPARSE)
    SiconosMatrixException::selfThrow("SparseMat* SimpleMatrix::sparse(): the current matrix is not a Sparse matrix");

  return mat.Sparse;
}

SparseCoordinateMat* SimpleMatrix::sparseCoordinate(unsigned int, unsigned int) const
{
  if(_num != Siconos::SPARSE_COORDINATE)
    SiconosMatrixException::selfThrow("SparseMat* SimpleMatrix::sparse(): the current matrix is not a Sparse matrix");

  return mat.SparseCoordinate;
}

BandedMat* SimpleMatrix::banded(unsigned int, unsigned int) const
{
  if(_num != Siconos::BANDED)
    SiconosMatrixException::selfThrow("BandedMat* SimpleMatrix::banded(): the current matrix is not a Banded matrix");

  return mat.Banded;
}

ZeroMat* SimpleMatrix::zero_mat(unsigned int, unsigned int) const
{
  if(_num != Siconos::ZERO)
    SiconosMatrixException::selfThrow("ZeroMat* SimpleMatrix::zero_mat(): the current matrix is not a Zero matrix");

  return mat.Zero;
}

IdentityMat* SimpleMatrix::identity(unsigned int, unsigned int) const
{
  if(_num != Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("IdentityMat* SimpleMatrix::identity(): the current matrix is not a Identity matrix");

  return mat.Identity;
}

double* SimpleMatrix::getArray(unsigned int, unsigned int) const
{
  if(_num == Siconos::SPARSE)
    SiconosMatrixException::selfThrow("SimpleMatrix::getArray(): not yet implemented for sparse matrix.");

  if(_num == Siconos::DENSE)
    return (((*mat.Dense).data()).data());
  else if(_num == Siconos::TRIANGULAR)
    return &(((*mat.Triang).data())[0]);
  else if(_num == Siconos::SYMMETRIC)
    return &(((*mat.Sym).data())[0]);
  else if(_num == Siconos::ZERO)
  {
    ZeroMat::iterator1 it = (*mat.Zero).begin1();
    return const_cast<double*>(&(*it));
  }
  else if(_num == Siconos::IDENTITY)
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
  if(_num == Siconos::DENSE)
    *mat.Dense = ublas::zero_matrix<double>(size1, size2);
  else if(_num == Siconos::TRIANGULAR)
    *mat.Triang = ublas::zero_matrix<double>(size1, size2);

  else if(_num == Siconos::SYMMETRIC)
    *mat.Sym = ublas::zero_matrix<double>(size1, size2);

  else if(_num == Siconos::SPARSE)
    *mat.Sparse = ublas::zero_matrix<double>(size1, size2);

  else if(_num == Siconos::SPARSE_COORDINATE)
    *mat.SparseCoordinate = ublas::zero_matrix<double>(size1, size2);

  else if(_num == Siconos::BANDED)
    *mat.Banded = ublas::zero_matrix<double>(size1, size2);

  else if(_num == Siconos::IDENTITY)
    SiconosMatrixException::selfThrow("SimpleMatrix::zero(): you can not set to zero a matrix of type Identity!.");
  resetFactorizationFlags();
  // if _num == Siconos::ZERO: nothing
}

void SimpleMatrix::randomize()
{
  if(_num == Siconos::DENSE)
    Siconos::algebra::fill(*mat.Dense);
  else
    SiconosMatrixException::selfThrow("SimpleMatrix::randomize(): only implemented for dense matrices.");
  resetFactorizationFlags();
}

void SimpleMatrix::randomize_sym()
{
  if(_num == Siconos::DENSE)
    Siconos::algebra::fill_sym(*mat.Dense);
  else
    SiconosMatrixException::selfThrow("SimpleMatrix::randomize_sym(): only implemented for dense matrices.");
  resetFactorizationFlags();
}

void SimpleMatrix::eye()
{
  unsigned int size1 = size(0);
  unsigned int size2 = size(1);
  if(_num == Siconos::DENSE)
    *mat.Dense = ublas::identity_matrix<double>(size1, size2);

  else if(_num == Siconos::TRIANGULAR)
    *mat.Triang = ublas::identity_matrix<double>(size1, size2);

  else if(_num == Siconos::SYMMETRIC)
    *mat.Sym = ublas::identity_matrix<double>(size1, size2);

  else if(_num == Siconos::SPARSE)
    *mat.Sparse = ublas::identity_matrix<double>(size1, size2);

  else if(_num == Siconos::BANDED)
    *mat.Banded = ublas::identity_matrix<double>(size1, size2);

  else if(_num == Siconos::ZERO)
    SiconosMatrixException::selfThrow("SimpleMatrix::eye(): you can not set to identity a matrix of type Zero!.");
  resetFactorizationFlags();
}



unsigned int SimpleMatrix::size(unsigned int index) const
{
  if(_num == Siconos::DENSE)
  {
    if(index == 0) return (*mat.Dense).size1();
    else  return (*mat.Dense).size2();
  }
  else if(_num == Siconos::TRIANGULAR)
  {
    if(index == 0) return (*mat.Triang).size1();
    else return (*mat.Triang).size2();
  }
  else if(_num == Siconos::SYMMETRIC)
  {
    if(index == 0) return (*mat.Sym).size1();
    else  return (*mat.Sym).size2();
  }
  else if(_num == Siconos::SPARSE)
  {
    if(index == 0) return (*mat.Sparse).size1();
    else return (*mat.Sparse).size2();
  }
  else if(_num == Siconos::SPARSE_COORDINATE)
  {
    if(index == 0) return (*mat.SparseCoordinate).size1();
    else return (*mat.SparseCoordinate).size2();
  }
  else if(_num == Siconos::BANDED)
  {
    if(index == 0) return (*mat.Banded).size1();
    else  return (*mat.Banded).size2();
  }
  else if(_num == Siconos::ZERO)
  {
    if(index == 0) return (*mat.Zero).size1();
    else  return (*mat.Zero).size2();
  }
  else if(_num == Siconos::IDENTITY)
  {
    if(index == 0) return (*mat.Identity).size1();
    else  return (*mat.Identity).size2();
  }
  else return 0;


};


//=======================
// set matrix dimension
//=======================

void SimpleMatrix::resize(unsigned int row, unsigned int col, unsigned int lower, unsigned int upper, bool preserve)
{

  if(_num == Siconos::DENSE)
  {
    (*mat.Dense).resize(row, col, preserve);
  }
  else if(_num == Siconos::TRIANGULAR)
  {
    (*mat.Triang).resize(row, col, preserve);
  }
  else if(_num == Siconos::SYMMETRIC)
  {
    (*mat.Sym).resize(row, col, preserve);
  }
  else if(_num == Siconos::SPARSE)
  {
    (*mat.Sparse).resize(row, col, preserve);
  }
  else if(_num == Siconos::SPARSE_COORDINATE)
  {
    (*mat.SparseCoordinate).resize(row, col, preserve);
  }
  else if(_num == Siconos::BANDED)
  {
    (*mat.Banded).resize(row, col, lower, upper, preserve);
  }
  else if(_num == Siconos::ZERO)
  {
    (*mat.Zero).resize(row, col, preserve);
  }
  else if(_num == Siconos::IDENTITY)
  {
    (*mat.Identity).resize(row, col, preserve);
  }
  resetFactorizationFlags();
}


//=====================
// screen display
//=====================

void SimpleMatrix::display() const
{
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  if(size(0) == 0 || size(1) ==0)
  {
    std::cout << "SimpleMatrix::display(): empty matrix" << std::endl;
  }
  std::cout << "SimpleMatrix storage type - num = " << _num << "\n";
  if(_num == Siconos::DENSE)
  {
    Siconos::algebra::print_m(*mat.Dense);
    //std::cout << *mat.Dense << std::endl;
  }
  else if(_num == Siconos::TRIANGULAR)
    std::cout << *mat.Triang << std::endl;
  else if(_num == Siconos::SYMMETRIC)
    std::cout << *mat.Sym << std::endl;
  else if(_num == Siconos::SPARSE)
  {
    std::cout << "non zero element (nnz) = " <<  mat.Sparse->nnz() << std::endl;

    std::cout << *mat.Sparse << std::endl;
  }
  else if(_num == Siconos::SPARSE_COORDINATE)
  {
    std::cout << *mat.SparseCoordinate << std::endl;
  }
  else if(_num == Siconos::BANDED)
    std::cout << *mat.Banded << std::endl;
  else if(_num == Siconos::ZERO)
    std::cout << *mat.Zero << std::endl;
  else if(_num == Siconos::IDENTITY)
    std::cout << *mat.Identity << std::endl;
}
void SimpleMatrix::displayExpert(bool brief) const
{
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  if(size(0) == 0 || size(1) ==0)
  {
    std::cout << "SimpleMatrix::display(): empty matrix" << std::endl;
  }
  std::cout << "SimpleMatrix storage type - num = " << _num << "\n";
  if(_num == Siconos::DENSE)
  {
    Siconos::algebra::print_m(*mat.Dense);
    //std::cout << *mat.Dense << std::endl;
  }
  else if(_num == Siconos::TRIANGULAR)
    std::cout << *mat.Triang << std::endl;
  else if(_num == Siconos::SYMMETRIC)
    std::cout << *mat.Sym << std::endl;
  else if(_num == Siconos::SPARSE)
  {
    std::cout << "non zero element (nnz) = " <<  mat.Sparse->nnz() << std::endl;
    std::cout << "non zero element (nnz_capacity) = " <<  mat.Sparse->nnz_capacity() << std::endl;
    std::cout << "filled1 = " <<  mat.Sparse->filled1() << std::endl;
    std::cout << "filled2 = " <<  mat.Sparse->filled2() << std::endl;

    std::cout << "index_data1 = [ " ;
    size_t i=0;
    for(i = 0; i < mat.Sparse->filled1()-1 ; i++)
    {
      std::cout << mat.Sparse->index1_data()[i] << ", " ;
    }
    std::cout << mat.Sparse->index1_data()[i] << "]" <<  std::endl;

    std::cout << "index_data2 = [" ;
    for(i = 0; i < mat.Sparse->filled2()-1 ; i++)
    {
      std::cout << mat.Sparse->index2_data()[i] << ", " ;
    }
    std::cout << mat.Sparse->index2_data()[i] << "]" << std::endl;

    std::cout << "value_data = [" ;
    for(i = 0; i < mat.Sparse->filled2()-1 ; i++)
    {
      std::cout << mat.Sparse->value_data()[i] << ", " ;
    }
    std::cout << mat.Sparse->value_data()[i] << "]" << std::endl;

    std::cout << *mat.Sparse << std::endl;
  }
  else if(_num == Siconos::SPARSE_COORDINATE)
  {
    std::cout << "non zero element (nnz) = " <<  mat.SparseCoordinate->nnz() << std::endl;



    for(size_t i = 0; i < mat.SparseCoordinate->nnz(); ++i)
    {
      //std::cout << i << std::endl;
      std::cout << "M(" << mat.SparseCoordinate->index1_data()[i] << ", " ;
      std::cout << mat.SparseCoordinate->index2_data()[i] << ") =  " ;
      std::cout << mat.SparseCoordinate->value_data()[i] << std::endl;
    }
  }
  else if(_num == Siconos::BANDED)
    std::cout << *mat.Banded << std::endl;
  else if(_num == Siconos::ZERO)
    std::cout << *mat.Zero << std::endl;
  else if(_num == Siconos::IDENTITY)
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
  if(sm._num == Siconos::DENSE)
    os << *sm.mat.Dense;
  else if(sm._num == Siconos::TRIANGULAR)
    os << *sm.mat.Triang;
  else if(sm._num == Siconos::SYMMETRIC)
    os << *sm.mat.Sym;
  else if(sm._num == Siconos::SPARSE)
    os << *sm.mat.Sparse;
  else if(sm._num == Siconos::BANDED)
    os << *sm.mat.Banded;
  else if(sm._num == Siconos::ZERO)
    os << *sm.mat.Zero;
  else if(sm._num == Siconos::IDENTITY)
    os << *sm.mat.Identity;
  return os;
}


void SimpleMatrix::assign(const SimpleMatrix &smat)
{

  switch(_num)
  {
  case Siconos::SPARSE:
  {


    switch(smat.num())
    {
    case Siconos::SPARSE:
    {
      mat.Sparse->assign(smat.getSparse());
      break;
    }
    default:
    {
    }

    }
  }
  default:
  {
    SiconosMatrixException::selfThrow("SimpleMatrix::assign(const SimpleMatrix& A) : do not know how to assign for the given storage type ");
  }
  }
}

// void prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, bool init)
// {
//   assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
//   if(init)
//     y.zero();
//   unsigned int startRow = 0;
//   unsigned int startCol = 0;
//   // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA x[i]), with subA a submatrix of A,
//   // starting from position startRow in rows and startCol in columns.
//   // private_prod takes also into account the fact that each block of x can also be a block.
//   VectorOfVectors::const_iterator it;
//   for(it = x.begin(); it != x.end(); ++it)
//   {
//     private_addprod(A, startRow, startCol, **it, y);
//     startCol += (*it)->size();
//   }
// }



// void private_addprod(const SiconosMatrix& A, unsigned int startRow, unsigned int startCol, const BlockVector& x, SiconosVector& y)
// {
//   assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
//   assert(!A.isBlock() && "private_addprod(A,start,x,y) error: not yet implemented for block matrix.");
//   VectorOfVectors::const_iterator it;
//   unsigned int startColBis = startCol;
//   for(it = x.begin(); it != x.end(); ++it)
//   {
//     private_addprod(A, startRow, startColBis, **it, y);
//     startColBis += (*it)->size();
//   }

// }

// // x block, y siconos
// void private_prod(const SiconosMatrix& A, unsigned int startRow, const BlockVector& x, SiconosVector& y, bool init)
// {
//   assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY
//   if(init)  // y = subA * x , else y += subA * x
//     y.zero();
//   private_addprod(A, startRow, 0, x, y);
// }

// // x and y blocks
// void private_prod(SPC::SiconosMatrix A, const unsigned int startRow, SPC::BlockVector x, SP::BlockVector y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   unsigned int row = startRow;
//   VectorOfVectors::const_iterator it;
//   for(it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(*A, row, *x, **it, init);
//     row += (*it)->size();
//   }
// }

// // x and y blocks
// void private_prod(SPC::SiconosMatrix A, const unsigned int startRow, SPC::SiconosVector x, SP::BlockVector y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   unsigned int row = startRow;
//   VectorOfVectors::const_iterator it;
//   for(it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(*A, row, *x, **it, init);
//     row += (*it)->size();
//   }
// }

// void private_addprod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SP::SiconosVector y)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
//   VectorOfVectors::const_iterator it;
//   unsigned int startColBis = startCol;
//   for(it = x->begin(); it != x->end(); ++it)
//   {
//     private_addprod((*it), A, startRow, startColBis, y);
//     startColBis += (*it)->size();
//   }

// }

// void private_prod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::BlockVector  y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   unsigned int col = startCol;
//   VectorOfVectors::const_iterator it;
//   for(it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(x, A, col, *it, init);
//     col += (*it)->size();
//   }
// }

// void private_prod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::SiconosVector  y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between el. of A of index (col) startCol and startCol + sizeY
//   if(init)  // y = subA * x , else y += subA * x
//     y->zero();
//   private_addprod(x, A, startCol, 0, y);

// }

// void private_prod(SPC::BlockVector x, SPC::SiconosMatrix A, unsigned int startCol, SP::BlockVector  y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   unsigned int col = startCol;
//   VectorOfVectors::const_iterator it;
//   for(it = y->begin(); it != y->end(); ++it)
//   {
//     private_prod(x, A, col, *it, init);
//     col += (*it)->size();
//   }
// }

// void private_addprod(double a, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SPC::SiconosVector x, SP::SiconosVector y)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   if(A->isBlock())
//     SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for block matrix.");

//   // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
//   // Then computation of y = subA*x + y.
//   unsigned int numA = A->num();
//   unsigned int numY = y->num();
//   unsigned int numX = x->num();
//   unsigned int sizeX = x->size();
//   unsigned int sizeY = y->size();

//   if(numX != numY)
//     SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x and y of different types.");

//   if(numY == 1 && numX == 1)
//   {

//     assert(y->dense() != x->dense());

//     if(numA == 1)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->dense(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else if(numA == 2)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->triang(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else if(numA == 3)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->sym(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else if(numA == 4)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//     else //if(numA==5)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->banded(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
//   }
//   else // x and y sparse
//   {
//     if(numA == 4)
//       *y->sparse() += a * prod(ublas::subrange(*A->sparse(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->sparse());
//     else
//       SiconosMatrixException::selfThrow("private_addprod(A,start,x,y) error: not yet implemented for x, y  sparse and A not sparse.");
//   }

// }

// void private_prod(double a, SPC::SiconosMatrix A, unsigned int startRow, SPC::SiconosVector x, SP::SiconosVector  y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");


//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el. of index (row) startRow and startRow + sizeY

//   if(init)  // y = subA * x , else y += subA * x
//     y->zero();
//   private_addprod(a, A, startRow, 0, x, y);

// }

unsigned SimpleMatrix::copyData(double* data) const
{
  assert((_num == Siconos::DENSE) && "SiconosMatrix::copyData : forbidden: the current matrix is not dense.");

  unsigned size = mat.Dense->size1() * mat.Dense->size2();
  siconosBindings::detail::copy(size, getArray(), 1, data, 1);
  return size;
}
