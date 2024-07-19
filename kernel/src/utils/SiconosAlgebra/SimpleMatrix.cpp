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
#include "SimpleMatrix.hpp"

// #include <assert.h>  // for assert

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>            // for opera...
#include <boost/numeric/ublas/matrix_proxy.hpp>  // for matri...
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "BlockMatrix.hpp"           // for Block...
#include "NumericsToolsNamespace.h"  // for NumericsMatrix and NumericsSparseMatrix
#include "SiconosAlgebraTools.hpp"   // for randomize ...
#include "SiconosAlgebraTypes.hpp"   // for UblasType
#include "SiconosException.hpp"      // for Sicon...
#include "SiconosMatrixOp.hpp"       // matrix operators declarations
#include "bindings_utils.hpp"        // for fill
#include "boost/numeric/bindings/blas.hpp"
#include "boost/numeric/bindings/ublas/matrix.hpp"
#include "io.hpp"  // for read
// #define DEBUG_MESSAGES
#include "siconos_debug.h"
#ifdef DEBUG_MESSAGES
#include <cs.h>
#endif

namespace ublas = boost::numeric::ublas;
namespace bindings_blas = boost::numeric::bindings::blas;

// =================================================
//                CONSTRUCTORS
// =================================================

// siconos::algebra::SimpleMatrix::SimpleMatrix() : SiconosMatrix(UblasType::DENSE)
// {
//   mat.Dense = new DenseMat(ublas::zero_matrix<double>());
// }

// parameters: dimensions and type.
siconos::algebra::SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, UblasType typ,
                                             unsigned int upper, unsigned int lower)
    : SiconosMatrix(UblasType::DENSE) {
  if (typ == UblasType::DENSE) {
    mat.Dense = new DenseMat(ublas::zero_matrix<double>(row, col));
    // _num = 1; default value
  } else if (typ == UblasType::TRIANGULAR) {
    mat.Triang = new TriangMat(ublas::zero_matrix<double>(row, col));
    _num = UblasType::TRIANGULAR;
  } else if (typ == UblasType::SYMMETRIC) {
    mat.Sym = new SymMat(ublas::zero_matrix<double>(row, col));
    _num = UblasType::SYMMETRIC;
  } else if (typ == UblasType::SPARSE) {
    mat.Sparse = new SparseMat(row, col, upper);
    _num = UblasType::SPARSE;
    zero();
  } else if (typ == UblasType::SPARSE_COORDINATE) {
    mat.SparseCoordinate = new SparseCoordinateMat(row, col, upper);
    _num = UblasType::SPARSE_COORDINATE;
    zero();
  } else if (typ == UblasType::BANDED) {
    mat.Banded = new BandedMat(row, col, upper, lower);
    _num = UblasType::BANDED;
    zero();
  } else if (typ == UblasType::ZERO) {
    mat.Zero = new ZeroMat(row, col);
    _num = UblasType::ZERO;
  } else if (typ == UblasType::IDENTITY) {
    mat.Identity = new IdentityMat(row, col);
    _num = UblasType::IDENTITY;
  } else
    THROW_EXCEPTION("invalid type.");
}

// parameters: dimensions, input value and type
siconos::algebra::SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col,
                                             double inputValue, UblasType typ,
                                             unsigned int upper, unsigned int lower)
    : SiconosMatrix(typ) {
  // This constructor has sense only for dense matrices ...
  if (typ == UblasType::DENSE) {
    mat.Dense = new DenseMat(ublas::scalar_matrix<double>(row, col, inputValue));
    // _num = UblasType::DENSE; default value
  } else
    THROW_EXCEPTION("invalid type.");
}

// // parameters: a vector (stl) of double and the type.
// siconos::algebra::SimpleMatrix::SimpleMatrix(const std::vector<double>& v, unsigned int row,
// unsigned int col, UblasType typ, unsigned int lower, unsigned int upper):
//   SiconosMatrix(1, row, col), _isPLUFactorized(false), _isQRFactorized(false),
//   _isPLUInversed(false), _isCholeskyFactorized(false), _isCholeskyFactorizedInPlace(false)
// {
//   if( (  (v.size() != row*col) && (typ != UblasType::SYMMETRIC && typ !=
//   UblasType::BANDED) )
//       || (v.size() != row*row && typ == UblasType::SYMMETRIC)
//       || (typ == UblasType::BANDED && ( (v.size()) != (unsigned int)(std::max)(row,
//       col)*(lower+1+upper) ) ))
//      THROW_EXCEPTION("invalid vector size");

//   if(typ == UblasType::DENSE)
//     {
//       mat.Dense = new DenseMat(row,col);
//       // _num = UblasType::DENSE; default value
//     }
//   else if(typ == UblasType::TRIANGULAR)
//     {
//       mat.Triang = new TriangMat(row,col);
//       _num = UblasType::TRIANGULAR;
//     }
//   else if(typ == UblasType::SYMMETRIC)
//     {
//       mat.Sym = new SymMat(row);
//       _num = UblasType::SYMMETRIC;
//     }
//   else if(typ == UblasType::SPARSE)
//     {
//        THROW_EXCEPTION("warning -- use constructor(const SparseMat &m) or
//        constructor(UblasType, int row, int col) with UblasType = UblasType::SPARSE");

//     }
//   else if(typ == UblasType::BANDED)
//     {
//       mat.Banded = new BandedMat(row, col, lower, upper);
//       _num = UblasType::UblasType::BANDED;
//     }
//   else
//      THROW_EXCEPTION("invalid type of matrix given");

//   std::copy(v.begin(), v.end(), (vect.Dense)->begin());

// }

// Copy constructors
siconos::algebra::SimpleMatrix::SimpleMatrix(const SimpleMatrix &m) : SiconosMatrix(m.num()) {
  _isSymmetric = m.isSymmetric();
  _isPositiveDefinite = m.isPositiveDefinite();

  _isPLUFactorized = m.isPLUFactorized();
  _isPLUFactorizedInPlace = m.isPLUFactorizedInPlace();
  _isPLUInversed = m.isPLUInversed();

  if (_num == UblasType::DENSE) {
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    noalias(*mat.Dense) = (*m.dense());
  }
  //   mat.Dense = new DenseMat(*m.dense());

  else if (_num == UblasType::TRIANGULAR)
    mat.Triang = new TriangMat(*m.triang());

  else if (_num == UblasType::SYMMETRIC)

    mat.Sym = new SymMat(*m.sym());

  else if (_num == UblasType::SPARSE)
    mat.Sparse = new SparseMat(*m.sparse());

  else if (_num == UblasType::SPARSE_COORDINATE)
    mat.SparseCoordinate = new SparseCoordinateMat(*m.sparseCoordinate());

  else if (_num == UblasType::BANDED)
    mat.Banded = new BandedMat(*m.banded());

  else if (_num == UblasType::ZERO)
    mat.Zero = new ZeroMat(m.size(0), m.size(1));

  else  // if(_num == UblasType::IDENTITY)
    mat.Identity = new IdentityMat(m.size(0), m.size(1));
}

/** copy constructor of a block given by the coord = [r0A r1A c0A c1A]
 *  \param A the matrix for extracting the block
 */
siconos::algebra::SimpleMatrix::SimpleMatrix(const SimpleMatrix &A,
                                             const std::vector<std::size_t> &coord)
    : SiconosMatrix(A.num()) {
  if (coord[0] >= coord[1]) THROW_EXCEPTION("Empty row range coord[0]>= coord[1]");
  if (coord[2] >= coord[3]) THROW_EXCEPTION("Empty column range coord[2]>= coord[3]");
  if (coord[1] > A.size(0)) THROW_EXCEPTION("row index too large.");
  if (coord[3] > A.size(1)) THROW_EXCEPTION("column index too large.");

  if (_num == UblasType::DENSE) {
    ublas::matrix_range<DenseMat> subA(*A.dense(), ublas::range(coord[0], coord[1]),
                                       ublas::range(coord[2], coord[3]));
    mat.Dense = new DenseMat(subA);
  } else if (_num == UblasType::TRIANGULAR) {
    ublas::matrix_range<TriangMat> subA(*A.triang(), ublas::range(coord[0], coord[1]),
                                        ublas::range(coord[2], coord[3]));
    mat.Triang = new TriangMat(subA);
  } else if (_num == UblasType::SYMMETRIC) {
    ublas::matrix_range<SymMat> subA(*A.sym(), ublas::range(coord[0], coord[1]),
                                     ublas::range(coord[2], coord[3]));
    mat.Sym = new SymMat(subA);
  } else if (_num == UblasType::SPARSE) {
    ublas::matrix_range<SparseMat> subA(*A.sparse(), ublas::range(coord[0], coord[1]),
                                        ublas::range(coord[2], coord[3]));
    mat.Sparse = new SparseMat(subA);
  } else if (_num == UblasType::SPARSE_COORDINATE) {
    ublas::matrix_range<SparseCoordinateMat> subA(*A.sparseCoordinate(),
                                                  ublas::range(coord[0], coord[1]),
                                                  ublas::range(coord[2], coord[3]));
    mat.SparseCoordinate = new SparseCoordinateMat(subA);
  } else if (_num == UblasType::BANDED) {
    ublas::matrix_range<BandedMat> subA(*A.banded(), ublas::range(coord[0], coord[1]),
                                        ublas::range(coord[2], coord[3]));
    mat.Banded = new BandedMat(subA);
  } else if (_num == UblasType::ZERO) {
    mat.Zero = new ZeroMat(coord[1] - coord[0], coord[3] - coord[2]);
  } else  // if(_num == UblasType::IDENTITY)
    mat.Identity = new IdentityMat(coord[1] - coord[0], coord[3] - coord[2]);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const SiconosMatrix &m) : SiconosMatrix(m.num()) {
  // _num is set in SiconosMatrix constructor with m.num() ... must be changed if m is Block
  auto numM = m.num();

  _isSymmetric = m.isSymmetric();
  _isPositiveDefinite = m.isPositiveDefinite();

  _isPLUFactorized = m.isPLUFactorized();
  _isPLUFactorizedInPlace = m.isPLUFactorizedInPlace();
  _isPLUInversed = m.isPLUInversed();

  if (numM != UblasType::BLOCK) {
    const SimpleMatrix &mm = static_cast<const SimpleMatrix &>(m);
    if (mm.ipiv()) {
      _ipiv = std::make_shared<VInt>(*(mm.ipiv()));
    }
  }

  if (numM == UblasType::BLOCK)  // ie if m is Block, this matrix is set to a dense.
  {
    const BlockMatrix &mB = static_cast<const BlockMatrix &>(m);
    _num = UblasType::DENSE;
    // get number of blocks in a row/col of m.
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    unsigned int posRow = 0;
    unsigned int posCol = 0;

    for (auto it = mB._mat->begin1(); it != mB._mat->end1(); ++it) {
      for (auto it2 = it.begin(); it2 != it.end(); ++it2) {
        setBlock(posRow, posCol, **it2);
        posCol += (*it2)->size(1);
      }
      posRow += (*it)->size(0);
      posCol = 0;
    }
  } else if (_num == UblasType::DENSE) {
    mat.Dense = new DenseMat(m.size(0), m.size(1));
    noalias(*mat.Dense) = (*m.dense());
  }

  else if (_num == UblasType::TRIANGULAR)
    mat.Triang = new TriangMat(*m.triang());

  else if (_num == UblasType::SYMMETRIC)
    mat.Sym = new SymMat(*m.sym());

  else if (_num == UblasType::SPARSE)
    mat.Sparse = new SparseMat(*m.sparse());

  else if (_num == UblasType::SPARSE_COORDINATE)
    mat.SparseCoordinate = new SparseCoordinateMat(*m.sparseCoordinate());

  else if (_num == UblasType::BANDED)
    mat.Banded = new BandedMat(*m.banded());

  else if (_num == UblasType::ZERO)
    mat.Zero = new ZeroMat(m.size(0), m.size(1));

  else  // if(_num == UblasType::IDENTITY)
    mat.Identity = new IdentityMat(m.size(0), m.size(1));
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const DenseMat &m)
    : SiconosMatrix(UblasType::DENSE) {
  mat.Dense = new DenseMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const TriangMat &m)
    : SiconosMatrix(UblasType::TRIANGULAR) {
  mat.Triang = new TriangMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const SymMat &m)
    : SiconosMatrix(UblasType::SYMMETRIC) {
  mat.Sym = new SymMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const SparseMat &m)
    : SiconosMatrix(UblasType::SPARSE) {
  mat.Sparse = new SparseMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const SparseCoordinateMat &m)
    : SiconosMatrix(UblasType::SPARSE_COORDINATE) {
  mat.SparseCoordinate = new SparseCoordinateMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const BandedMat &m)
    : SiconosMatrix(UblasType::BANDED) {
  mat.Banded = new BandedMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const ZeroMat &m)
    : SiconosMatrix(UblasType::ZERO) {
  mat.Zero = new ZeroMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const IdentityMat &m)
    : SiconosMatrix(UblasType::IDENTITY) {
  mat.Identity = new IdentityMat(m);
}

siconos::algebra::SimpleMatrix::SimpleMatrix(const std::string &file, bool ascii)
    : SiconosMatrix(UblasType::DENSE) {
  mat.Dense = new DenseMat();
  if (ascii) {
    io::read(file, *this, io::ASCII_IN);
  } else {
    io::read(file, *this, io::BINARY_IN);
  }
}

siconos::algebra::SimpleMatrix::~SimpleMatrix() noexcept {
  if (_num == UblasType::DENSE) {
    delete (mat.Dense);
    if (_numericsMatrix) {
      // _numericsMatrix->matrix0 points to the array contained in the ublas matrix
      // To avoid double free on this pointer, we set it to NULL before deletion
      if (_numericsMatrix->matrix0) _numericsMatrix->matrix0 = nullptr;
    }
  } else if (_num == UblasType::TRIANGULAR)
    delete (mat.Triang);
  else if (_num == UblasType::SYMMETRIC)
    delete (mat.Sym);
  else if (_num == UblasType::SPARSE_COORDINATE)
    delete (mat.SparseCoordinate);
  else if (_num == UblasType::SPARSE)
    delete (mat.Sparse);
  else if (_num == UblasType::BANDED)
    delete (mat.Banded);
  else if (_num == UblasType::ZERO)
    delete (mat.Zero);
  else if (_num == UblasType::IDENTITY)
    delete (mat.Identity);
}

void siconos::algebra::SimpleMatrix::updateNumericsMatrix() {
  /* set the numericsMatrix */
  NumericsMatrix *NM;
  if (_num == UblasType::DENSE) {
    _numericsMatrix.reset(NM_new(),
                          NM_free_not_dense);  // When we reset, we do not free the matrix0
                                               // that is linked to the array of the boost
                                               // container
    NM = _numericsMatrix.get();
    double *data = (double *)(getArray());
    DEBUG_EXPR(NV_display(data, size(0) * size(1)););
    NM_fill(NM, NM_DENSE, size(0), size(1),
            data);  // Pointer link
  } else {
    // For all the other cases, we build a sparse matrix and we call numerics for the
    // factorization of a sparse matrix.
    _numericsMatrix.reset(NM_create(NM_SPARSE, size(0), size(1)), NM_free);
    NM = _numericsMatrix.get();
    _numericsMatrix->matrix2->origin = NSM_CSC;
    NM_csc_alloc(NM, nnz());
    fillCSC(numericsSparseMatrix(NM)->csc, std::numeric_limits<double>::epsilon());
    DEBUG_EXPR(cs_print(numericsSparseMatrix(NM)->csc, 0););
  }
}

bool siconos::algebra::SimpleMatrix::checkSymmetry(double tol) const {
  auto m_trans = std::make_shared<SimpleMatrix>(*this);
  m_trans->trans();
  double err = (*this - *m_trans).normInf();
  if ((*m_trans).normInf() > 0.0) {
    err /= (*m_trans).normInf();
  }
  // std::cout << "err_rel  ="<< err <<"\n";
  return (err < tol);
}
//======================================
// get Ublas component (dense, sym ...)
//======================================

const siconos::algebra::DenseMat siconos::algebra::SimpleMatrix::getDense(unsigned int,
                                                                          unsigned int) const {
  if (_num != UblasType::DENSE) THROW_EXCEPTION(" the current matrix is not a Dense matrix");

  return *mat.Dense;
}

const siconos::algebra::TriangMat siconos::algebra::SimpleMatrix::getTriang(
    unsigned int, unsigned int) const {
  if (_num != UblasType::TRIANGULAR)
    THROW_EXCEPTION("the current matrix is not a Triangular matrix");

  return *mat.Triang;
}

const siconos::algebra::SymMat siconos::algebra::SimpleMatrix::getSym(unsigned int,
                                                                      unsigned int) const {
  if (_num != UblasType::SYMMETRIC)
    THROW_EXCEPTION("he current matrix is not a Symmetric matrix");

  return *mat.Sym;
}

const siconos::algebra::SparseMat siconos::algebra::SimpleMatrix::getSparse(
    unsigned int, unsigned int) const {
  if (_num != UblasType::SPARSE) THROW_EXCEPTION("the current matrix is not a Sparse matrix");

  return *mat.Sparse;
}

const siconos::algebra::SparseCoordinateMat
siconos::algebra::SimpleMatrix::getSparseCoordinate(unsigned int, unsigned int) const {
  if (_num != UblasType::SPARSE_COORDINATE)
    THROW_EXCEPTION("the current matrix is not a Sparse Coordinate matrix");

  return *mat.SparseCoordinate;
}
const siconos::algebra::BandedMat siconos::algebra::SimpleMatrix::getBanded(
    unsigned int, unsigned int) const {
  if (_num != UblasType::BANDED) THROW_EXCEPTION("the current matrix is not a Banded matrix");

  return *mat.Banded;
}

const siconos::algebra::ZeroMat siconos::algebra::SimpleMatrix::getZero(unsigned int,
                                                                        unsigned int) const {
  if (_num != UblasType::ZERO) THROW_EXCEPTION("the current matrix is not a Zero matrix");

  return *mat.Zero;
}

const siconos::algebra::IdentityMat siconos::algebra::SimpleMatrix::getIdentity(
    unsigned int, unsigned int) const {
  if (_num != UblasType::IDENTITY)
    THROW_EXCEPTION("the current matrix is not a Identity matrix");

  return *mat.Identity;
}

siconos::algebra::DenseMat *siconos::algebra::SimpleMatrix::dense(unsigned int,
                                                                  unsigned int) const {
  if (_num != UblasType::DENSE) THROW_EXCEPTION("the current matrix is not a Dense matrix");

  return mat.Dense;
}

siconos::algebra::TriangMat *siconos::algebra::SimpleMatrix::triang(unsigned int,
                                                                    unsigned int) const {
  if (_num != UblasType::TRIANGULAR)
    THROW_EXCEPTION("the current matrix is not a Triangular matrix");

  return mat.Triang;
}

siconos::algebra::SymMat *siconos::algebra::SimpleMatrix::sym(unsigned int,
                                                              unsigned int) const {
  if (_num != UblasType::SYMMETRIC)
    THROW_EXCEPTION("the current matrix is not a Symmetric matrix");

  return mat.Sym;
}

siconos::algebra::SparseMat *siconos::algebra::SimpleMatrix::sparse(unsigned int,
                                                                    unsigned int) const {
  if (_num != UblasType::SPARSE) THROW_EXCEPTION("the current matrix is not a Sparse matrix");

  return mat.Sparse;
}

siconos::algebra::SparseCoordinateMat *siconos::algebra::SimpleMatrix::sparseCoordinate(
    unsigned int, unsigned int) const {
  if (_num != UblasType::SPARSE_COORDINATE)
    THROW_EXCEPTION("the current matrix is not a Sparse matrix");

  return mat.SparseCoordinate;
}

siconos::algebra::BandedMat *siconos::algebra::SimpleMatrix::banded(unsigned int,
                                                                    unsigned int) const {
  if (_num != UblasType::BANDED) THROW_EXCEPTION("the current matrix is not a Banded matrix");

  return mat.Banded;
}

siconos::algebra::ZeroMat *siconos::algebra::SimpleMatrix::zero_mat(unsigned int,
                                                                    unsigned int) const {
  if (_num != UblasType::ZERO) THROW_EXCEPTION("the current matrix is not a Zero matrix");

  return mat.Zero;
}

siconos::algebra::IdentityMat *siconos::algebra::SimpleMatrix::identity(unsigned int,
                                                                        unsigned int) const {
  if (_num != UblasType::IDENTITY)
    THROW_EXCEPTION("the current matrix is not a Identity matrix");

  return mat.Identity;
}

double *siconos::algebra::SimpleMatrix::getArray(unsigned int, unsigned int) const {
  if (_num == UblasType::SPARSE) THROW_EXCEPTION("not yet implemented for sparse matrix.");

  if (_num == UblasType::DENSE)
    return (((*mat.Dense).data()).data());
  else if (_num == UblasType::TRIANGULAR)
    return &(((*mat.Triang).data())[0]);
  else if (_num == UblasType::SYMMETRIC)
    return &(((*mat.Sym).data())[0]);
  else if (_num == UblasType::ZERO) {
    ZeroMat::iterator1 it = (*mat.Zero).begin1();
    return const_cast<double *>(&(*it));
  } else if (_num == UblasType::IDENTITY) {
    IdentityMat::iterator1 it = (*mat.Identity).begin1();
    return const_cast<double *>(&(*it));
  } else
    return &(((*mat.Banded).data())[0]);
}

// ===========================
//       fill matrix
// ===========================

void siconos::algebra::SimpleMatrix::zero() {
  unsigned int size1 = size(0);
  unsigned int size2 = size(1);
  if (_num == UblasType::DENSE)
    *mat.Dense = ublas::zero_matrix<double>(size1, size2);
  else if (_num == UblasType::TRIANGULAR)
    *mat.Triang = ublas::zero_matrix<double>(size1, size2);

  else if (_num == UblasType::SYMMETRIC)
    *mat.Sym = ublas::zero_matrix<double>(size1, size2);

  else if (_num == UblasType::SPARSE)
    *mat.Sparse = ublas::zero_matrix<double>(size1, size2);

  else if (_num == UblasType::SPARSE_COORDINATE)
    *mat.SparseCoordinate = ublas::zero_matrix<double>(size1, size2);

  else if (_num == UblasType::BANDED)
    *mat.Banded = ublas::zero_matrix<double>(size1, size2);

  else if (_num == UblasType::IDENTITY)
    THROW_EXCEPTION("you can not set to zero a matrix of type Identity!.");
  resetFactorizationFlags();
  // if _num == UblasType::ZERO: nothing
}

void siconos::algebra::SimpleMatrix::randomize() {
  if (_num == UblasType::DENSE)
    siconos::algebra::internal::randomize(*mat.Dense);
  else if (_num == UblasType::TRIANGULAR)
    siconos::algebra::internal::randomize(*mat.Triang);
  else if (_num == UblasType::SYMMETRIC)
    siconos::algebra::internal::randomize(*mat.Sym);
  else if (_num == UblasType::SPARSE)
    siconos::algebra::internal::randomize(*mat.Sparse);
  else if (_num == UblasType::BANDED)
    siconos::algebra::internal::randomize(*mat.Banded);
  else
    THROW_EXCEPTION("only implemented for dense matrices.");
  resetFactorizationFlags();
}

void siconos::algebra::SimpleMatrix::eye() {
  unsigned int size1 = size(0);
  unsigned int size2 = size(1);
  if (_num == UblasType::DENSE)
    *mat.Dense = ublas::identity_matrix<double>(size1, size2);

  else if (_num == UblasType::TRIANGULAR)
    *mat.Triang = ublas::identity_matrix<double>(size1, size2);

  else if (_num == UblasType::SYMMETRIC)
    *mat.Sym = ublas::identity_matrix<double>(size1, size2);

  else if (_num == UblasType::SPARSE)
    *mat.Sparse = ublas::identity_matrix<double>(size1, size2);

  else if (_num == UblasType::BANDED)
    *mat.Banded = ublas::identity_matrix<double>(size1, size2);

  else if (_num == UblasType::ZERO)
    THROW_EXCEPTION("you can not set to identity a matrix of type Zero!.");
  resetFactorizationFlags();
}

unsigned int siconos::algebra::SimpleMatrix::size(unsigned int index) const {
  if (_num == UblasType::DENSE) {
    if (index == 0)
      return (*mat.Dense).size1();
    else
      return (*mat.Dense).size2();
  } else if (_num == UblasType::TRIANGULAR) {
    if (index == 0)
      return (*mat.Triang).size1();
    else
      return (*mat.Triang).size2();
  } else if (_num == UblasType::SYMMETRIC) {
    if (index == 0)
      return (*mat.Sym).size1();
    else
      return (*mat.Sym).size2();
  } else if (_num == UblasType::SPARSE) {
    if (index == 0)
      return (*mat.Sparse).size1();
    else
      return (*mat.Sparse).size2();
  } else if (_num == UblasType::SPARSE_COORDINATE) {
    if (index == 0)
      return (*mat.SparseCoordinate).size1();
    else
      return (*mat.SparseCoordinate).size2();
  } else if (_num == UblasType::BANDED) {
    if (index == 0)
      return (*mat.Banded).size1();
    else
      return (*mat.Banded).size2();
  } else if (_num == UblasType::ZERO) {
    if (index == 0)
      return (*mat.Zero).size1();
    else
      return (*mat.Zero).size2();
  } else if (_num == UblasType::IDENTITY) {
    if (index == 0)
      return (*mat.Identity).size1();
    else
      return (*mat.Identity).size2();
  } else
    return 0;
};

//=======================
// set matrix dimension
//=======================

void siconos::algebra::SimpleMatrix::resize(unsigned int row, unsigned int col,
                                            unsigned int lower, unsigned int upper,
                                            bool preserve) {
  if (_num == UblasType::DENSE) {
    (*mat.Dense).resize(row, col, preserve);
  } else if (_num == UblasType::TRIANGULAR) {
    (*mat.Triang).resize(row, col, preserve);
  } else if (_num == UblasType::SYMMETRIC) {
    (*mat.Sym).resize(row, col, preserve);
  } else if (_num == UblasType::SPARSE) {
    (*mat.Sparse).resize(row, col, preserve);
  } else if (_num == UblasType::SPARSE_COORDINATE) {
    (*mat.SparseCoordinate).resize(row, col, preserve);
  } else if (_num == UblasType::BANDED) {
    (*mat.Banded).resize(row, col, lower, upper, preserve);
  } else if (_num == UblasType::ZERO) {
    (*mat.Zero).resize(row, col, preserve);
  } else if (_num == UblasType::IDENTITY) {
    (*mat.Identity).resize(row, col, preserve);
  }
  resetFactorizationFlags();
}

//=====================
// screen display
//=====================

void siconos::algebra::SimpleMatrix::display() const {
  if (size(0) == 0 || size(1) == 0) {
    std::cout << "siconos::algebra::SimpleMatrix::display(): empty matrix" << "\n";
  }
  std::cout << "SimpleMatrix storage type - num = "
            << static_cast<std::underlying_type<UblasType>::type>(_num) << "\n";
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  if (_num == UblasType::DENSE) {
    siconos::algebra::boost_bindings::print_m(*mat.Dense);
    // std::cout << *mat.Dense << "\n";
  } else if (_num == UblasType::TRIANGULAR)
    std::cout << *mat.Triang << "\n";
  else if (_num == UblasType::SYMMETRIC)
    std::cout << *mat.Sym << "\n";
  else if (_num == UblasType::SPARSE) {
    std::cout << "non zero element (nnz) = " << mat.Sparse->nnz() << "\n";

    std::cout << *mat.Sparse << "\n";
  } else if (_num == UblasType::SPARSE_COORDINATE) {
    std::cout << *mat.SparseCoordinate << "\n";
  } else if (_num == UblasType::BANDED)
    std::cout << *mat.Banded << "\n";
  else if (_num == UblasType::ZERO)
    std::cout << *mat.Zero << "\n";
  else if (_num == UblasType::IDENTITY)
    std::cout << *mat.Identity << "\n";
}
void siconos::algebra::SimpleMatrix::displayExpert(bool brief) const {
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  if (size(0) == 0 || size(1) == 0) {
    std::cout << "siconos::algebra::SimpleMatrix::display(): empty matrix" << "\n";
  }
  std::cout << "SimpleMatrix storage type - num = "
            << static_cast<std::underlying_type<UblasType>::type>(_num) << "\n";
  if (_num == UblasType::DENSE) {
    siconos::algebra::boost_bindings::print_m(*mat.Dense);
    // std::cout << *mat.Dense << "\n";
  } else if (_num == UblasType::TRIANGULAR)
    std::cout << *mat.Triang << "\n";
  else if (_num == UblasType::SYMMETRIC)
    std::cout << *mat.Sym << "\n";
  else if (_num == UblasType::SPARSE) {
    std::cout << "non zero element (nnz) = " << mat.Sparse->nnz() << "\n";
    std::cout << "non zero element (nnz_capacity) = " << mat.Sparse->nnz_capacity() << "\n";
    std::cout << "filled1 = " << mat.Sparse->filled1() << "\n";
    std::cout << "filled2 = " << mat.Sparse->filled2() << "\n";

    std::cout << "index_data1 = [ ";
    size_t i = 0;
    for (i = 0; i < mat.Sparse->filled1() - 1; i++) {
      std::cout << mat.Sparse->index1_data()[i] << ", ";
    }
    std::cout << mat.Sparse->index1_data()[i] << "]" << "\n";

    std::cout << "index_data2 = [";
    for (i = 0; i < mat.Sparse->filled2() - 1; i++) {
      std::cout << mat.Sparse->index2_data()[i] << ", ";
    }
    std::cout << mat.Sparse->index2_data()[i] << "]" << "\n";

    std::cout << "value_data = [";
    for (i = 0; i < mat.Sparse->filled2() - 1; i++) {
      std::cout << mat.Sparse->value_data()[i] << ", ";
    }
    std::cout << mat.Sparse->value_data()[i] << "]" << "\n";

    std::cout << *mat.Sparse << "\n";
  } else if (_num == UblasType::SPARSE_COORDINATE) {
    std::cout << "non zero element (nnz) = " << mat.SparseCoordinate->nnz() << "\n";

    for (size_t i = 0; i < mat.SparseCoordinate->nnz(); ++i) {
      // std::cout << i << "\n";
      std::cout << "M(" << mat.SparseCoordinate->index1_data()[i] << ", ";
      std::cout << mat.SparseCoordinate->index2_data()[i] << ") =  ";
      std::cout << mat.SparseCoordinate->value_data()[i] << "\n";
    }
  } else if (_num == UblasType::BANDED)
    std::cout << *mat.Banded << "\n";
  else if (_num == UblasType::ZERO)
    std::cout << *mat.Zero << "\n";
  else if (_num == UblasType::IDENTITY)
    std::cout << *mat.Identity << "\n";
}

//=====================
// convert to an ostream
//=====================

std::ostream &siconos::algebra::operator<<(std::ostream &os, const SimpleMatrix &sm) {
  if (sm._num == UblasType::DENSE)
    os << *sm.mat.Dense;
  else if (sm._num == UblasType::TRIANGULAR)
    os << *sm.mat.Triang;
  else if (sm._num == UblasType::SYMMETRIC)
    os << *sm.mat.Sym;
  else if (sm._num == UblasType::SPARSE)
    os << *sm.mat.Sparse;
  else if (sm._num == UblasType::BANDED)
    os << *sm.mat.Banded;
  else if (sm._num == UblasType::ZERO)
    os << *sm.mat.Zero;
  else if (sm._num == UblasType::IDENTITY)
    os << *sm.mat.Identity;
  return os;
}

void siconos::algebra::SimpleMatrix::assign(const SimpleMatrix &smat) {
  if (_num != UblasType::SPARSE)
    THROW_EXCEPTION("Assign only implemented for UblasType::SPARSE matrix.");

  // switch (_num) {
  //   case UblasType::SPARSE:

  switch (smat.num()) {
    case UblasType::SPARSE:
      mat.Sparse->assign(smat.getSparse());
      break;
    default:
      break;
  }
  // default:
  //     THROW_EXCEPTION("do not know how to assign for the given storage type ");
  //     break;
  // }
}

// void prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, bool init)
// {
//   assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
//   if(init)
//     y.zero();
//   unsigned int startRow = 0;
//   unsigned int startCol = 0;
//   // In private_addprod, the sum of all blocks of x, x[i], is computed: y = Sum_i (subA
//   x[i]), with subA a submatrix of A,
//   // starting from position startRow in rows and startCol in columns.
//   // private_prod takes also into account the fact that each block of x can also be a block.
//   VectorOfVectors::const_iterator it;
//   for(it = x.begin(); it != x.end(); ++it)
//   {
//     private_addprod(A, startRow, startCol, **it, y);
//     startCol += (*it)->size();
//   }
// }

// void private_addprod(const SiconosMatrix& A, unsigned int startRow, unsigned int startCol,
// const BlockVector& x, SiconosVector& y)
// {
//   assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
//   assert(!A.isBlock() && "private_addprod(A,start,x,y) error: not yet implemented for block
//   matrix."); VectorOfVectors::const_iterator it; unsigned int startColBis = startCol; for(it
//   = x.begin(); it != x.end(); ++it)
//   {
//     private_addprod(A, startRow, startColBis, **it, y);
//     startColBis += (*it)->size();
//   }

// }

// // x block, y siconos
// void private_prod(const SiconosMatrix& A, unsigned int startRow, const BlockVector& x,
// SiconosVector& y, bool init)
// {
//   assert(!(A.isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");
//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el.
//   of index (row) startRow and startRow + sizeY if(init)  // y = subA * x , else y += subA *
//   x
//     y.zero();
//   private_addprod(A, startRow, 0, x, y);
// }

// // x and y blocks
// void private_prod(std::shared_ptr<const SiconosMatrix> A, const unsigned int startRow,
// std::shared_ptr<const BlockVector> x, std::shared_ptr<BlockVector> y, bool init)
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
// void private_prod(std::shared_ptr<const SiconosMatrix> A, const unsigned int startRow,
// std::shared_ptr<const SiconosVector> x, std::shared_ptr<BlockVector> y, bool init)
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

// void private_addprod(std::shared_ptr<const BlockVector> x, std::shared_ptr<const
// SiconosMatrix> A, unsigned int startRow, unsigned int startCol,
// std::shared_ptr<SiconosVector> y)
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

// void private_prod(std::shared_ptr<const SiconosVector> x, std::shared_ptr<const
// SiconosMatrix> A, unsigned int startCol, std::shared_ptr<BlockVector>  y, bool init)
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

// void private_prod(std::shared_ptr<const BlockVector> x, std::shared_ptr<const SiconosMatrix>
// A, unsigned int startCol, std::shared_ptr<SiconosVector>  y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A),
//   between el. of A of index (col) startCol and startCol + sizeY if(init)  // y = subA * x ,
//   else y += subA * x
//     y->zero();
//   private_addprod(x, A, startCol, 0, y);

// }

// void private_prod(std::shared_ptr<const BlockVector> x, std::shared_ptr<const SiconosMatrix>
// A, unsigned int startCol, std::shared_ptr<BlockVector>  y, bool init)
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

// void private_addprod(double a, std::shared_ptr<const SiconosMatrix> A, unsigned int
// startRow, unsigned int startCol, std::shared_ptr<const SiconosVector> x,
// std::shared_ptr<SiconosVector> y)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   if(A->isBlock())
//      THROW_EXCEPTION("not yet implemented for block matrix.");

//   // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and
//   between columns startCol and (startCol+sizeX).
//   // Then computation of y = subA*x + y.
//   auto numA = A->num();
//   auto numY = y->num();
//   auto numX = x->num();
//   unsigned int sizeX = x->size();
//   unsigned int sizeY = y->size();

//   if(numX != numY)
//      THROW_EXCEPTION("not yet implemented for x and y of different types.");

//   if(numY == 1 && numX == 1)
//   {

//     assert(y->dense() != x->dense());

//     if(numA == 1)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->dense(), startRow, startRow +
//       sizeY, startCol, startCol + sizeX), *x->dense());
//     else if(numA == 2)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->triang(), startRow, startRow +
//       sizeY, startCol, startCol + sizeX), *x->dense());
//     else if(numA == 3)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->sym(), startRow, startRow +
//       sizeY, startCol, startCol + sizeX), *x->dense());
//     else if(numA == 4)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->sparse(), startRow, startRow +
//       sizeY, startCol, startCol + sizeX), *x->dense());
//     else //if(numA==5)
//       noalias(*y->dense()) += a * prod(ublas::subrange(*A->banded(), startRow, startRow +
//       sizeY, startCol, startCol + sizeX), *x->dense());
//   }
//   else // x and y sparse
//   {
//     if(numA == 4)
//       *y->sparse() += a * prod(ublas::subrange(*A->sparse(), startRow, startRow + sizeY,
//       startCol, startCol + sizeX), *x->sparse());
//     else
//        THROW_EXCEPTION("not yet implemented for x, y  sparse and A not sparse.");
//   }

// }

// void private_prod(double a, std::shared_ptr<const SiconosMatrix> A, unsigned int startRow,
// std::shared_ptr<const SiconosVector> x, std::shared_ptr<SiconosVector>  y, bool init)
// {
//   assert(!(A->isPLUFactorizedInPlace()) && "A is PLUFactorizedInPlace in prod !!");

//   // Computes y = subA *x (or += if init = false), subA being a sub-matrix of A, between el.
//   of index (row) startRow and startRow + sizeY

//   if(init)  // y = subA * x , else y += subA * x
//     y->zero();
//   private_addprod(a, A, startRow, 0, x, y);

// }

unsigned siconos::algebra::SimpleMatrix::copyData(double *data) const {
  assert((_num == UblasType::DENSE) &&
         "SiconosMatrix::copyData : forbidden: the current matrix is not dense.");

  unsigned size = mat.Dense->size1() * mat.Dense->size2();
  bindings_blas::detail::copy(size, getArray(), 1, data, 1);
  return size;
}
