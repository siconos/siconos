/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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


/** ! \file SiconosAlgebraTypeDef.hpp
    \brief List of typedefs related to vectors and matrices types (dense, sparse ...)
    used in Siconos.

*/

#ifndef SiconosAlgebraTypeDef
#define SiconosAlgebraTypeDef

// Make sure Fortran function have a calling convention compatible with gfortran
#ifndef BIND_FORTRAN_LOWERCASE_UNDERSCORE
#define BIND_FORTRAN_LOWERCASE_UNDERSCORE
#endif

// We do not want to link to any Boost lib, we use header only parts
#ifndef BOOST_ALL_NO_LIB
#define BOOST_ALL_NO_LIB
#endif

#include <vector>
#include <boost/numeric/ublas/fwd.hpp> // Boost forward declarations 

// #include <array>
#include <complex>

//#include "SiconosPointers.hpp"
#include "SiconosFwd.hpp"

// #include "SiconosVisitor.hpp"


namespace ublas = boost::numeric::ublas;

namespace Siconos
{
  /** Siconos::UBLAS_TYPE is an enumerated type of Siconos::DENSE, TRIANGULAR, SYMMETRIC,
      SPARSE, BANDED. It is used to describe the type of matrix or
      vector we want to construct.
  */
  enum UBLAS_TYPE {
    /** id for dense matrix or vector */
    DENSE = 1, 
    /** id for triangular matrix */
    TRIANGULAR,
    /** id for symmetric matrix */
    SYMMETRIC,
    /** id for sparse matrix or vector */
    SPARSE,
    /** id for banded matrix */
    BANDED,
    /** id for zero matrix */
    ZERO,
    /** id for identity matrix */
    IDENTITY,
    /** id for sparse matrix or vector */
    SPARSE_COORDINATE};
}
// Notes:
// Vector definition in boost: vector<T,A> see http://www.boost.org/libs/numeric/ublas/doc/vector.htm
// T: value type
// A: storage type

/** Some containers for vectors - Used for example to handle x and its
    derivatives in DynamicalSystem. */
typedef std::vector<SP::SiconosVector> VectorOfVectors;
TYPEDEF_SPTR(VectorOfVectors)

/** Some containers for vectors - Used for example in Relation to compute y and r
 * when there are 2 DS*/
typedef std::vector<SP::BlockVector> VectorOfBlockVectors;
TYPEDEF_SPTR(VectorOfBlockVectors)

/** Some containers for matrices - Used for example to handle the
    various jacobian in LagrangianDS. */
typedef std::vector<SP::SiconosMatrix> VectorOfMatrices;
TYPEDEF_SPTR(VectorOfMatrices)

/** Some containers for matrices - Used for example to handle the
    various jacobian in LagrangianDS. */
typedef std::vector<SP::SimpleMatrix> VectorOfSMatrices;
TYPEDEF_SPTR(VectorOfSMatrices)

/** Vector of indices (size_t type) */
typedef std::vector<std::size_t> Index;
TYPEDEF_SPTR(Index)

/** Vector of indices (int type) */
typedef std::vector<int> VInt;
TYPEDEF_SPTR(VInt)

/* ---  Definition of the matrices types available in Siconos  --- */

/** Dense matrix */
typedef ublas::matrix<double, ublas::column_major, std::vector<double> > DenseMat;
TYPEDEF_SPTR(DenseMat)

/** Triangular matrix */
typedef ublas::triangular_matrix<double, ublas::upper, ublas::column_major> TriangMat;
TYPEDEF_SPTR(TriangMat)

/** Symmetric matrix */
typedef ublas::symmetric_matrix<double, ublas::upper, ublas::column_major> SymMat;

TYPEDEF_SPTR(SymMat)

/** Banded matrix */
typedef ublas::banded_matrix<double, ublas::column_major > BandedMat;
TYPEDEF_SPTR(BandedMat)

/** Sparse matrix */
typedef ublas::compressed_matrix<double, ublas::column_major, 0, Index > SparseMat;
TYPEDEF_SPTR(SparseMat)

/** Sparse matrix (coordinate) */
typedef ublas::coordinate_matrix<double, ublas::column_major, 0, Index > SparseCoordinateMat;
TYPEDEF_SPTR(SparseCoordinateMat)

/** Null matrix */
typedef ublas::zero_matrix<double> ZeroMat;
TYPEDEF_SPTR(ZeroMat)

/** Identity matrix */
typedef ublas::identity_matrix<double> IdentityMat;
TYPEDEF_SPTR(IdentityMat)

/** Sparse block matrix */
typedef ublas::compressed_matrix<SP::SiconosMatrix> BlocksMat;
TYPEDEF_SPTR(BlocksMat)

/** Complex dense matrix  */
typedef ublas::matrix<std::complex<double>, ublas::column_major> complex_matrix;
TYPEDEF_SPTR(complex_matrix)

/* ---  Definition of the vectors types available in Siconos  --- */
/** Dense vector */
typedef ublas::vector<double, std::vector<double> > DenseVect;
TYPEDEF_SPTR(DenseVect)

/** Sparse vector */
typedef ublas::compressed_vector<double> SparseVect;
TYPEDEF_SPTR(SparseVect)

/** Complex dense vector */
typedef ublas::vector<std::complex<double> > complex_vector;
TYPEDEF_SPTR(complex_vector)

// Set this to use lapack::optimal_workspace where required in lapack routines.
#define USE_OPTIMAL_WORKSPACE

#endif
