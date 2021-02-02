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


/** ! \file SiconosAlgebraTypeDef.hpp
    \brief Header file for Siconos Algebra objects

    This file provides typedef for matrix and vector objects, const values and so on ...
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
#include <limits>
#include <boost/numeric/ublas/fwd.hpp>

#include "SiconosConfig.h"
#include <array>

#include <complex>
#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

#include "SiconosVisitor.hpp"

/** Const from old version of SiconosVector - To be reviewed */
const char N_DOUBLE_PRECISION[] = "%1.52e "; // double mantisse precision /!\ DEPENDS ON MACHINE
const unsigned int M_MAXSIZEFORDISPLAY = 10;
const std::string DEFAULT_FORMAT = "ascii";

namespace Siconos
{
  /** Siconos::UBLAS_TYPE is an enumerated type of Siconos::DENSE, TRIANGULAR, SYMMETRIC,
      SPARSE, BANDED. It is used to describe the type of matrix or
      vector we want to construct.
  */
  enum UBLAS_TYPE {
    BLOCK =0,
    /** id for dense matrix or vector */
    DENSE = 1, 
    /** id for triangular matrix */
    TRIANGULAR=2,
    /** id for symmetric matrix */
    SYMMETRIC=3,
    /** id for sparse matrix or vector */
    SPARSE=4,
    /** id for banded matrix */
    BANDED=5,
    /** id for zero matrix */
    ZERO=6,
    /** id for identity matrix */
    IDENTITY=7,
    /** id for sparse matrix or vector */
    SPARSE_COORDINATE=8};
}
// Notes:
// Vector definition in boost: vector<T,A> see http://www.boost.org/libs/numeric/ublas/doc/vector.htm
// T: value type
// A: storage type

/** Objects used to define block matrices and vectors:*/



/** Some containers for vectors - Used for example to handle x and its
    derivatives in DynamicalSystem. */
typedef std::vector<SP::SiconosVector> VectorOfVectors;

/** Some containers for vectors - Used for example in Relation to compute y and r
 * when there are 2 DS*/
typedef std::vector<SP::BlockVector> VectorOfBlockVectors;

/** Some containers for matrices - Used for example to handle the
    various jacobian in LagrangianDS. */
typedef std::vector<SP::SiconosMatrix> VectorOfMatrices;


/** Some containers for matrices - Used for example to handle the
    various jacobian in LagrangianDS. */
typedef std::vector<SP::SimpleMatrix> VectorOfSMatrices;

//typedef std::vector<SP::SimpleMatrix> VectorOfSimpleMatrices;

/** Iterator through vector of matrices */
typedef VectorOfMatrices::iterator VectorOfMatricesIterator;

/** const Iterator through vector of matrices */
typedef VectorOfMatrices::const_iterator VectorOfMatricesConstIterator;

/** type of object used to save indices */
typedef std::vector<std::size_t> Index;
TYPEDEF_SPTR(Index)
TYPEDEF_SPTR(VectorOfBlockVectors)
TYPEDEF_SPTR(VectorOfVectors)
TYPEDEF_SPTR(VectorOfMatrices)
TYPEDEF_SPTR(VectorOfSMatrices)

namespace ublas = boost::numeric::ublas;

/* Various matrix types available in Siconos */

/** DenseMat is a typedef of boost::ublas::numeric::matrix<double, column_major, std::vector<double> >  */
typedef ublas::matrix<double, ublas::column_major, std::vector<double> > DenseMat;

TYPEDEF_SPTR(DenseMat)

/** TriangMat is a typedef of boost::ublas::numeric::triangular_matrix<double, upper, column_major, std::vector<double> >
 */
typedef ublas::triangular_matrix<double, ublas::upper, ublas::column_major> TriangMat;

TYPEDEF_SPTR(TriangMat)

/** SymMat is a typedef of boost::ublas::numeric::symmetric_matrix<double, upper, column_major, std::vector<double> >
 */
typedef ublas::symmetric_matrix<double, ublas::upper, ublas::column_major> SymMat;

TYPEDEF_SPTR(SymMat)

/** BandedMat is a typedef of boost::ublas::numeric::banded_matrix<double, column_major, std::vector<double> >
 */
typedef ublas::banded_matrix<double, ublas::column_major > BandedMat;
TYPEDEF_SPTR(BandedMat)

/** SparseMat is a typedef of boost::ublas::numeric::compressed_matrix<double>
 */
typedef ublas::compressed_matrix<double, ublas::column_major, 0, Index > SparseMat;
TYPEDEF_SPTR(SparseMat)

/** SparseMat is a typedef of boost::ublas::numeric::mapped_matrix<double>
 */
typedef ublas::coordinate_matrix<double, ublas::column_major, 0, Index > SparseCoordinateMat;
TYPEDEF_SPTR(SparseCoordinateMat)

/** ZeroMat is a typedef of boost::ublas::numeric::zero_matrix, ie null matrix.
 */
typedef ublas::zero_matrix<double> ZeroMat;
TYPEDEF_SPTR(ZeroMat)

/** IdentityMat is a typedef of boost::ublas::identity_matrix ie identity matrix.
 */
typedef ublas::identity_matrix<double> IdentityMat;
TYPEDEF_SPTR(IdentityMat)

/** A collection of pointers to matrices ; blocksMat is a typedef of
    boost::ublas::numeric::mapped_matrix<SiconosMatrix* > */
typedef ublas::compressed_matrix<SP::SiconosMatrix> BlocksMat;
TYPEDEF_SPTR(BlocksMat)

/** Complex matrix
 */
typedef ublas::matrix<std::complex<double>, ublas::column_major> complex_matrix;
TYPEDEF_SPTR(complex_matrix)

/** Various vector types available in Siconos **/

/** DenseVect is a typedef of boost::ublas::numeric::vector<double, std::vector<double> >
 */
typedef ublas::vector<double, std::vector<double> > DenseVect;
TYPEDEF_SPTR(DenseVect)

/** SparseVect is a typedef of boost::ublas::numeric::mapped<double>
 */
typedef ublas::compressed_vector<double> SparseVect;
TYPEDEF_SPTR(SparseVect)

/** Complex vector
 */
typedef ublas::vector<std::complex<double> > complex_vector;
TYPEDEF_SPTR(complex_vector)


// Set this to use lapack::optimal_workspace where required in lapack routines.
#define USE_OPTIMAL_WORKSPACE

#endif
