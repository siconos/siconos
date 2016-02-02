/* Siconos-Kernel, Copyright INRIA 2005-2012.
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


/** ! \file SiconosAlgebraTypeDef.hpp
    \brief Header file for Siconos Algebra objects

    This file provides typedef for matrix and vector objects, const values and so on ...
    \author SICONOS Development Team - copyright INRIA
    \date (creation) 10/17/2006
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
#ifndef SICONOS_CXXVERSION
#error SICONOS_CXXVERSION is not defined in SiconosConfig.h
#endif

#if SICONOS_CXXVERSION >= 201103L
#include <array>
#else
#include <boost/array.hpp>
#endif

#include <complex>
#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

#include "SiconosVisitor.hpp"

/** Const from old version of SiconosVector - To be reviewed */
const char N_DOUBLE_PRECISION[] = "%1.52e "; // double mantisse precision /!\ DEPENDS ON MACHINE
const unsigned int M_MAXSIZEFORDISPLAY = 10;
const std::string DEFAULT_FORMAT = "ascii";

/** Siconos::UBLAS_TYPE is an enumerated type of Siconos::DENSE, TRIANGULAR, SYMMETRIC,
    SPARSE, BANDED. It is used to describe the type of matrix or
    vector we want to construct.
 */
namespace Siconos
{
enum UBLAS_TYPE {DENSE = 1, TRIANGULAR, SYMMETRIC, SPARSE, BANDED, ZERO, IDENTITY};
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

namespace ublas = boost::numeric::ublas;

/** Various matrix types available in Siconos **/

/** DenseMat is a typedef of boost::ublas::numeric::matrix<double, column_major, std::vector<double> >
 */
typedef ublas::matrix<double, ublas::column_major, std::vector<double> > DenseMat;
TYPEDEF_SPTR(DenseMat)

//typedef ublas::matrix<double, ublas::column_major, ublas::bounded_array<double, 10000> > DenseMat;
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

/** SparseMat is a typedef of boost::ublas::numeric::mapped_matrix<double>
 */
typedef ublas::compressed_matrix<double, ublas::column_major, 0, Index > SparseMat;
TYPEDEF_SPTR(SparseMat)

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
