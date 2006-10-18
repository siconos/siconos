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


/** Header file for Siconos Algebra objects
 *  \brief This file provides typedef for matrix and vector objects, const values and so on ...
 *  \author SICONOS Development Team - copyright INRIA
 *  \date (creation) 10/17/2006
 */

#ifndef SiconosAlgebra
#define SiconosAlgebra


#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <cassert>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>


using namespace boost::numeric::ublas;

/** type of object used to save indices */
typedef std::vector<unsigned int> Index;

/** Various matrix types available in Siconos **/

/**\brief DenseMat is a typedef of boost::ublas::numeric::matrix<double, row_major, std::vector<double> >
 */
typedef matrix<double, row_major, std::vector<double> > DenseMat;

/**\brief TriangMat is a typedef of boost::ublas::numeric::triangular_matrix<double, upper, row_major, std::vector<double> >
 */
typedef triangular_matrix<double, upper, row_major, std::vector<double> > TriangMat;

/**\brief SymMat is a typedef of boost::ublas::numeric::symmetric_matrix<double, upper, row_major, std::vector<double> >
 */
typedef symmetric_matrix<double, upper, row_major, std::vector<double> > SymMat;

/**\brief BandedMat is a typedef of boost::ublas::numeric::banded_matrix<double, row_major, std::vector<double> >
 */
typedef banded_matrix<double, row_major, std::vector<double> > BandedMat;

/**\brief SparseMat is a typedef of boost::ublas::numeric::mapped_matrix<double>
 */
typedef compressed_matrix<double, row_major, 0, Index, std::vector<double> > SparseMat;

/** Various vector types available in Siconos **/

/** \brief DenseVect is a typedef of boost::ublas::numeric::vector<double, std::vector<double> >
 */
typedef boost::numeric::ublas::vector<double, std::vector<double> > DenseVect;

/** \brief SparseVect is a typedef of boost::ublas::numeric::mapped<double>
 */
typedef compressed_vector<double, 0, Index, std::vector<double> > SparseVect;

/** value used to compare matrices. Matrices A and B are equal when (A-B).normInf()<tolerance. */
const double tolerance = 1e-14;

/** \brief TYP is an enumerated type of DENSE, TRIANGULAR, SYMMETRIC, SPARSE, BANDED. TYP is used to describe the type of matrix or vector we want to construct.
 */
enum TYP {DENSE = 1, TRIANGULAR, SYMMETRIC, SPARSE, BANDED};

// Notes:
// Vector definition in boost: vector<T,A> see http://www.boost.org/libs/numeric/ublas/doc/vector.htm
// T: value type
// A: storage type

/** Objects used to define block matrices and vectors:*/

class MySiconosMatrix;
/**\brief block-matrix: a collection of pointers to matrices ; blocksMat is a typedef of boost::ublas::numeric::mapped_matrix<MySiconosMatrix* > */
typedef compressed_matrix<MySiconosMatrix*> BlocksMat;

class MySiconosVector;
/**\brief block-vector: a collection of pointers to vectors;  blocksVect is a typedef of boost::ublas::numeric::mapped_matrix<MySiconosMatrix* > */
typedef std::vector<MySiconosVector*> BlocksVect;
/** the relative iterators */
typedef BlocksVect::iterator BlockVectIterator;
typedef BlocksVect::const_iterator ConstBlockVectIterator;
#endif

