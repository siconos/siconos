//
// Copyright (C) INRIA 1999-2008
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//%
// @file kernel/Utils.hpp
// @author Rémy MOZUL
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): mozul@inria.fr
//
// @brief reducing and expanding functions for vecotr and matrices
//

#ifndef __Kernel_Utils_hpp
#define __Kernel_Utils_hpp

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

/**
 * Reduce the input matrix to the given indexes
 *
 * @param[in] mat (matrix)
 * @param[in] index1 (int)
 * @param[in] index2 (int)
 */
void reduceMatrix(matrix<double, column_major> & mat, vector<int> & index1, vector<int> & index2);

/**
 * Reduce the input matrix to the given indexes
 *
 * @param[in] mat (matrix)
 * @param[in] index1 (bool) <em> dim number of rows of mat </em>
 * @param[in] index2 (int)
 */
void reduceMatrix(matrix<double, column_major> & mat, vector<bool> & index1, vector<int> & index2);

/**
 * Reduce the input vector to the given index
 *
 * @param[in] vec (vector)
 * @param[in] index (int)
 *
 * @return reduced vector <em> dim size of index </em>
 */
vector<double> reduceVector(vector<double> & vec, vector<int> & index);

/**
 * Reduce the input vector to the given index
 *
 * @param[in/out] vec (vector)
 * @param[in] index (int)
 *
 * @return reduced vector <em> dim size of index </em>
 */
vector<double> reduceVector(vector<double, array_adaptor<double> > & vec, vector<int> & index);

/**
 * Reduce the input vector to the given index
 *
 * @param[in/out] vec (vector)
 * @param[in] index (bool) <em> dim size of vec </em>
 *
 */
void reduceVector(vector<double> & vec, vector<bool> & index);

/**
 * Reduce the input vector to the given index
 *
 * @param[in/out] vec (vector)
 * @param[in] index (bool) <em> dim size of vec </em>
 *
 */
void reduceVector(vector<bool> & vec, vector<bool> & index);

/**
 * Expand a vector2 into the given indexes of vector1
 *
 * @param[in/out] tab (double)
 * @param[in] vector (vector)
 * @param[in] index (int)
 */
void expandVector(double * tab, vector<double> & vector2, vector<int> & index);

/**
 * Expand a vector2 into the given indexes of vector1
 *
 * @param[in/out] tab (double)
 * @param[in] vector (vector)
 * @param[in] index (int) <em> dim size of tab </em>
 */
void expandVector(double * tab, vector<double> & vector, vector<bool> & index);

#endif /* __Kernel_Utils_hpp */
