/* Siconos-Numerics, Copyright INRIA 2005-2010.
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

#ifndef SparseMatrix_H
#define SparseMatrix_H

/*! \page SparseMatrixPage Matrix Storage in Numerics

*/




/*!\file SparseMatrix.h
  \brief Structure definition and functions related to sparse matrix storage in Numerics
*/
#include <stdio.h>

typedef struct   /* matrix in compressed row/column or triplet form */
{
  int nzmax ;     /* maximum number of entries */
  int m ;     /* number of rows */
  int n ;     /* number of columns */
  int *p ;      /* compressed: row (size m+1) or column (size n+1) pointers; triplet: row indices (size nz) */
  int *i ;      /* compressed: column or row indices, size nzmax; triplet: column indices (size nz) */
  double *x ;     /* numerical values, size nzmax */
  int nz ;      /* # of entries in triplet matrix;
           -1 for compressed columns;
           -2 for compressed rows */
} SparseMatrix;

#ifdef __cplusplus
extern "C"
{
#endif
  /** Screen display of the matrix content
      \param M the matrix to be displayed
   */
  void printSparseMat(const SparseMatrix* const M);

#ifdef __cplusplus
}
#endif

#endif
