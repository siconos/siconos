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

/*! \file SimpleLCP.h

  \brief Interface for Path(Ferris) Solver
  \author Olivier Bonnefond
*/

/*! \page FerrisPath Path (Ferris) Solver Interface

  The problem solved is:

  \f$ Mx + q  \perp l <= x <= u \f$

  M is a square matrix.  M should also be a positive semidefinite matrix.

  The LCP uses the (i, j, data) format for inputing sparse matrice (M).
  We adopt the FORTRAN convention that the indices begin at 1.

  Here is a description of each of the arguments to the LCP_Create functions:

  - variables - the number of variables in the problem
  - m_nnz - the number of nonzeros in the M matrix
  - m_i - a vector of size m_nnz containing the row indices for M
  - m_j - a vector of size m_nnz containing the column indices for M
  - m_ij - a vector of size m_nnz containing the data for M
  - q - a vector of size variables
  - lb - a vector of size variables containing the lower bounds on x
  - ub - a vector of size variables containing the upper bounds on x
*/

#ifndef SIMPLELCP_H
#define SIMPLELCP_H

#include "include/Types.h"
/*because libpath46.so*/
const unsigned short int *__ctype_b;
const __int32_t *__ctype_tolower ;

/**


 */
void SimpleLCP(int variables,
               int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
               double *lb, double *ub,
               MCP_Termination *status, double *z);



void printLCP(int variables,
              int m_nnz, int *m_i, int *m_j, double *m_ij, double *q,
              double *lb, double *ub);

int nbNonNulElems(int n, double *M, double tol);
void FortranToPathSparse(int n, double *M, double tol, int *m_i, int *m_j, double *m_ij);
void ABCDtoM(int n , int m, double *A , double *B , double *C , double *D , double *a, double *b, double *M, double *q);
#endif

