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

/*! \file SimpleLCP.h

  \brief Interface for Path(Ferris) Solver
  \author Olivier Bonnefond
*/

/* \page FerrisPath Path (Ferris) Solver Interface

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

#include "SiconosConfig.h"
#include <stdint.h>

#ifdef HAVE_PATHFERRIS
#include "PATH_SDK/include/Types.h"
#else
typedef void MCP_Termination;
#endif

//const unsigned short int *__ctype_b;
//const int32_t *__ctype_tolower ;

/**


 */
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
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

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

