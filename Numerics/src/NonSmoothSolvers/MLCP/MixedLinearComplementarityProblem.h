/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#ifndef MLCP_PROBLEM_H
#define MLCP_PROBLEM_H

/*!\file MixedLinearComplementarityProblem.h
  \brief Structure used to define a Mixed Linear Complementarity Problem

  \author Olivier Bonnefon
*/

/*! \page MLCProblem Mixed Linear Complementarity problems (MLCP)
  \section mlcpIntro The problem
  Find \f$(z,w)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ z + q = w \\
  \end{array}
  \right.
  \f$

  \f$
  z=
  \left\lbrace
  \begin{array}{c}
  u\\
  v\\
  \end{array}
  \right.
  \f$

  \f$
  w=
  \left\lbrace
  \begin{array}{c}
  w_{1}\\
  w_{2}\\
  \end{array}
  \right.
  \f$

  \f$ u, w_{1}\f$ are vectors of size n.\n
  \f$ v, w_{2}\f$ are vectors of size m.

  \f$ w_1=0\f$.\n
  \f$ 0 \le w_{2} \perp v \ge 0 \f$

  Old version (?):

  Try \f$(u,v,w)\f$ such that:\n
  \f$
  \left\lbrace
  \begin{array}{l}
  A u + Cv +a =0\\
  D u + Bv +b = w \\
  0 \le v \perp  w \ge 0\\
  \end{array}
  \right.
  \f$

  where  A is an (\f$ n \times n\f$ ) matrix, B is an (\f$ m \times m\f$ ) matrix,  C is an (\f$ n \times m\f$ ) matrix,\n
  D is an (\f$ m \times n\f$ ) matrix,    a and u is an (\f$ n \f$ ) vectors b,v and w is an (\f$ m \f$ ) vectors.

  \section mlcpSolversList Available solvers

  The solvers and their parameters are described in \ref MLCPSolvers . \n

  Use the generic function mlcp_solver(), to call one the the specific solvers listed below:
  - mlcp_pgs(), (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver
  - mlcp_rpgs(), (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver
  - mlcp_psor(), projected successive overrelaxation method
  - mlcp_rpsor(), regularized projected successive overrelaxation method
  - mlcp_path(), path solver
  - mlcp_enum(), enumeratif solver
  - mlcp_simplex(), solver based on the simplex algorithm
  - mlcp_direct_path(), use the last solution to find the solution. If it failed the mlcp_path is called.
  - mlcp_direct_enum(), use the last solution to find the solution. If it failed the mlcp_enum is called.
  - mlcp_direct_simplex(), use the last solution to find the solution. If it failed the mlcp_simplex is called.

  (see the functions/solvers list in MLCP_solvers.h)

*/

#include "NumericsMatrix.h"

/** Mixed Linear Complementarity Problem elements
    \param n dim of the linear constraints
    \param m dim of the Complementarity constraints
    \param M matrix of the MLCP
    \param q vector
    \param A matrix of the MLCP
    \param B matrix of the MLCP
    \param C matrix of the MLCP
    \param D matrix of the MLCP
    \param a vector
    \param b vector
    \param problemType 0 if the problem is saved using (M) and 1 if (A,B,C,D)
    \param blocksLine The lines blocksLine[i] to blocksLine[i+1] are equalities iif !bloksIsComp[i], else is a complementarity block.
    blocksLine[0]=0 ... blocksLine[last] = n+m
 */
typedef struct
{
  int n;
  int m;
  int * blocksLine;
  int * blocksIsComp;
  NumericsMatrix* M;
  NumericsMatrix* Bblock;
  double * q;
  double *A;
  double *B;
  double *C;
  double *D;
  double *a;
  double *b;
  int problemType;
} MixedLinearComplementarityProblem;

void displayMat(double * M, int Nblin, int Nbcol, int incLin);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  void displayMLCP(MixedLinearComplementarityProblem* p);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

