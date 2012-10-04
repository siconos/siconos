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
  \left\{
  \begin{array}{l}
  M \ z + q = w \\
  w_1=0 \\
  0 \le w_{2} \perp v \ge 0
  \end{array}
  \right.
  \text{ with }
  z=
  \left[
  \begin{array}{c}
  u\\
  v\\
  \end{array}
  \right]
  \text{ and }
  w=
  \left[
  \begin{array}{c}
  w_{1}\\
  w_{2}\\
  \end{array}
  \right]
  \f$

  \f$ u, w_{1}\f$ are vectors of size n.\n
  \f$ v, w_{2}\f$ are vectors of size m.


  Another storage is also possible for the MLCP problem:

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

  Note that all the algorithms are not available for the two options of storage M,q or A,B,C,D,a,b

  (see the functions/solvers list in MLCP_solvers.h)

*/

#include "NumericsMatrix.h"

/**  \struct MixedLinearComplementarityProblem MixedLinearComplementarityProblem.h
 *  \brief Structure that contains and defines  \ref MLCProblem
 *
 Find \f$(z,w)\f$ such that:\n

  \f$
  \left\{
  \begin{array}{l}
  M \ z + q = w \\
  w_1=0 \\
  0 \le w_{2} \perp v \ge 0
  \end{array}
  \right.
  \text{ with }
  z=
  \left[
  \begin{array}{c}
  u\\
  v\\
  \end{array}
  \right]
  \text{ and }
  w=
  \left[
  \begin{array}{c}
  w_{1}\\
  w_{2}\\
  \end{array}
  \right]
  \f$

  \f$ u, w_{1}\f$ are vectors of size n.\n
  \f$ v, w_{2}\f$ are vectors of size m.

 * See \ref MLCProblem for more details.
 */
typedef struct
{
  int isStorageType1; /**< boolean for storageType1 1 if the problem is saved using (M,q), 0 otherwise */
  int isStorageType2; /**< boolean for storageType1 1 if the problem is saved using (A,B,C,D,a,b), 0 otherwise*/
  int n; /**< number of the linear constraints*/
  int m; /**< number of the complementarity constraints*/
  int * blocksRows;  /**< The rows from blocksRows[i] to blocksRows[i+1]-1 forms a block of equalities iif bloksIsComp[i]=0,
                       else the block is a complementarity block. */
  int * blocksIsComp; /**< if bloksIsComp[i]=0, then block i formed by the  rows from blocksRows[i] to blocksRows[i+1]-1 is an equality block
                       else the block is a complementarity block.
                       The number of total blocks is given by NbBlocks such that blocksRows[NbBlocks] = n+m*/
  NumericsMatrix* M; /**< M matrix of the MLCP */
  NumericsMatrix* Bblock; /**< Bblock  ?*/
  double *q; /**< q vector of the MLCP */
  double *A; /**< A matrix of the MLCP */
  double *B; /**< B matrix of the MLCP */
  double *C; /**< C matrix of the MLCP */
  double *D; /**< D matrix of the MLCP */
  double *a; /**< a vector of the MLCP */
  double *b; /**< b vector of the MLCP */
} MixedLinearComplementarityProblem;


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void displayMLCP(MixedLinearComplementarityProblem* p);

  /** \fn  int mixedLinearComplementarity_newFromFile(MixedLinearComplementarityProblem* problem, FILE* file)
   *  \brief function to read and create a MixedLinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a MixedLinearComplementarityProblem to create
   *  \param file pointer to a FILE
   */
  int mixedLinearComplementarity_newFromFile(MixedLinearComplementarityProblem* problem, FILE* MLCPfile);

  /** \fn  int mixedLinearComplementarity_newFromFilename(MixedLinearComplementarityProblem* problem, FILE* MLCPfile)
   *  \brief function to read and create a MixedLinearComplementarityProblem
   *   from a file
   *  \param problem pointer to a MixedLinearComplementarityProblem to create
   *  \param filename that contains the mlcp
   */
  int mixedLinearComplementarity_newFromFilename(MixedLinearComplementarityProblem* problem, char* filename);

  /** \fn  void freeMixedLinearComplementarityProblem(LinearComplementarityProblem* problem)
   *  \brief function to delete a LinearComplementarityProblem
   *  \param problem  pointer to a LinearComplementarityProblem to delete
   */
  void freeMixedLinearComplementarityProblem(MixedLinearComplementarityProblem* problem);



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

