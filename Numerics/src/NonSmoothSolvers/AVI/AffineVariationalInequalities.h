/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#ifndef AVI_PROBLEM_H
#define AVI_PROBLEM_H

/*!\file AffineVariationalInequalities.h
 * \brief Definitions for AVI
 *
 * \authors Olivier Huber
*/

/*! \page AVI Affine Variational Inequalities (AVI)
 * \section aviIntro The problem
 *  The Affine Variational Inequalities (AVI) is defined by
 *
 *   Find \f$(z,w)\f$ such that:\n
 *   \f{equation*}{
 *   \begin{cases}
 *   M \ z + q = w \\
 *   0 \le w \perp z \ge 0 \\
 *   \end{cases}
 *   \f}
 *
 * where \f$ w, z, q\f$ are vectors of size \f$n\f$ and \f$ M \f$ is a \f$n\times n\f$ matrix.
 *
 * The notation \f$x \perp y\f$ means that \f$x^Ty =0\f$. Inequalities involving vectors
 * are understood to hold component-wise.
 *
 * From more details on theory and analysis of AVI, we refer to
 *
 * R.W. Cottle, J.S. Pang, and R.E. Stone. <i>The Linear Complementarity Problem.</i> Academic Press, Inc., Boston, MA, 1992.
 *
 *  \section aviSolversList Available solvers

  The solvers and their parameters are described in \ref AVISolvers . \n
  Use the generic function AVI_driver() to call one the the specific solvers listed below:

    - AVI_caoferris(), direct solver for AVI based on pivoting method principle for degenerate problem.\n
  Choice of pivot variable is performed via lexicographic ordering
  Ref: "The Linear Complementarity Problem" Cottle, Pang, Stone (1992)\n

  (see also the functions/solvers list in AVI_Solvers.h and numbering in AVI_cst.h)

*/

#include "NumericsMatrix.h"
#include "polyhedron.h"

/** \struct AffineVariationalInequalities AffineVariationalInequalities.h
 *  \brief Structure that contains and defines  \ref AVI
 *
 *   Find \f$(z)\f$ such that:\n
 *   \f{equation*}{
 *   \begin{cases}
 *   M \ z + q = w \\
 *   0 \le w \perp z \ge 0 \\
 *   \end{cases}
 *   \f}
 *
 * where \f$ w, z, q\f$ are vectors of size \f$n\f$ and \f$ M \f$ is a \f$n\times n\f$ matrix.
 * See \ref AVI for more details.
 */
typedef struct
{
  unsigned int size; /**<  size of the problem */
  NumericsMatrix* M; /**< M matrix of the AVI (see the mathematical description)*/
  double* q; /**< vector of the AVI (see the mathematical description)*/
  double* d; /** Covering vector (optional) */
  Polyhedron* poly; /** Polyhedra where the solution has to belong */
} AffineVariationalInequalities;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** \fn void AVI_display(AffineVariationalInequalities* problem)
   *  \brief function to display a AffineVariationalInequalities
   *  \param  problem pointer to a AffineVariationalInequalities to display
   */
  void AVI_display(AffineVariationalInequalities* problem);

  /** \fn int AVI_printInFile(AffineVariationalInequalities*  problem, FILE* file)
   *  \brief function to write in a file a AffineVariationalInequalities
   *  \param problem pointer to a AffineVariationalInequalities to print
   *  \param file pointer to a FILE
   *  \return 1 if successfull
   */
  int AVI_printInFile(AffineVariationalInequalities* problem, FILE* file);

  /** \fn  int AVI_newFromFile(AffineVariationalInequalities* problem, FILE* file)
   *  \brief function to read and create a AffineVariationalInequalities
   *   from a file
   *  \param problem pointer to a AffineVariationalInequalities to create
   *  \param file pointer to a FILE
   *  \return 1 if successfull
   */
  int AVI_newFromFile(AffineVariationalInequalities* problem, FILE* file);

  /** \fn  int AVI_newFromFilename(AffineVariationalInequalities* problem, FILE* file)
   *  \brief function to read and create a AffineVariationalInequalities
   *   from a file
   *  \param problem pointer to a AffineVariationalInequalities to create
   *  \param filename that contains the AVI
   *  \return 1 if successfull
   */
  int AVI_newFromFilename(AffineVariationalInequalities* problem, char* filename);

  /** \fn  void freeAVI(AffineVariationalInequalities* problem)
   *  \brief function to delete a AffineVariationalInequalities
   *  \param problem  pointer to a AffineVariationalInequalities to delete
   */
  void freeAVI(AffineVariationalInequalities* problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

