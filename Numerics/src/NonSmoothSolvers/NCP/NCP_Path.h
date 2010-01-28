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

#ifndef NCPPath_H
#define NCPPath_H

#include "Standalone_Path.h"
#include "PathAlgebra.h"
/*!\file NCP_Path.h

  \brief Interface to Path Solver for NCP problems

  Solves the following Nonlinear Complementarity Problem with Path Solver:
  \f{eqnarray*}
  0 \le F(z) \perp z \ge 0 \\
  \f}

  \author Franck Perignon, 21/05/2008
*/

#ifdef __cplusplus
extern "C"
{
#endif

  /** Path solver for NCP problem
      \param n size of the vector z
      \param z vector
      \param F pointer to function used to compute \f$ F(z) \f$
      \param jacobianF pointer to function used to compute \f$ \nabla_zF(z) \f$
      \param iparam vector of int parameters (useless at the time)
      \param dparam vector of double parameters (useless at the time)
  */
  int NCP_Path(int n, double* z, FuncEvalPtr F, JacEvalPtr jacobianF, int* iparam, double* dparam);

#ifdef __cplusplus
}
#endif

#endif
