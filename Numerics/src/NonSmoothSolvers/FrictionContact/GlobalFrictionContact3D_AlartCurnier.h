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

#ifndef GlobalFrictionContact3D_AlartCurnier_h
#define GlobalFrictionContact3D_AlartCurnier_h

#include "SolverOptions.h"
#include "GlobalFrictionContactProblem.h"

/** Compute y = alpha * a * x + beta * y 
like classical blas axpy but for cs-sparse
matrices.
Note FP : strange to find such function here, is it?
\param[in] alpha a scalar
\param[in] A address of matrix of double
\param[in] x address of vector of double
\param[in] beta a scalar
\param[in,out] y address of vector of double
*/
int cs_aaxpy(const double alpha, const cs *A, const double *x,
             const double beta, double *y);


int globalFrictionContact3D_AlartCurnier_setDefaultSolverOptions(
  SolverOptions* options);


void globalFrictionContact3D_sparseGlobalAlartCurnierInit(
  SolverOptions *SO);

void globalFrictionContact3D_AlartCurnier(
  GlobalFrictionContactProblem* problem,
  double *reaction, 
  double *velocity, 
  double *globalVelocity,
  int *info, 
  SolverOptions* options);

#endif
