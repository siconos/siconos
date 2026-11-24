// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NCPPath_H
#define NCPPath_H

#include "PathAlgebra.h"
#include "Standalone_Path.h"
/*!\file NCP_Path.h

  \brief Interface to Path Solver for NCP problems

  Solves the following Nonlinear Complementarity Problem with Path Solver:
  \f{eqnarray*}
  0 \le F(z) \perp z \ge 0 \\
  \f}

*/

#if defined(__cplusplus)
extern "C" {
#endif

/** Path solver for NCP problem
    \param n size of the vector z
    \param z vector
    \param F pointer to function used to compute \f$ F(z) \f$
    \param jacobianF pointer to function used to compute \f$ \nabla_zF(z) \f$
    \param iparam vector of int parameters (useless at the time)
    \param dparam vector of double parameters (useless at the time)
    \return 0 if successfull
*/
int NCP_Path(int n, double* z, FuncEvalPtr F, JacEvalPtr jacobianF, int* iparam,
             double* dparam);

#if defined(__cplusplus)
}
#endif

#endif
