/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef NCP_PROBLEM_H
#define NCP_PROBLEM_H

#include "NumericsFwd.h"  // for NonlinearComplementarityProblem, NumericsMa...

/*!\file NonlinearComplementarityProblem.h
 * Data structure to formalize a Nonlinear Complementarity Problem (NCP)
 *
 */

/** type for user defined function used to compute F and its jacobian.
 */
typedef void (*ptrFunctionNCP)(void* env, int n, double* z, double* F);
typedef void (*ptrFunctionJacNCP)(void* env, int n, double* z, NumericsMatrix* jacF);

/**
    The structure that defines a Nonlinear Complementarity Problem (NCP) :
    Find two vectors \f$ (z,w \in {{\mathrm{I\!R}}}^{n}) \f$ such that:

    \f[
    w &= F(z) \\ 0 &\le w \perp z \ge 0
    \f]
*/
struct NonlinearComplementarityProblem {
  unsigned n;               /**< size of the problem */
  ptrFunctionNCP compute_F; /**< pointer to the function used to compute \f$ F(z) \f$ */
  ptrFunctionJacNCP compute_nabla_F; /**< pointer to the function used to
                                        compute \f$ \nabla_z F(z) \f$ */
  NumericsMatrix* nabla_F;           /**< storage for \f$ \nabla_z F \f$*/
  void* env; /**< environment for the compute_Fmcp and compute_nabla_F function.
               When called from Python, it contains an object with compute_F and
               compute_nabla_F as methods.
               When called from C, it can reference a data struct containing
               variables needed for the computations.*/
};

#if defined(__cplusplus)
extern "C" {
#endif
/** free an NCP problem
 *
 *  \param ncp structure to free
 */
void NCP_free(NonlinearComplementarityProblem* ncp);

/** create an empty NCP problem
 *
 *  \return an MixedComplementarityProblem instance
 */
NonlinearComplementarityProblem* newNCP(void);

#if defined(__cplusplus)
}
#endif

#endif
