/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MCP_PROBLEM_H
#define MCP_PROBLEM_H

/*!\file MixedComplementarityProblem.h
 */

#include "NumericsFwd.h"  // for MixedComplementarityProblem, NumericsMatrix

/**
    type for user defined function used to compute Fmcp and its jacobian.
    TODO : set properly the list of arguments for this function, when
    things will be clearer ...
 */
typedef void (*ptrFunctionMCP)(int size, double* z, double* F);
typedef void (*ptrFunctionMCP2)(void* env, int n, double* z, double* F);
typedef void (*ptrFunctionMCP_nabla)(void* env, int n, double* z, NumericsMatrix* F);

/**
    Structure that contains and defines a MixedComplementarityProblem
*/
struct MixedComplementarityProblem {
  int n1;                                  /**< number of equalities constraints */
  int n2;                                  /**< size of complementarity variables */
  ptrFunctionMCP2 compute_Fmcp;            /**< pointer to the function used to compute
                                              \f$ F_{mcp}(z) = (G(z), H(z)) \f$ */
  ptrFunctionMCP_nabla compute_nabla_Fmcp; /**< pointer to the function used to
                                              compute \f$ \nabla_z F_{mcp} \f$ */
  NumericsMatrix* nabla_Fmcp;              /**< storage for \f$ \nabla_z F_{mcp} \f$*/
  void* env; /**< environment for the compute_Fmcp and compute_nabla_Fmcp
               function. When called from Python, it contains an object with
               compute_Fmcp and compute_nabla_Fmcp as methods.
               When called from C, it can reference a data struct containing
               variables needed for the computations.*/
};

/**
    The structure that defines a Mixed Complementarity problem (MCP) : find two
    vectors \f$ (z,w \in {{\mathrm{I\!R}}}^{n+m}) \f$ such that:

    \f[
    w &= \begin{pmatrix}w_e\\w_i\end{pmatrix} = F(z) \ \
    w_e &=0 \					       \
    0 &\le w_i \perp z_i \ge 0
    \f]

    where "i" (resp. "e") stands for inequalities (resp. equalities). The vector
    \f$ z \f$ is splitted like \f$ w \f$:

    \f[
    z =\begin{pmatrix}z_e\\z_i\end{pmatrix}
    \f]

    \f$ z_i,w_i \f$ are vectors of size sizeEqualities, \f$ z_e,w_e \f$ vectors
    of size sizeInequalities

    and \f$ F \f$ is a non linear function that must be user-defined.
*/
struct MixedComplementarityProblem_old {
  int sizeEqualities;              /**< size of equalities \f$ z_e, w_e \f$ size */
  int sizeInequalities;            /**< size of inequalities \f$ z_i,w_i \f$ size */
  ptrFunctionMCP computeFmcp;      /**< pointer to the function to compute F(z) */
  ptrFunctionMCP computeNablaFmcp; /** pointer to the function to compute the
                                      jacobian of F(z) */
  /** The value F(z) */
  double* Fmcp;
  /** jacobian of F(z) */
  double* nablaFmcp;
};

#if defined(__cplusplus)
extern "C" {
#endif

/** free a MCP problem
 *
 *  \param mcp structure to free
 */
void mixedComplementarityProblem_free(MixedComplementarityProblem* mcp);

/** \return an empty MixedComplementarityProblem instance
 */
MixedComplementarityProblem* mixedComplementarityProblem_new(void);

/** Delete a MixedComplementarityProblem_old
 *
 *  \param problem pointer to a MixedComplementarityProblem_old to delete
 */
void mixedComplementarityProblem_old_free(MixedComplementarityProblem_old* problem);

#if defined(__cplusplus)
}
#endif

#endif
