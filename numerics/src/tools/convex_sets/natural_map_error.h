// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef natural_map_error_H
#define natural_map_error_H

/*!\file natural_map_error.h
  \brief Scaled error of a complementarity problem over a product of convex sets,
  measured with the natural map.

  For a stack of sets \f$ K \f$ (convex_set.h), \f$ r \f$ the primal variable
  (reaction, stress...) and \f$ v \f$ the dual variable (modified velocity, strain
  rate..., already multiplied by \f$ \rho \f$), the residual is the natural map
  \f$ F = r - P_K(r - v) \f$ and the error is
  - no scaling: \f$ \|F\| \f$;
  - per block (::NATURAL_MAP_ERROR_SCOPE_PER_BLOCK):
    \f$ \sqrt{\sum_i \|F_i\|^2 / \sigma(r_i, v_i)^2} \f$;
  - per vector (::NATURAL_MAP_ERROR_SCOPE_PER_VECTOR): \f$ \|F\| / \sigma(r, v) \f$;

  where the scale \f$ \sigma \f$ combines a norm (::NATURAL_MAP_ERROR_NORM) of r and of v
  with a reduction (::NATURAL_MAP_ERROR_REDUCTION):
  - MAX: \f$ \min(\max(1, \|r\|), \max(1, \|v\|)) \f$, i.e. the error is the max of the
    two relative errors \f$ \|F\|/\max(1,\|r\|) \f$ and \f$ \|F\|/\max(1,\|v\|) \f$;
  - MEAN: \f$ \sqrt{\max(1, \|r\|) \max(1, \|v\|)} \f$;
  - MIXED: \f$ \max(1, \|r\| + \|v\|) \f$.

  The norms of a block \f$ x = (x_N, x_T) \f$ are
  - EUCLIDEAN: \f$ \|x\|_2 \f$;
  - CONE (spectral norm of the Jordan algebra of the second-order cone):
    \f$ |x_N| + \|x_T\| \f$;
  - CONE_MU: \f$ |x_N| + \mu \|x_T\| \f$, \f$ \mu \f$ the coefficient of a cone block
    (1 for the other sets);
  - CONE_WEIGHTED: \f$ \sqrt{x_N^2 + \alpha \|x_T\|^2} \f$;

  and over a stack: the sum of the block norms (CONE, CONE_MU), the square root of the sum
  of the squares (EUCLIDEAN, CONE_WEIGHTED). The theory (scale invariance, behaviour on
  the three regions of the cone) is in docs/devel/Error_SOC_natural_map.tex; the use in
  the solvers in docs/devel/NATURAL_MAP_ERROR.md.

  The same option indices are used by all the problems (FrictionContact 2D/3D, global
  FrictionContact, Plasticity, MixedProblem): iparam[::NATURAL_MAP_ERROR_IPARAM_NORM],
  ..., dparam[::NATURAL_MAP_ERROR_DPARAM_ALPHA].
*/

#include "NumericsFwd.h"  // for SolverOptions
#include "convex_set.h"   // for ConvexSet

/** Norm used in the scale of the error */
enum NATURAL_MAP_ERROR_NORM {
  /** Euclidean norm (default: value 0 of zero-initialized options) */
  NATURAL_MAP_ERROR_NORM_EUCLIDEAN = 0,
  /** no scaling: absolute error \f$ \|F\| \f$ */
  NATURAL_MAP_ERROR_NORM_NO = 1,
  /** \f$ |x_N| + \|x_T\| \f$ */
  NATURAL_MAP_ERROR_NORM_CONE = 3,
  /** \f$ |x_N| + \mu \|x_T\| \f$ */
  NATURAL_MAP_ERROR_NORM_CONE_MU = 5,
  /** \f$ \sqrt{x_N^2 + \alpha \|x_T\|^2} \f$ */
  NATURAL_MAP_ERROR_NORM_CONE_WEIGHTED = 6
};

/** Scope of the scale */
enum NATURAL_MAP_ERROR_SCOPE {
  /** one scale per block (contact, element) */
  NATURAL_MAP_ERROR_SCOPE_PER_BLOCK = 0,
  /** one scale for the whole vector */
  NATURAL_MAP_ERROR_SCOPE_PER_VECTOR = 1
};

/** Reduction of the norms of r and v into the scale */
enum NATURAL_MAP_ERROR_REDUCTION {
  /** \f$ \min(\max(1, \|r\|), \max(1, \|v\|)) \f$ */
  NATURAL_MAP_ERROR_REDUCTION_MAX = 0,
  /** \f$ \sqrt{\max(1, \|r\|) \max(1, \|v\|)} \f$ */
  NATURAL_MAP_ERROR_REDUCTION_MEAN = 1,
  /** \f$ \max(1, \|r\| + \|v\|) \f$ */
  NATURAL_MAP_ERROR_REDUCTION_MIXED = 2
};

/** Indices of the scaling options in SolverOptions, common to all the problems */
enum NATURAL_MAP_ERROR_OPTIONS {
  /** iparam: norm (::NATURAL_MAP_ERROR_NORM) */
  NATURAL_MAP_ERROR_IPARAM_NORM = 24,
  /** iparam: scope (::NATURAL_MAP_ERROR_SCOPE) */
  NATURAL_MAP_ERROR_IPARAM_SCOPE = 25,
  /** iparam: reduction (::NATURAL_MAP_ERROR_REDUCTION) */
  NATURAL_MAP_ERROR_IPARAM_REDUCTION = 26,
  /** dparam: weight alpha of ::NATURAL_MAP_ERROR_NORM_CONE_WEIGHTED */
  NATURAL_MAP_ERROR_DPARAM_ALPHA = 5
};

/** Scaling of the error */
typedef struct {
  /** norm (::NATURAL_MAP_ERROR_NORM) */
  int norm;
  /** scope (::NATURAL_MAP_ERROR_SCOPE) */
  int scope;
  /** reduction (::NATURAL_MAP_ERROR_REDUCTION) */
  int reduction;
  /** weight of ::NATURAL_MAP_ERROR_NORM_CONE_WEIGHTED */
  double alpha;
} NaturalMapErrorScaling;

/** Partial sums of an error over several stacks of sets (product of sets of different
    types), see natural_map_error_add() */
typedef struct {
  /** sum of the squared residuals (divided by the squared block scales per block) */
  double F2;
  /** accumulated norm of r (sum, or sum of squares, depending on the norm) */
  double nr;
  /** accumulated norm of v */
  double nv;
} NaturalMapErrorSum;

#if defined(__cplusplus)
extern "C" {
#endif

/** Scaling read from the options (indices ::NATURAL_MAP_ERROR_OPTIONS).

    \param options solver options, may be NULL
    \param defaults scaling returned when options is NULL
    \return the scaling
*/
NaturalMapErrorScaling natural_map_error_scaling(const SolverOptions* options,
                                                 NaturalMapErrorScaling defaults);

/** Norm of a stack of m blocks of dimension dim (see natural_map_error.h).

    \param norm the norm (::NATURAL_MAP_ERROR_NORM; NO gives 1)
    \param x vector of size m * dim
    \param m number of blocks
    \param dim dimension of a block
    \param mu coefficients of the blocks for ::NATURAL_MAP_ERROR_NORM_CONE_MU (may be NULL:
           1)
    \param alpha weight of ::NATURAL_MAP_ERROR_NORM_CONE_WEIGHTED
    \return the norm
*/
double natural_map_error_norm(int norm, const double* x, int m, int dim, const double* mu,
                              double alpha);

/** Scale \f$ \sigma(x, y) \geq 1 \f$ combining the norms of two stacks (see
    natural_map_error.h).

    \param norm the norm (::NATURAL_MAP_ERROR_NORM; NO gives 1)
    \param reduction the reduction (::NATURAL_MAP_ERROR_REDUCTION)
    \param x first vector, size m * dim
    \param y second vector, size m * dim
    \param m number of blocks
    \param dim dimension of a block
    \param mu coefficients of the blocks for CONE_MU (may be NULL: 1)
    \param alpha weight of CONE_WEIGHTED
    \return the safeguarded scale
*/
double natural_map_error_scale(int norm, int reduction, const double* x, const double* y,
                               int m, int dim, const double* mu, double alpha);

/** Start an error sum.

    \param[out] sum the sum to initialize
*/
void natural_map_error_init(NaturalMapErrorSum* sum);

/** Add the residuals of a stack of sets to an error sum.

    \param s the scaling
    \param K the stack of sets
    \param r primal variable, size K->dim * K->nblocks
    \param v dual variable (modified, multiplied by rho), same size
    \param[in,out] sum the error sum
*/
void natural_map_error_add(const NaturalMapErrorScaling* s, const ConvexSet* K,
                           const double* r, const double* v, NaturalMapErrorSum* sum);

/** Error from an error sum.

    \param s the scaling
    \param sum the error sum
    \return the error
*/
double natural_map_error_value(const NaturalMapErrorScaling* s, const NaturalMapErrorSum* sum);

/** Scaled error of one stack of sets: natural_map_error_init(), natural_map_error_add(),
    natural_map_error_value().

    \param s the scaling
    \param K the stack of sets
    \param r primal variable, size K->dim * K->nblocks
    \param v dual variable (modified, multiplied by rho), same size
    \return the error
*/
double natural_map_error(const NaturalMapErrorScaling* s, const ConvexSet* K, const double* r,
                         const double* v);

#if defined(__cplusplus)
}
#endif

#endif
