// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef convex_set_H
#define convex_set_H

/*!\file convex_set.h
  \brief Stacks of closed convex sets: projection, Jacobian of the projection, natural map.

  Most problems of numerics are complementarity problems over a product of closed convex
  cones or sets \f$ K = K_0 \times \cdots \times K_{m-1} \f$ in \f$ \mathbb{R}^{m \cdot
  dim} \f$, written with the natural map
  \f[
     F(r, v) = r - P_K(r - v) = 0
  \f]
  (\f$ v \f$ already contains the problem-specific terms: De Saxcé shift
  \f$ \mu \|u_T\| e_N \f$, dilatancy, parameter \f$ \rho \f$, change of variables...).
  The problems differ only by the set \f$ K \f$; a ::ConvexSet describes a stack of
  \p nblocks sets of the same type and dimension, with per-block parameters:

  | type                          | set of one block                           | p1     | p2 |
  |-------------------------------|--------------------------------------------|--------|------|
  | ::CONVEX_SET_FREE             | \f$ \mathbb{R}^{dim} \f$ (equalities)      | -      | - |
  | ::CONVEX_SET_ORTHANT          | \f$ \mathbb{R}^{dim}_+ \f$ (LCP)           | -      | - |
  | ::CONVEX_SET_BOX              | \f$ [p_1, p_2] \f$ (Relay)                 | lower  |
  upper| | ::CONVEX_SET_BALL             | \f$ B(0, p_1) \f$ (Von Mises plane stress) | radius
  | -    | | ::CONVEX_SET_CYLINDER         | \f$ \mathbb{R} \times B(0, p_1) \f$ (Von Mises
  plane strain) | radius | - | | ::CONVEX_SET_POSITIVE_CYLINDER| \f$ \mathbb{R}_+ \times B(0,
  p_1) \f$ (Tresca) | radius | -  | | ::CONVEX_SET_SOC              | \f$ \{ \|x_T\| \leq p_1
  x_N \} \f$ (Coulomb, Drucker-Prager) | mu | - | | ::CONVEX_SET_ROLLING_CONE     | rolling
  cone, dim 3 or 5                   | mu     | mu_r |

  p1 and p2 are indexed by block (p1[i]), except for ::CONVEX_SET_BOX where they are
  indexed by component (p1[i * dim + k]).

  The scaled error built on the natural map is in natural_map_error.h; the projections
  themselves in projection_orthant.h, projection_ball.h, projection_soc.h,
  projection_rolling_cone.h.
*/

#include "projection_ball.h"
#include "projection_orthant.h"
#include "projection_rolling_cone.h"
#include "projection_soc.h"

/** Type of the sets of a ::ConvexSet stack */
typedef enum {
  /** whole space: no constraint (equality rows) */
  CONVEX_SET_FREE = 0,
  /** nonnegative orthant */
  CONVEX_SET_ORTHANT = 1,
  /** box [p1, p2] (bounds per component) */
  CONVEX_SET_BOX = 2,
  /** ball of radius p1 */
  CONVEX_SET_BALL = 3,
  /** cylinder R x B(0, p1) */
  CONVEX_SET_CYLINDER = 4,
  /** positive cylinder R+ x B(0, p1) */
  CONVEX_SET_POSITIVE_CYLINDER = 5,
  /** second-order cone of coefficient p1 */
  CONVEX_SET_SOC = 6,
  /** rolling friction cone of coefficients p1 (sliding), p2 (rolling); dim 3 or 5 */
  CONVEX_SET_ROLLING_CONE = 7
} ConvexSetType;

/** A stack of nblocks convex sets of the same type and dimension (see convex_set.h) */
typedef struct {
  /** type of the sets */
  ConvexSetType type;
  /** dimension of one block */
  int dim;
  /** number of blocks */
  int nblocks;
  /** first parameter, per block (per component for a box); may be NULL if unused */
  const double* p1;
  /** second parameter, per block (per component for a box); may be NULL if unused */
  const double* p2;
} ConvexSet;

/** Build a ::ConvexSet.

    \param type type of the sets
    \param dim dimension of one block
    \param nblocks number of blocks
    \param p1 first parameter (see convex_set.h), may be NULL if unused
    \param p2 second parameter, may be NULL if unused
    \return the stack descriptor (it does not own p1, p2)
*/
static inline ConvexSet convex_set(ConvexSetType type, int dim, int nblocks, const double* p1,
                                   const double* p2) {
  ConvexSet K = {type, dim, nblocks, p1, p2};
  return K;
}

/** Projection of one block on its set, in place.

    \param K the stack of sets
    \param i block index
    \param[in,out] x the block (size K->dim)
    \return a status, whose meaning depends on the type (see the projection_*.h headers)
*/
static inline unsigned convex_set_project_block(const ConvexSet* K, int i, double* x) {
  switch (K->type) {
    case CONVEX_SET_SOC:
      return (K->dim == 3) ? soc_projection_3d(x, K->p1[i])
                           : soc_projection(x, K->p1[i], K->dim);
    case CONVEX_SET_ORTHANT:
      orthant_projection(x, K->dim);
      return 0;
    case CONVEX_SET_BOX:
      box_projection(x, K->dim, K->p1 + i * K->dim, K->p2 + i * K->dim);
      return 0;
    case CONVEX_SET_BALL:
      return ball_projection(x, K->dim, K->p1[i]);
    case CONVEX_SET_CYLINDER:
      return cylinder_projection(x, K->dim, K->p1[i]);
    case CONVEX_SET_POSITIVE_CYLINDER:
      return positive_cylinder_projection(x, K->dim, K->p1[i]);
    case CONVEX_SET_ROLLING_CONE:
      return (K->dim == 5) ? rolling_cone_projection(x, K->p1[i], K->p2[i])
                           : rolling_cone_2d_projection(x, K->p1[i], K->p2[i]);
    case CONVEX_SET_FREE:
    default:
      return 0;
  }
}

#if defined(__cplusplus)
extern "C" {
#endif

/** Natural map of one block, \f$ F_i = r_i - P_{K_i}(r_i - v_i) \f$, and its squared norm.

    \param K the stack of sets
    \param i block index
    \param r the block of r (size K->dim)
    \param v the block of v (size K->dim)
    \param[out] F the natural map of the block (size K->dim), may be NULL
    \return \f$ \|F_i\|^2 \f$
*/
double convex_set_natural_map_block(const ConvexSet* K, int i, const double* r,
                                    const double* v, double* F);

/** Projection on the stack of sets, in place.

    \param K the stack of sets
    \param[in,out] x vector of size K->dim * K->nblocks
*/
void convex_set_project(const ConvexSet* K, double* x);

/** An element of the generalized Jacobian of the projection of one block, dense
    K->dim x K->dim, column-major.

    \param K the stack of sets
    \param i block index
    \param[out] J the Jacobian
    \param x the block (not modified)
    \return 0 on success, 1 if no Jacobian is implemented for this type and dimension
*/
int convex_set_projection_jacobian_block(const ConvexSet* K, int i, double* J,
                                         const double* x);

/** Natural map on the stack: \f$ F = r - P_K(r - v) \f$.

    \param K the stack of sets
    \param r vector of size K->dim * K->nblocks
    \param v vector of the same size
    \param[out] F the natural map, same size (may alias neither r nor v)
    \return \f$ \|F\|^2 \f$
*/
double convex_set_natural_map(const ConvexSet* K, const double* r, const double* v, double* F);

#if defined(__cplusplus)
}
#endif

#endif
