// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file ProjectionOnConstraintsTools.hpp
  \brief Position update of the projection on constraints, shared by
  TimeSteppingDirectProjection and TimeSteppingCombinedProjection.

  See docs/devel/PROJECTION_SCHEMES.md for the schemes, the variants and their validation.
*/
#ifndef PROJECTIONONCONSTRAINTSTOOLS_H
#define PROJECTIONONCONSTRAINTSTOOLS_H

#include <unordered_map>

#include "SiconosVector.hpp"
#include "SimulationGraphs.hpp"

namespace siconos::modeling {
class NonSmoothDynamicalSystem;
class DynamicalSystem;
}  // namespace siconos::modeling

namespace siconos::simulation {

/** How the positions are updated in the projection loop of the projected Moreau-Jean
    schemes (TimeSteppingDirectProjection, Algorithm 1 of Acary, CMAME 2013, and
    TimeSteppingCombinedProjection, Algorithm 2).

    After the velocity step, \f$ a = q_k + h(\theta v_{k+1} + (1-\theta) v_k) \f$ may
    violate the constraints \f$ g(q) \geq 0 \f$. The projection is

    \f[ \min_q \tfrac12 \|q - a\|^2 \quad \text{s.t.} \quad g_\alpha(q) = 0 \ (\text{Alg. 2,
    active set}) \quad \text{or} \quad g(q) \geq 0 \ (\text{Alg. 1}), \f]

    whose first-order (KKT) conditions are
    \f$ q = a + G(q)^\top \mu,\ g_\alpha(q) = 0 \f$
    (resp. \f$ 0 \leq g(q) \perp \mu \geq 0 \f$), with \f$ G = \nabla_q^\top g \f$.
    It is solved by a sequence of linearized projections (one MLCP per iteration,
    MLCPProjectOnConstraints, matrix \f$ G_j G_j^\top \f$ with \f$ G_j = G(q_j) \f$,
    \f$ q_0 = a \f$). The two variants differ by the right-hand side and the update:
 */
enum class ProjectionUpdate {
  /** \f$ G_j G_j^\top \lambda_j = -g(q_j) \f$ (complementarity for Alg. 1),
      \f$ q_{j+1} = q_j + G_j^\top \lambda_j \f$.

      Gauss-Newton on \f$ g(q) = 0 \f$: it converges to a point of the constraint
      manifold, but \f$ q - a = \sum_j G(q_j)^\top \lambda_j \f$ is not normal to the
      manifold at q, so it is not the projection of a. The difference with the
      projection is \f$ O(\kappa \|\mu\|^2) \f$ (\f$ \kappa \f$: curvature of g), zero for
      linear constraints. Historical Siconos behaviour (default). */
  incremental,
  /** \f$ G_j G_j^\top \mu_j = -\left[g(q_j) + G_j (a - q_j)\right] \f$
      (complementarity for Alg. 1), \f$ q_{j+1} = a + G_j^\top \mu_j \f$.

      Fixed-point iteration on the KKT conditions \f$ q = a + G(q)^\top \mu \f$: at
      convergence, q is the (Euclidean) projection of a. This is the projection of the
      Scilab code Hots (ProjectedMoreau.sci, CombinedProjectedMoreau.sci) used for
      Acary, CMAME 2013. Only for Lagrangian systems (LagrangianDS, LagrangianSparseDS). */
  from_origin
};

/** origin a of the projection, for each dynamical system */
using ProjectionOrigin = std::unordered_map<const siconos::modeling::DynamicalSystem*,
                                            siconos::algebra::SiconosVector>;

/** Save the origin of the projection, a[ds] = q(ds), for all the dynamical systems.

    Called once per time step, before the projection loop, when the positions are the
    result of the velocity step (q = a).

    \param nsds the nonsmooth dynamical system
    \param[out] a the saved positions
    \throw if a NewtonEulerDS is present (ProjectionUpdate::from_origin is not implemented
    for quaternion-based coordinates)
 */
void saveProjectionOrigin(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          ProjectionOrigin& a);

/** Shift the right-hand side of the linearized projection for
    ProjectionUpdate::from_origin:

    \f[ y^0 \leftarrow g(q_j) + G_j (a - q_j) \f]

    for each interaction of indexSet, where \f$ (a - q_j) \f$ is the concatenation over
    the source and target dynamical systems of the interaction (same order as the
    columns of LagrangianR::contactMatrix()). MLCPProjectOnConstraints::computeqBlock uses
    \f$ y^0 \f$ as the vector q of the MLCP.

    Must be called right after computeCriteria() (which computes \f$ y^0 = g(q_j) \f$
    and \f$ H = G(q_j) \f$) and before the projection OSNS is solved. \f$ y^0 \f$ is
    recomputed by the next computeCriteria().

    \param indexSet the index set of the projection problem
    \param a the origin of the projection (saveProjectionOrigin)
 */
void shiftProjectionRhs(siconos::graphs::InteractionsGraph& indexSet,
                        const ProjectionOrigin& a);

/** Update the positions with the solution of the linearized projection,
    \f$ p^0 = G_j^\top \lambda \f$ (input of level 0, computed by updateInput(0)):

    - ProjectionUpdate::incremental: \f$ q \leftarrow q + p^0 \f$;
    - ProjectionUpdate::from_origin: \f$ q \leftarrow a + p^0 \f$.

    For a NewtonEulerDS (incremental only): \f$ q \leftarrow q + p^0 \f$, then the
    quaternion is normalized and T(q) recomputed.

    \param nsds the nonsmooth dynamical system
    \param a the origin of the projection (only used for from_origin)
    \param mode the update rule
    \return \f$ \max_{ds} \|q_{j+1} - q_j\|_\infty \f$, used as stopping test of the
    from_origin fixed point
 */
double updatePositionsAfterProjection(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                                      const ProjectionOrigin& a, ProjectionUpdate mode);

}  // namespace siconos::simulation

#endif  // PROJECTIONONCONSTRAINTSTOOLS_H
