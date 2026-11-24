/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file
  Combined projected Moreau-Jean time-stepping scheme (Acary, CMAME 2013, Algorithm 2)
*/
#ifndef TIMESTEPPINGCOMBINEDPROJECTION_H
#define TIMESTEPPINGCOMBINEDPROJECTION_H

#include "ProjectionOnConstraintsTools.hpp"
#include "TimeStepping.hpp"
#include "Topology.hpp"

namespace siconos::simulation {

/** Combined projected Moreau-Jean time-stepping scheme, Algorithm 2 of

    V. Acary, Projected event-capturing time-stepping schemes for nonsmooth mechanical
    systems with unilateral contact and Coulomb's friction, Computer Methods in Applied
    Mechanics and Engineering 256 (2013) 224-250.

    Unlike Algorithm 1 (TimeSteppingDirectProjection), the active set is not guessed by a
    forecast: it is built within the step, by enlarging it with the constraints violated at
    the end of the step and restarting the step. The velocity problem and the position
    projection are then done on the same active set, the projection being an equality
    projection (the active contacts are held closed). One time step from
    \f$ (q_k, v_k) \f$, with \f$ \mathcal{A}_0 = \emptyset \f$ (indexSet 1 is cleared
    at the beginning of advanceToEvent()):

    for \f$ i = 0, 1, \ldots \f$ (index set iterations, at most _kIndexSetMax):
    1. Velocity step on \f$ \mathcal{A}_i \f$ (TimeStepping::newtonSolve, OSNS
       SICONOS_OSNSP_TS_VELOCITY on indexSet 1): Moreau-Jean \f$ \theta \f$-method
       with the Newton impact law
       \f$ 0 \leq G_\alpha(q_k)(v_{k+1} + e v_k) \perp \lambda_\alpha \geq 0,\
       \alpha \in \mathcal{A}_i \f$; then
       \f$ a = q_k + h(\theta v_{k+1} + (1-\theta) v_k) \f$.
       Each iteration restarts the step: for \f$ i \geq 1 \f$, updateAndSwapAllOutput()
       resets the DS states to \f$ (q_k, v_k) \f$ and recomputes the outputs and
       Jacobians there, without overwriting the memories \f$ y_k \f$ used by the
       Newton impact law.
    2. Projection on \f$ \mathcal{A}_i \f$ (indexSet 2 = indexSet 1, OSNS
       SICONOS_OSNSP_TS_POS, MLCPProjectOnConstraints with
       _doCombinedProjOnEquality = true: equality rows): while
       \f$ |g_\alpha(q_j)| \geq \f$ _constraintTolUnilateral for an active contact, and at
       most _projectionMaxIteration times, solve
       \f$ G_j G_j^\top \lambda = -g(q_j) \f$ and update the positions with the rule
       chosen by setProjectionUpdate() (ProjectionUpdate::incremental, default, or
       ProjectionUpdate::from_origin, the exact projection of a). With from_origin, the
       loop also continues while \f$ \|q_{j+1} - q_j\|_\infty > \f$ _projectionTolq.
    3. Activation (updateIndexSet(1), MoreauJeanCombinedProjectionOSI::
       addInteractionInIndexSet): \f$ \mathcal{A}_{i+1} = \mathcal{A}_i \cup
       \{\alpha : g_\alpha(q_{k+1}) < -\f$ _constraintTolUnilateral \f$\} \f$.
       Only violated constraints are added: a closed contact (\f$ g = 0 \f$, e.g. after
       the projection) is not, otherwise the equality projection would hold it closed.
       If no constraint is added, the index sets are stable and the step is over.

    The velocities are not modified by the projection. Only
    MoreauJeanCombinedProjectionOSI and indexSet level 2 for the projection are supported.

    See docs/devel/PROJECTION_SCHEMES.md for the derivation, the variants, the
    differences with the reference implementation (Scilab Hots) and the validation.
 */
class TimeSteppingCombinedProjection : public TimeStepping {
 protected:
  ACCEPT_SERIALIZATION(TimeSteppingCombinedProjection);

  /** level of IndexSet on which we project (only 2 is implemented). indexSet 2 is
   *  rebuilt from indexSet 1 before each projection
   *  (MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet(inter, 2) returns true:
   *  all the active contacts are projected, whatever the sign of their impulse)
   */
  unsigned int _indexSetLevelForProjection{2};

  /** Cumulated Number of steps performed is the Newton Loop */
  unsigned int _cumulatedNewtonNbIterations{0};

  /** Number of iteration of projection
   */
  unsigned int _nbProjectionIteration{0};

  /** Number of cumulated iteration of projection
   */
  unsigned int _nbCumulatedProjectionIteration{0};

  /** Number of iteration for stabilizating indexsets
   */
  unsigned int _nbIndexSetsIteration{0};

  /** tolerance for the violation of the equality
   *  constraints at the  position level.
   */
  double _constraintTol{1.e-8};

  /** tolerance for the violation of the unilateral
   *  constraints at the  position level. Used both as the stopping test of the
   *  projection loop (\f$ |g_\alpha| < \f$ tol on the active contacts) and as the
   *  activation threshold of indexSet 1 (\f$ g_\alpha < -\f$ tol).
   */
  double _constraintTolUnilateral{1.e-8};

  /** maximum violation for the violation of the unilateral
   *  constraints at the  position level.
   */
  double _maxViolationUnilateral{0.};

  /** maximum violation for the violation of the equality
   *  constraints at the  position level.
   */
  double _maxViolationEquality{0.};

  /** Default maximum number of projection iteration*/
  unsigned int _projectionMaxIteration{50};

  /** Default maximum number of index set activation iteration (an exception is thrown
   *  beyond)*/
  unsigned int _kIndexSetMax{50};

  /** disabled or enabled projection (Debug Projection) */
  bool _doCombinedProj{true};

  /** disabled or enabled projection On Equality (or Unilateral) for unilateral
   * constraints */
  bool _doCombinedProjOnEquality{true};

  /** update rule of the positions in the projection loop (see ProjectionUpdate) */
  ProjectionUpdate _projectionUpdate{ProjectionUpdate::incremental};

  /** with ProjectionUpdate::from_origin, the projection loop also continues while
      \f$ \|q_{j+1} - q_j\|_\infty \f$ is larger than this tolerance */
  double _projectionTolq{1.e-12};

  /** origin a of the projection (positions after the velocity step), saved for
      ProjectionUpdate::from_origin */
  ProjectionOrigin _projectionOrigin;

  /** Boolean to check if the index sets are stabilized in the Combined
   * Projection Algorithm */
  bool _isIndexSetsStable{false};

  /** update indexSets[i] of the topology, using current y and lambda values of Interactions
   *
   *  \param level unsigned int: the level of the set to be updated
   */
  void updateIndexSet(siconos::simulation::Topology::size_type level) override;

  /** compute the outputs at the beginning of the step and save them in memory (y_k).
   *
   *  newtonSolve() is called once per index set iteration, and each iteration restarts the
   *  step with the enlarged index set (Acary, CMAME 2013, Algorithm 2). The first call is
   *  done at the beginning of the step. For the next ones, the DS states are reset to their
   *  beginning-of-step values (memory) and the outputs and Jacobians are recomputed there,
   *  without saving them in memory: otherwise y_k would be overwritten by the end-of-step
   *  values of the previous iteration, and the velocity problem would be built with the
   *  Jacobians at the end of the step.
   */
  void updateAndSwapAllOutput() override;

  // struct _SimulationEffectOnOSNSP;

 public:
  void initializeOneStepNSProblem() override;

  /** Constructor with the time-discretisation.
   *
   *  \param nsds the nsds that we want to simulate
   *  \param td a pointer to a timeDiscretisation (linked to the model
   *  that owns this simulation)
   *  \param osi a one step integrator
   *  \param osnspb_velo a one step non smooth problem for the velocity
   *  formulation \param osnspb_pos a one step non smooth problem for the
   *  position formulation \param _level
   */
  TimeSteppingCombinedProjection(
      std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
      std::shared_ptr<TimeDiscretisation> td,
      std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
      std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_velo,
      std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_pos,
      unsigned int _level = 2);

  // /** default constructor
  //  */
  // TimeSteppingCombinedProjection() = default;

  virtual ~TimeSteppingCombinedProjection() noexcept = default;

  /** get the Number of iteration of projection
   *
   *  \return unsigned int nbProjectionIteration
   */
  inline unsigned int nbProjectionIteration() { return _nbProjectionIteration; }
  /** get the Number of cumulated iteration of projection
   *
   *  \return unsigned int
   */
  inline unsigned int nbCumulatedProjectionIteration() {
    return _nbCumulatedProjectionIteration;
  }

  /** get the  Cumulated Number of steps performed in the Newton Loop
   *
   *  \return unsigned int
   */
  inline unsigned int cumulatedNewtonNbIterations() { return _cumulatedNewtonNbIterations; }

  /** get the Number of iteration for stabilizating indexsets
   *
   *  \return unsigned int
   */
  inline unsigned int nbIndexSetsIteration() { return _nbIndexSetsIteration; }

  inline void setConstraintTol(double v) { _constraintTol = v; }

  inline void setConstraintTolUnilateral(double v) { _constraintTolUnilateral = v; }
  inline double constraintTolUnilateral() const { return _constraintTolUnilateral; }

  inline double maxViolationUnilateral() { return _maxViolationUnilateral; }
  inline double maxViolationEquality() { return _maxViolationEquality; }

  inline void setProjectionMaxIteration(unsigned int v) { _projectionMaxIteration = v; }

  inline void setDoCombinedProj(unsigned int v) { _doCombinedProj = v; }

  /** set the update rule of the positions in the projection loop
   *  \param v ProjectionUpdate::incremental (default) or ProjectionUpdate::from_origin
   *  (exact projection, Lagrangian systems only)
   */
  inline void setProjectionUpdate(ProjectionUpdate v) { _projectionUpdate = v; }
  inline ProjectionUpdate projectionUpdate() const { return _projectionUpdate; }

  /** set the tolerance on \f$ \|q_{j+1} - q_j\|_\infty \f$ of the from_origin fixed
   *  point (default 1e-12) */
  inline void setProjectionTolq(double v) { _projectionTolq = v; }
  inline double projectionTolq() const { return _projectionTolq; }

  inline bool doCombinedProjOnEquality() { return _doCombinedProjOnEquality; }

  void advanceToEvent() override;

  void computeCriteria(bool* runningProjection);
};
}  // namespace siconos::simulation
#endif  // TIMESTEPPINGCOMBINEDPROJECTION_H
