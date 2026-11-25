/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file
  Time-Stepping simulation with projections on constraints
*/
#ifndef TIMESTEPPINGDIRECTPROJECTION_H
#define TIMESTEPPINGDIRECTPROJECTION_H

#include "ProjectionOnConstraintsTools.hpp"
#include "TimeStepping.hpp"

namespace siconos::simulation {
/**
   Projected Moreau-Jean time-stepping scheme, Algorithm 1 ("direct projection") of

   V. Acary, Projected event-capturing time-stepping schemes for nonsmooth mechanical
   systems with unilateral contact and Coulomb's friction, Computer Methods in Applied
   Mechanics and Engineering 256 (2013) 224-250.

   The Moreau-Jean scheme satisfies the constraints at the velocity level only, so the
   positions drift (penetration). Algorithm 1 adds, at the end of each step, a
   projection of the positions onto the admissible set. One time step from
   \f$ (q_k, v_k) \f$:

   1. Activation (MoreauJeanDirectProjectionOSI::addInteractionInIndexSet, indexSet 1):
      a contact \f$ \alpha \f$ is active if its forecast gap is small,
      \f$ g_\alpha(q_k) + \gamma h \dot g_\alpha(q_k, v_k) \leq \epsilon_{act} \f$
      (\f$ \gamma = 1/2 \f$ by default, MoreauJeanOSI::setGamma;
      \f$ \epsilon_{act} \f$ = MoreauJeanDirectProjectionOSI::activateYPosThreshold(),
      1e-7 by default).
   2. Velocity step (TimeStepping::newtonSolve, OSNS SICONOS_OSNSP_TS_VELOCITY, level 1):
      Moreau-Jean \f$ \theta \f$-method with the Newton impact law on the active
      contacts, \f$ 0 \leq G_\alpha (v_{k+1} + e v_k) \perp \lambda_\alpha \geq 0 \f$,
      with \f$ G = G(q_k) \f$; then
      \f$ a = q_k + h(\theta v_{k+1} + (1-\theta) v_k) \f$ (the DS positions).
   3. Projection (OSNS SICONOS_OSNSP_TS_POS, MLCPProjectOnConstraints, input/output
      level 0, on the interactions of indexSet(_indexSetLevelForProjection)): while a
      constraint is violated by more than _constraintTolUnilateral (unilateral) or
      _constraintTol (equality), and at most _projectionMaxIteration times, solve the
      linearized projection
      \f$ 0 \leq g(q_j) + G_j G_j^\top \lambda \perp \lambda \geq 0 \f$
      (equality rows for equality constraints) and update the positions. The update rule
      is selected by setProjectionUpdate() (ProjectionUpdate::incremental, default, or
      ProjectionUpdate::from_origin, the exact projection of a; see ProjectionUpdate).
      With from_origin, the loop also continues while
      \f$ \|q_{j+1} - q_j\|_\infty > \f$ _projectionTolq (convergence of the fixed point).
   4. The velocities are not modified by the projection.

   With _indexSetLevelForProjection = 0 the projection involves all the constraints
   (the setting of Algorithm 1); with 1, only the contacts active in the velocity step.

   See docs/devel/PROJECTION_SCHEMES.md for the derivation, the variants, the
   differences with the reference implementation (Scilab Hots) and the validation.
*/
class TimeSteppingDirectProjection : public TimeStepping {
 protected:
  ACCEPT_SERIALIZATION(TimeSteppingDirectProjection);

  /** level of IndexSet on which we project (default =1 (activated contact))
   */
  unsigned int _indexSetLevelForProjection{1};

  /** Number of iteration of projection
   */
  unsigned int _nbProjectionIteration{0};

  /** tolerance for the violation of the equality
   *  constraints at the  position level.
   */
  double _constraintTol{1.e-08};

  /** tolerance for the violation of the unilateral
   *  constraints at the  position level.
   */
  double _constraintTolUnilateral{1.e-08};

  /** maximum violation for the violation of the unilateral
   *  constraints at the  position level.
   */
  double _maxViolationUnilateral{0.};

  /** maximum violation for the violation of the equality
   *  constraints at the  position level.
   */
  double _maxViolationEquality{0.};

  /** Default maximum number of projection iteration*/
  unsigned int _projectionMaxIteration{10};

  /** disabled or enabled projection (Debug Projection) */
  bool _doProj{true};
  bool _doOnlyProj{false};

  /** update rule of the positions in the projection loop (see ProjectionUpdate) */
  ProjectionUpdate _projectionUpdate{ProjectionUpdate::incremental};

  /** with ProjectionUpdate::from_origin, the projection loop also continues while
      \f$ \|q_{j+1} - q_j\|_\infty \f$ is larger than this tolerance */
  double _projectionTolq{1.e-12};

  /** origin a of the projection (positions after the velocity step), saved for
      ProjectionUpdate::from_origin */
  ProjectionOrigin _projectionOrigin;

 public:
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
  TimeSteppingDirectProjection(
      std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
      std::shared_ptr<TimeDiscretisation> td,
      std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
      std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_velo,
      std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_pos,
      unsigned int _level = 1);

  void initializeOneStepNSProblem() override;

  // /** default constructor
  //  */
  // TimeSteppingDirectProjection() = default;

  virtual ~TimeSteppingDirectProjection() noexcept = default;

  /** \return the Number of iteration of projection
   */
  inline unsigned int nbProjectionIteration() { return _nbProjectionIteration; }

  inline void setConstraintTol(double v) { _constraintTol = v; }

  inline void setConstraintTolUnilateral(double v) { _constraintTolUnilateral = v; }

  inline double maxViolationUnilateral() { return _maxViolationUnilateral; }

  inline double maxViolationEquality() { return _maxViolationEquality; }

  inline void setProjectionMaxIteration(unsigned int v) { _projectionMaxIteration = v; }

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

  inline void setDoProj(bool v) { _doProj = v; }
  inline void setDoOnlyProj(bool v) { _doOnlyProj = v; }

  void advanceToEvent() override;

  void nextStep() override;

  void computeCriteria(bool* runningProjection);

  void newtonSolve(double criterion, unsigned int maxStep) override;
};
}  // namespace siconos::simulation

#endif  // TIMESTEPPINGDIRECTPROJECTION_H
