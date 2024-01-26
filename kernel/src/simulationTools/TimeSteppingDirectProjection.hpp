/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*! \file
  Time-Stepping simulation with projections on constraints
*/
#ifndef TIMESTEPPINGDIRECTPROJECTION_H
#define TIMESTEPPINGDIRECTPROJECTION_H

#include "TimeStepping.hpp"

namespace siconos::simulation {
/**
   Time-Stepping scheme with a direct projection onto the constraint
   thanks to the GGL augmentation of the system

   For details, have a look on
   Projected event-capturing time-stepping schemes for nonsmooth mechanical
   systems with unilateral contact and Coulomb's friction Vincent Acary Computer
   Methods in Applied Mechanics and Engineering, Elsevier, 2013, 256, pp.
   224-250
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

  void initOSNS() override;

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

  inline void setDoProj(bool v) { _doProj = v; }
  inline void setDoOnlyProj(bool v) { _doOnlyProj = v; }

  void advanceToEvent() override;

  void nextStep() override;

  void computeCriteria(bool *runningProjection);

  void newtonSolve(double criterion, unsigned int maxStep) override;
};
}  // namespace siconos::simulation

#endif  // TIMESTEPPINGDIRECTPROJECTION_H
