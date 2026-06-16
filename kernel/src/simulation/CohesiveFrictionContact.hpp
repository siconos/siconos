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
  Cohesive Friction-Contact Non-Smooth Problem
*/
#ifndef COHESIVEFRICTIONCONTACT_H
#define COHESIVEFRICTIONCONTACT_H

#include "FrictionContact.hpp"

struct FrictionContactProblem;
struct SolverOptions;

namespace siconos::nonsmooth_formulations {

/**
   Formalization and Resolution of a Cohesive Friction-Contact Problem

   This class extends FrictionContact to handle cohesive zone models.
   It adds support for cohesive forces that act before contact occurs,
   modeling the progressive damage and failure of interfaces.

   The cohesive contribution is added to the q vector of the LCP/LCP-like
   problem, shifting the contact condition to account for the cohesion force.

   With:
   - \f$ velocity = q + M reaction \f$
   - Cohesion modifies q: \f$ q_{eff} = q + V q_{cohesion} \f$
   - Standard friction law applies

   The dimension of the problem (2D or 3D) is given by the variable
   contactProblemDim and the proper Numerics driver will be called according to
   this value.

   \b Construction: just set Numerics Solver id

   Main functions:
   - compute(time) formalize, solve and post-process the problem.

   pre- and post-pro are common to all LinearOSNS and defined in this class.
 */
class CohesiveFrictionContact : public FrictionContact {
 protected:
  ACCEPT_SERIALIZATION(CohesiveFrictionContact);

  /** Matrix V for cohesive contribution mapping
   *  Maps cohesive forces from indexSet0 to the OSNS problem
   */
  std::shared_ptr<OSNSMatrix> _V{nullptr};

  /** Matrix H0 for direct assembly */
  std::shared_ptr<OSNSMatrix> _H0{nullptr};

  /** Cohesion force contribution to q vector */
  std::shared_ptr<siconos::algebra::SiconosVector> _q_cohesion{nullptr};

  /** Size of the cohesive problem (number of contacts in indexSet0) */
  siconos::algebra::Index _sizeOutput_cohesion{0};

 public:
  /** constructor (solver id and dimension)
   *
   *  \param dimPb dimension (2D or 3D) of the friction-contact problem (default: 3)
   *  \param numericsSolverId id of the solver to be used (default: SICONOS_FRICTION_3D_NSGS)
   */
  CohesiveFrictionContact(int dimPb = 3, int numericsSolverId = SICONOS_FRICTION_3D_NSGS);

  /** constructor from a pre-defined solver options set
   *
   *  \param dimPb dimension (2D or 3D) of the friction-contact problem
   *  \param options the options set
   */
  CohesiveFrictionContact(int dimPb, std::shared_ptr<SolverOptions> options);

  /** destructor
   */
  ~CohesiveFrictionContact() noexcept override = default;

  /** initialize the CohesiveFrictionContact problem
   *  \param simulation the simulation, owner of this OSNS problem
   */
  void initialize(std::shared_ptr<siconos::simulation::Simulation> simulation) override;

  /** Compute cohesion contribution for a single interaction
   *  \param vertex_inter vertex descriptor for the interaction
   *  \param pos position in the global vector
   */
  void computeQCohesionBlock(
      siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
      siconos::algebra::Index pos);

  /** Update q vector with cohesion contribution
   *  \param time the current time
   */
  void updateQWithQCohesion(double time);

  /** Compute the V matrix for cohesive contribution
   */
  void computeV();

  /** build problem coefficients (including cohesion)
   *  \param time the current time
   *  \return true if succeeded
   */
  bool preCompute(double time) override;

  /** Post-process after solve
   */
  void postCompute() override;

  /** Check compatibility of NS law with this OSNS problem
   *  \param nslaw the non-smooth law to check
   *  \return true if compatible
   */
  bool checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw) override;

  /** print the data to the screen */
  void display() const override;
};

}  // namespace siconos::nonsmooth_formulations

#endif  // COHESIVEFRICTIONCONTACT_H
