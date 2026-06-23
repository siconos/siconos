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
/*! \file CohesiveFrictionContact.hpp
  \brief Cohesive Friction-Contact Non-Smooth Problem Solver

  \section sec_cfc_overview Overview

  CohesiveFrictionContact extends the standard FrictionContact solver to handle
  cohesive zone models. While standard contact only acts when bodies overlap 
  (gap <= 0), cohesive forces act when the interface is intact and the gap is
  positive but small.

  This is essential for simulating:
  - Delamination of composite materials
  - Fracture mechanics with crack propagation
  - Adhesive contact and debonding
  - Rock mechanics with joint failure

  \section sec_cfc_mathematical Mathematical Formulation

  The standard friction-contact problem is:
  \f[
  \begin{cases}
    velocity = q + M \cdot reaction \\
    gap \leq 0 \perp reaction_n \geq 0 \\
    \|reaction_t\| \leq \mu \cdot reaction_n
  \end{cases}
  \f]

  With cohesive forces, the velocity equation becomes:
  \f[
  velocity = (q + q_{cohesion}) + M \cdot reaction
  \f]

  where \f$ q_{cohesion} \f$ represents the effect of cohesive traction forces
  on the relative velocity. This shifts the contact condition to account for
  the fact that bodies may be "pulled" together by cohesive forces even when
  not in physical contact.

  \section sec_cfc_matrixV The V Matrix

  The cohesive contribution is mapped to the OSNS problem via matrix V:
  \f[
  q_{eff} = q + V \cdot q_{cohesion}
  \f]

  where:
  - \f$ q \f$ is the original OSNS vector (local velocity at contact)
  - \f$ q_{cohesion} \f$ contains cohesive forces from all interactions
  - V maps cohesive forces to the contact space

  Matrix V is computed from the relation's Jacobian (H matrix) and accounts
  for the different contact points where cohesive and contact forces act.

  \section sec_cfc_algorithm Algorithm

  The solution process at each time step:
  1. Update contact detection (broadphase/narrowphase)
  2. For each cohesive interaction, compute cohesive force \f$ r_{cohesion} \f$
  3. Assemble \f$ q_{cohesion} \f$ vector from all cohesive interactions
  4. Compute V matrix mapping cohesive to contact space
  5. Update effective q: \f$ q_{eff} = q + V \cdot q_{cohesion} \f$
  6. Solve standard friction-contact problem with modified q
  7. Update internal cohesive state (damage parameters)

  \section sec_cfc_usage Usage Example

  \code
  #include "CohesiveFrictionContact.hpp"
  #include "BinaryCohesiveNSL.hpp"

  // Create the cohesive friction-contact solver
  auto osnsp = std::make_shared<siconos::nonsmooth_formulations::CohesiveFrictionContact>(
      3,                          // dimension: 3D
      SICONOS_FRICTION_3D_NSGS    // solver: Non-Smooth Gauss-Seidel
  );

  // Create simulation with cohesive solver
  auto simulation = std::make_shared<TimeStepping>(td, osi, osnsp);

  // Create cohesive interactions
  auto nslaw = std::make_shared<BinaryCohesiveNSL>(
      0.0, 0.0, 0.3,  // en, et, mu
      10.0e6,         // sigma_c: 10 MPa
      0.1e-3,         // delta_c: 0.1 mm
      3               // 3D
  );
  auto inter = std::make_shared<Interaction>(nslaw, relation);
  \endcode

  \see FrictionContact for the base friction-contact solver
  \see CohesiveZoneModelNIFNSL for cohesive law interface
  \see BinaryCohesiveNSL for a concrete cohesive law implementation
*/

#ifndef COHESIVEFRICTIONCONTACT_H
#define COHESIVEFRICTIONCONTACT_H

#include "FrictionContact.hpp"

struct FrictionContactProblem;
struct SolverOptions;

namespace siconos::nonsmooth_formulations {

/** \class CohesiveFrictionContact
 * \brief Solver for friction-contact problems with cohesive zone models
 *
 * This class extends FrictionContact to handle cohesive forces that act
 * before physical contact occurs. It is used with cohesive non-smooth laws
 * derived from CohesiveZoneModelNIFNSL.
 *
 * \par Key Differences from FrictionContact
 * - Additional q_cohesion vector for cohesive force contribution
 * - Matrix V for mapping cohesive forces to contact space
 * - Special handling of interactions in indexSet0 (distant but cohesive)
 *
 * \par OSNS Problem Structure
 * The standard friction-contact OSNS problem:
 * \f$ velocity = q + M \cdot reaction \f$
 *
 * Becomes:
 * \f$ velocity = (q + V \cdot q_{cohesion}) + M \cdot reaction \f$
 *
 * where the term \f$ V \cdot q_{cohesion} \f$ represents the effect of
 * cohesive forces on the relative velocity at contact points.
 *
 * \par Compatible Nonsmooth Laws
 * - CohesiveZoneModelNIFNSL and derived classes (BinaryCohesiveNSL, etc.)
 * - NewtonImpactFrictionNSL (treated as fully broken cohesive law)
 *
 * \par Simulation Flow
 * 1. preCompute(): Assembles q_cohesion and updates q with cohesive contribution
 * 2. compute(): Solves the modified friction-contact problem
 * 3. postCompute(): Updates internal cohesive state via Interaction::swapInMemory()
 */
class CohesiveFrictionContact : public FrictionContact {
 protected:
  /** \cond DEVEL */
  ACCEPT_SERIALIZATION(CohesiveFrictionContact);
  /** \endcond */

  /** Matrix V for cohesive contribution mapping.
   * 
   * Maps cohesive forces from interactions in indexSet0 (distant but cohesive)
   * to the contact space of the OSNS problem. This is needed because cohesive
   * forces act at different points than the contact points used in the OSNS.
   * 
   * Dimensions: (size of contact problem) x (size of cohesive problem)
   */
  std::shared_ptr<OSNSMatrix> _V{nullptr};

  /** Matrix H0 for direct assembly of cohesive contribution.
   * 
   * Stores the H matrices (Jacobians) for interactions in indexSet0,
   * used to compute the V matrix.
   */
  std::shared_ptr<OSNSMatrix> _H0{nullptr};

  /** Cohesion force contribution vector.
   * 
   * Contains the cohesive forces from all interactions in indexSet0,
   * assembled in the order they appear in the simulation graph.
   * 
   * The cohesive force for each interaction is obtained via
   * CohesiveZoneModelNIFNSL::r_cohesion().
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _q_cohesion{nullptr};

  /** Size of the cohesive problem.
   * 
   * Number of interactions in indexSet0 (distant but still cohesive).
   * This determines the dimension of _q_cohesion.
   */
  siconos::algebra::Index _sizeOutput_cohesion{0};

 public:
  /** \brief Constructor with dimension and solver id
   *
   * \param dimPb dimension of the problem: 2 for 2D, 3 for 3D friction-contact
   * \param numericsSolverId id of the solver to use (e.g., SICONOS_FRICTION_3D_NSGS)
   */
  CohesiveFrictionContact(int dimPb = 3, int numericsSolverId = SICONOS_FRICTION_3D_NSGS);

  /** \brief Constructor with pre-defined solver options
   *
   * \param dimPb dimension of the problem (2 or 3)
   * \param options the solver options structure
   */
  CohesiveFrictionContact(int dimPb, std::shared_ptr<SolverOptions> options);

  /** \brief Destructor */
  ~CohesiveFrictionContact() noexcept override = default;

  /** \brief Initialize the cohesive friction-contact problem
   *
   * Extends FrictionContact::initialize() to allocate additional matrices
   * and vectors needed for cohesive contributions:
   * - Allocates _V matrix
   * - Allocates _H0 matrix  
   * - Allocates _q_cohesion vector
   *
   * \param simulation the simulation that owns this OSNS problem
   */
  void initialize(std::shared_ptr<siconos::simulation::Simulation> simulation) override;

  /** \brief Compute cohesive contribution for a single interaction
   *
   * Computes the contribution of cohesive forces from a single interaction
   * to the global q_cohesion vector. Called during preCompute() for each
   * interaction in indexSet0.
   *
   * \param vertex_inter vertex descriptor for the interaction in the graph
   * \param pos starting position in the global q_cohesion vector
   */
  void computeQCohesionBlock(
      siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
      siconos::algebra::Index pos);

  /** \brief Update q vector with cohesive contribution
   *
   * Modifies the OSNS q vector to include the effect of cohesive forces:
   * \f$ q_{new} = q + V \cdot q_{cohesion} \f$
   *
   * This is called during preCompute() after q_cohesion is assembled.
   *
   * \param time the current simulation time
   */
  void updateQWithQCohesion(double time);

  /** \brief Compute the V matrix for cohesive mapping
   *
   * Builds matrix V that maps cohesive forces from indexSet0 interactions
   * to the contact space. This involves:
   * 1. Gathering H matrices (Jacobians) for all cohesive interactions
   * 2. Computing the mapping to contact space coordinates
   *
   * The V matrix structure depends on the graph of interactions and
   * their geometric relations.
   */
  void computeV();

  /** \brief Build problem coefficients including cohesive contribution
   *
   * Overrides LinearOSNS::preCompute() to:
   * 1. Call parent preCompute() for standard contact terms
   * 2. Compute q_cohesion for all indexSet0 interactions
   * 3. Update q with cohesive contribution
   *
   * \param time the current simulation time
   * \return true if pre-computation succeeded
   */
  bool preCompute(double time) override;

  /** \brief Post-process after OSNS solve
   *
   * Extends FrictionContact::postCompute() to handle cohesive state updates.
   * Ensures that internal variables are properly swapped after convergence.
   */
  void postCompute() override;

  /** \brief Check if a non-smooth law is compatible with this solver
   *
   * Verifies that the given law can be used with CohesiveFrictionContact.
   * Compatible laws include:
   * - CohesiveZoneModelNIFNSL and derived classes
   * - NewtonImpactFrictionNSL (as fallback for broken interfaces)
   *
   * \param nslaw the non-smooth law to check
   * \return true if the law is compatible
   */
  bool checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw) override;

  /** \brief Print problem data to stdout
   *
   * Displays problem statistics including:
   * - Number of contact interactions (indexSet1)
   * - Number of cohesive interactions (indexSet0)
   * - Matrix dimensions
   * - Solver information
   */
  void display() const override;
};

}  // namespace siconos::nonsmooth_formulations

#endif  // COHESIVEFRICTIONCONTACT_H
