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
  cohesive zone models (CZM) in a unified velocity-displacement formulation.
  While standard contact only acts when bodies overlap (gap <= 0), cohesive
  forces act when the interface is intact and can sustain traction before
  complete failure.

  This is essential for simulating:
  - Delamination of composite materials
  - Fracture mechanics with crack propagation
  - Adhesive contact and debonding processes
  - Rock mechanics with joint failure
  - Concrete cracking and damage

  \section sec_cfc_mathematical Mathematical Formulation

  The cohesive friction-contact problem uses a coupled velocity-displacement
  formulation with block matrix structure:

  \f[
  \underbrace{\begin{bmatrix} M & V \\ U & X \end{bmatrix}}_{\text{Mass/Stiffness matrix}}
  \underbrace{\begin{bmatrix} r_v \\ r_u \end{bmatrix}}_{\text{Reaction}}
  +
  \underbrace{\begin{bmatrix} q_v \\ q_u \end{bmatrix}}_{\text{Input}}
  =
  \underbrace{\begin{bmatrix} u_v \\ u_u \end{bmatrix}}_{\text{Velocity}}
  \f]

  Subject to constraints:
  \f[
  \begin{cases}
    \text{Contact:} & \text{gap} \leq 0 \perp r_{v,n} \geq 0, \quad \|r_{v,t}\| \leq \mu \cdot
r_{v,n} \\
    \text{Cohesion:} & r_{u} = f(\text{damage}, \text{history}) \quad \text{(from CZM law)}
  \end{cases}
  \f]

  where:
  - \f$ r_v \f$: Contact reactions (velocity-level, at active contacts)
  - \f$ r_u \f$: Cohesive reactions (displacement-level, at potential cohesive zones)
  - \f$ u_v \f$: Relative velocities at contact points
  - \f$ u_u \f$: Relative displacements at cohesive points
  - \f$ q_v \f$: Input velocities (from predictor step)
  - \f$ q_u \f$: Input displacements (history variables)

  \section sec_cfc_blocks Block Matrix Structure

  The block matrices have specific physical meanings:

  \par W (contact-contact block)
  Mass/stiffness coupling between contact degrees of freedom.
  Equivalent to the M matrix in standard FrictionContact.
  Size: \f$ (d \cdot n_c) \times (d \cdot n_c) \f$

  \par V (contact-cohesive coupling)
  Coupling from cohesive zones to contact points. Represents how
  cohesive deformation affects the contact kinematics.
  Size: \f$ (d \cdot n_c) \times (d \cdot n_{coh}) \f$

  \par U (cohesive-contact coupling)
  Coupling from contact points to cohesive zones. Represents how
  contact deformation affects cohesive zone kinematics.
  Size: \f$ (d \cdot n_{coh}) \times (d \cdot n_c) \f$

  \par X (cohesive-cohesive block)
  Stiffness of cohesive zones. Related to the cohesive law stiffness
  and damage state.
  Size: \f$ (d \cdot n_{coh}) \times (d \cdot n_{coh}) \f$

  \section sec_cfc_algorithm Algorithm

  The solution process at each time step:

  \par 1. Update contact detection
  - Run broadphase/narrowphase collision detection
  - Update indexSet1 (active contacts)
  - Update indexSet0 (distant but potentially cohesive)

  \par 2. Compute cohesive intensity/threshold
  For each interaction in indexSet0:
  - Get internal cohesive state (damage, history) c_n c_t

  \par 3. Assemble block matrices
  Build W, V, U, X from:
  - System mass matrix (H matrices from relations)
  - Cohesive law stiffness
  - Time integration parameters

  \par 4. Form global problem
  Construct M = [[W, V], [U, X]] and q = [q_v, q_u]

  \par 5. Solve coupled problem
  Solve for r = [r_v, r_u] with:
  - Friction cone constraints on r_v
  - Cohesive law constraints on r_u

  \par 6. Update cohesive state
  - Update damage parameters based on displacement/velocity
  - Store history via Interaction::swapInMemory()

  \section sec_cfc_index_sets Index Sets

  The solver manages two index sets:

  \par indexSet1 (Contact Index Set)
  Interactions with physical contact (gap <= 0).
  These contribute to the contact part of the problem (r_v).

  \par indexSet0 (Cohesive Index Set)
  Interactions without physical contact but with potential cohesion.
  These contribute to the cohesive part (r_u).
  An interaction can be in both sets simultaneously.

  \section sec_cfc_nslaws Compatible Nonsmooth Laws

  \par Primary: CohesiveZoneModelNIFNSL
  Base class for cohesive zone models. Derived classes include:
  - BinaryCohesiveNSL: Bilinear/exponential cohesive law
  - Can handle damage evolution and softening

  \par Fallback: NewtonImpactFrictionNSL To be implemented !!
  Standard contact law used when cohesion is completely broken
  (damage = 1, zero traction capacity).

  \section sec_cfc_usage Usage Example

  \code
  #include "CohesiveFrictionContact.hpp"
  #include "BinaryCohesiveNSL.hpp"

  // Create the cohesive friction-contact solver
  auto osnsp = std::make_shared<siconos::nonsmooth_formulations::CohesiveFrictionContact>(
      3,                          // dimension: 3D
      SICONOS_COHESIVE_FRICTION_3D_NSGS  // solver: NSGS for cohesive
  );

  // Configure solver options
  osnsp->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  osnsp->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;

  // Create simulation with cohesive solver
  auto simulation = std::make_shared<TimeStepping>(td, osi, osnsp);

  // Create cohesive interactions
  auto nslaw = std::make_shared<BinaryCohesiveNSL>(
      0.0,     // en: restitution coefficient (normal)
      0.0,     // et: restitution coefficient (tangent)
      0.3,     // mu: friction coefficient
      10.0e6,  // sigma_c: critical stress (10 MPa)
      0.1e-3,  // delta_c: critical opening (0.1 mm)
      3        // dimension: 3D
  );

  // Create relation and interaction
  auto relation = std::make_shared<NewtonEulerVelocityDisplacementR>();
  auto inter = std::make_shared<Interaction>(nslaw, relation);

  // Add to dynamical systems
  // The cohesive force will be computed automatically each time step
  \endcode

  \section sec_cfc_references References

  - Camacho, G. T., & Ortiz, M. (1996). Computational modelling of impact damage
    in brittle materials. International Journal of Solids and Structures.
  - Geubelle, P. H., & Baylor, J. S. (1998). Impact-induced delamination of
    composites: A 2D simulation. Composites Part B.
  - Acary, V., & Monerie, Y. (2006). Nonsmooth fracture dynamics using a cohesive
    zone approach (Technical report, INRIA).
  - Collins-Craft, N. A., Bourrier, F., & Acary, V. (2022). On the formulation and
    implementation of extrinsic cohesive zone models with contact. Computer Methods
    in Applied Mechanics and Engineering, 400, 115545.
  - Collins-Craft, N. A., & Acary, V. (2025). On the formulation and implementation of mixed mode I and mode II extrinsic cohesive zone models with contact and friction.

  \see FrictionContact for the base friction-contact solver
  \see CohesiveZoneModelNIFNSL for cohesive law interface
  \see BinaryCohesiveNSL for bilinear cohesive law implementation
  \see NewtonEulerVelocityDisplacementR for compatible relation
*/

#ifndef COHESIVEFRICTIONCONTACT_H
#define COHESIVEFRICTIONCONTACT_H

#include "CohesiveFrictionContact_options.h"
#include "FrictionContact.hpp"

//struct FrictionContactProblem;
struct CohesiveFrictionContactProblem;
struct SolverOptions;

namespace siconos::nonsmooth_formulations {

/** \class CohesiveFrictionContact
 * \brief Solver for friction-contact problems with cohesive zone models
 *
 * This class extends FrictionContact to handle cohesive forces that act
 * in the pre-failure regime of interfaces. It uses a coupled velocity-
 * displacement formulation with block matrix structure M = [[W,V],[U,X]].
 *
 * \par Problem Structure
 * The OSNS problem is formulated as:
 * \f$ M \cdot r + q = u \f$ with
 * - M = [[W, V], [U, X]] (block mass/stiffness matrix)
 * - r = [r_v, r_u] (contact and cohesive reactions)
 * - q = [q_v, q_u] (input velocity and displacement)
 * - u = [u_v, u_u] (output velocity and displacement rate)
 *
 * \par Key Features
 * - Unified treatment of contact and cohesive forces
 * - Block matrix structure for coupled velocity-displacement
 * - Damage evolution via internal cohesive state
 * - Compatible with standard contact laws (as limiting case)
 *
 * \par Compatible Nonsmooth Laws
 * - CohesiveZoneModelNIFNSL and derived classes (BinaryCohesiveNSL, etc.)
 * - NewtonImpactFrictionNSL (treated as fully damaged/broken interface)
 *
 * \par Simulation Flow
 * 1. initialize(): Allocate matrices W, V, U, X and vectors
 * 2. preCompute(): Assemble block matrices and form global problem
 * 3. compute(): Solve coupled problem for r = [r_v, r_u]
 * 4. postCompute(): Update damage state via swapInMemory()
 */
class CohesiveFrictionContact : public FrictionContact {
 protected:
  /** \cond DEVEL */
  ACCEPT_SERIALIZATION(CohesiveFrictionContact);
  /** \endcond */

  typedef int (*Driver)(CohesiveFrictionContactProblem *, double *, double *, SolverOptions *);

  /** Pointer to the Numerics driver function for cohesive friction-contact
   *
   * The driver solves the problem:
   * M * r + q = velocity, subject to friction and cohesive constraints
   */
  Driver _cohesiveFrictionContact_driver;

  /** Normal cohesion intensity vector (size: numberOfCohesivePoints).
   *
   * Contains the normal component of cohesive traction for each cohesive
   * interaction in indexSet0. Updated each time step from the cohesive law.
   */
  std::shared_ptr<std::vector<double>> _c_n{nullptr};

  /** Tangential cohesion intensity vector (size: numberOfCohesivePoints).
   *
   * Contains the tangential component of cohesive traction for each cohesive
   * interaction in indexSet0. Updated each time step from the cohesive law.
   */
  std::shared_ptr<std::vector<double>> _c_t{nullptr};

  /** Matrix V: contact-cohesive coupling block.
   *
   * Maps cohesive reactions to contact velocity space.
   * Part of the global matrix M = [[W, V], [U, X]].
   *
   * Dimensions: (d*n_c) x (d*n_coh)
   * where d=dimension, n_c=contacts, n_coh=cohesive points
   */
  std::shared_ptr<OSNSMatrix> _V{nullptr};

  /** Matrix U: cohesive-contact coupling block.
   *
   * Maps contact reactions to cohesive displacement space.
   * Part of the global matrix M = [[W, V], [U, X]].
   *
   * Dimensions: (d*n_coh) x (d*n_c)
   */
  std::shared_ptr<OSNSMatrix> _U{nullptr};

  /** Matrix X: cohesive-cohesive block.
   *
   * Stiffness of cohesive zones (from cohesive law).
   * Part of the global matrix M = [[W, V], [U, X]].
   *
   * Dimensions: (d*n_coh) x (d*n_coh)
   */
  std::shared_ptr<OSNSMatrix> _X{nullptr};

  /** Matrix H0: Jacobian storage for cohesive interactions.
   *
   * Stores the H matrices (Jacobians) for interactions in indexSet0,
   * used to compute the coupling blocks V, U, and X.
   */
  std::shared_ptr<OSNSMatrix> _H0{nullptr};

  /** Cohesive reaction vector r_u (displacement-level).
   *
   * Contains the cohesive reactions from all interactions in indexSet0,
   * assembled as r_u = [r_{u,1}, r_{u,2}, ..., r_{u,n_coh}].
   * Each block has dimension components (3 for 3D, 2 for 2D).
   *
   * These reactions are computed by the cohesive law based on
   * damage state and opening displacement.
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _q_cohesion{nullptr};

  /** Size of the cohesive problem (number of cohesive blocks).
   *
   * Number of interactions in indexSet0 (potential cohesive zones).
   * Determines the dimension of r_u: d * _sizeOutput_cohesion.
   */
  siconos::algebra::Index _sizeOutput_cohesion{0};

  bool _scaling_as_percussion = false;



 public:
  /** \brief Constructor with dimension and solver id
   *
   * \param dimPb dimension of the problem: 2 for 2D, 3 for 3D friction-contact
   * \param numericsSolverId id of the solver (e.g., SICONOS_COHESIVE_FRICTION_3D_NSGS)
   */
  CohesiveFrictionContact(int dimPb = 3, int numericsSolverId = SICONOS_COHESIVE_FRICTION_3D_NSGS);

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
   * Extends FrictionContact::initialize() to allocate block matrices
   * and vectors needed for the coupled formulation:
   * - Allocates _V, _U, _X coupling matrices
   * - Allocates _H0 for cohesive Jacobians
   * - Allocates _q_cohesion for cohesive reactions
   *
   * \param simulation the simulation that owns this OSNS problem
   */
  void initialize(std::shared_ptr<siconos::simulation::Simulation> simulation) override;

  /** \brief Update friction and cohesion coefficients
   *
   * Updates the mu (friction), c_n (normal cohesion), and c_t (tangential cohesion)
   * vectors from the non-smooth laws of all interactions.
   */
  void updateCoefficients();

  /** \return the cohesive friction contact problem from Numerics
   *
   * Returns the C-struct problem used by the Numerics solver,
   * containing M, q, mu, c_n, c_t, and solver options.
   */
  std::shared_ptr<CohesiveFrictionContactProblem> cohesiveFrictionContactProblem();

  /** \brief Compute cohesive reaction for a single interaction
   *
   * Computes the cohesive reaction contribution from one interaction
   * in indexSet0 and stores it in _q_cohesion at the specified position.
   * Called during preCompute() for each cohesive interaction.
   *
   * \param vertex_inter vertex descriptor for the interaction in the graph
   * \param pos starting position in _q_cohesion vector (in blocks)
   */
  void compute_q_cohesion_block(
      siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
      siconos::algebra::Index pos);

  /** \brief Assemble and update q vector with cohesive contribution
   *
   * This is the main pre-computation step that:
   * 1. Computes r_u for all cohesive interactions
   * 2. Assembles the global reaction vector r = [r_v, r_u]
   * 3. Forms the effective q = [q_v, q_u]
   *
   * \param time the current simulation time
   */
  void compute_q_cohesion(double time);

  /** \brief Compute the block matrices W, V, U, X
   *
   * Builds the block matrices of the global system:
   * - W: contact-contact (extends standard M)
   * - V: contact-cohesive coupling
   * - U: cohesive-contact coupling
   * - X: cohesive-cohesive stiffness
   *
   * These are assembled from H matrices (Jacobians) of all interactions.
   */
  void computeMatrices();

  /** \brief Pre-computation step
   *
   * Overrides LinearOSNS::preCompute() to:
   * 1. Call parent for standard contact terms (W block)
   * 2. Compute block matrices V, U, X
   * 3. Compute cohesive reactions r_u
   * 4. Assemble global problem M*r + q = u
   *
   * \param time the current simulation time
   * \return true if pre-computation succeeded
   */
  bool preCompute(double time) override;

  /** \brief Post-process after OSNS solve
   *
   * Extends FrictionContact::postCompute() to update cohesive state.
   * Ensures damage variables are properly updated via swapInMemory().
   */
  void postCompute() override;

  /** \brief Solve the cohesive friction-contact problem
   *
   * Calls the Numerics driver to solve:
   * M * r + q = u with friction and cohesive constraints
   *
   * \return solver convergence info (0: success, >0: failure)
   */
  int solve();

  /** \brief Main compute function
   *
   * Solves for unknown reactions r = [r_v, r_u] and velocities u = [u_v, u_u],
   * then updates the Interaction y (velocity) and lambda (reaction) variables.
   *
   * \param time the current time
   * \return solver convergence info (0: ok, >0: problem)
   */
  int compute(double time) override;

  /** \brief Check if a non-smooth law is compatible
   *
   * Verifies that the given law can be used with this solver.
   * Compatible laws:
   * - CohesiveZoneModelNIFNSL and derived classes
   * - NewtonImpactFrictionNSL (fallback for broken interfaces)
   *
   * \param nslaw the non-smooth law to check
   * \return true if compatible
   */
  bool checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw) override;

  /** \brief Print problem statistics
   *
   * Displays:
   * - Number of contacts (indexSet1 size)
   * - Number of cohesive zones (indexSet0 size)
   * - Block matrix dimensions
   * - Solver configuration
   */
  void display() const override;
};

}  // namespace siconos::nonsmooth_formulations

#endif  // COHESIVEFRICTIONCONTACT_H
