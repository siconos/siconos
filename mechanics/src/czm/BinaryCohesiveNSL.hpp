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
/*! \file BinaryCohesiveNSL.hpp
  \brief Binary cohesive zone model with impact and friction

  \section sec_bczm_overview Overview

  BinaryCohesiveNSL implements a simple binary (intact/broken) cohesive zone model.
  The interface can be in one of two states:
  - <b>Intact</b>: Can sustain traction up to a critical value \f$ \sigma_c \f$
  - <b>Broken</b>: No cohesive traction, standard contact behavior

  This is a simplified cohesive model suitable for brittle fracture simulations
  where progressive damage evolution is not the primary focus.

  \section sec_bczm_criteria Failure Criteria

  The interface transitions from intact to broken when:
  \f[
  \delta_n > \delta_c
  \f]
  where:
  - \f$ \delta_n \f$ is the normal displacement jump across the interface
  - \f$ \delta_c \f$ is the critical displacement for failure

  \section sec_bczm_shapes Cohesive Law Shapes

  Two shape types are available:

  \par DOOR_SHAPE (default)
  Constant traction until failure:
  \f[
  t_c = \begin{cases}
    \sigma_c & \text{if } \delta_n \leq \delta_c \text{ and intact} \\
    0 & \text{if broken}
  \end{cases}
  \f]

  \par TRIANGLE_SHAPE
  Linear softening after initial elastic phase:
  \f[
  t_c = \begin{cases}
    \sigma_c (1 - \delta_n/\delta_c) & \text{if } \delta_n \leq \delta_c \text{ and intact} \\
    0 & \text{if broken}
  \end{cases}
  \f]

  \section sec_bczm_variables Internal Variables

  The model stores 14 internal variables per interaction:
  - BETA_SURFACE: Damage parameter \f$ \beta \f$ (1=intact, 0=broken) and surface area
  - R_COHESION: Cohesive force vector (3 components)
  - DISPLACEMENT_JUMP: Current displacement jump across interface (3 components)
  - COHESIVE_POINT_1/2: Contact points on each body (3 components each)
  - NORMAL/TANGENT_1/TANGENT_2: Local coordinate frame (3 components each)
  - INITIAL_*: Initial values for restart and persistency

  \section sec_bczm_usage Usage Example

  \code
  #include "BinaryCohesiveNSL.hpp"

  // Create a 3D cohesive law with door-shaped profile
  auto nslaw = std::make_shared<siconos::mechanics::czm::BinaryCohesiveNSL>(
      0.0,                    // en: normal restitution
      0.0,                    // et: tangent restitution
      0.3,                    // mu: friction coefficient
      10.0e6,                 // sigma_c: 10 MPa critical traction
      0.1e-3,                 // delta_c: 0.1 mm critical displacement
      3,                      // size: 3D
      BinaryCohesiveNSL::ShapeType::DOOR_SHAPE
  );

  // Use with a contact relation
  auto relation = std::make_shared<NewtonEuler1DR>(...);
  auto inter = std::make_shared<Interaction>(nslaw, relation);
  \endcode

  \see CohesiveZoneModelNIFNSL for the base class
  \see CohesiveFrictionContact for the OSNS solver
*/

#ifndef BINARYCOHESIVENSL_H
#define BINARYCOHESIVENSL_H

#include "CohesiveZoneModelNIFNSL.hpp"
#include "SiconosVector.hpp"

namespace siconos::mechanics::czm {

// Forward declaration
class Interaction;

/** \class BinaryCohesiveNSL
 * \brief Binary (intact/broken) cohesive zone model with impact and friction
 *
 * Implements a simple cohesive zone model where the interface is either fully
 * intact (can sustain traction up to \f$ \sigma_c \f$) or fully broken (no
 * cohesive traction). The transition occurs when the displacement jump exceeds
 * the critical value \f$ \delta_c \f$.
 *
 * \par Cohesive Parameters
 * - \f$ \sigma_c \f$ (sigma_c): Critical traction stress (Pa). Maximum traction
 *   the interface can sustain when intact.
 * - \f$ \delta_c \f$ (delta_c): Critical displacement (m). Displacement jump at
 *   which complete failure occurs.
 *
 * \par Shape Types
 * - DOOR_SHAPE: Constant traction \f$ \sigma_c \f$ until failure (default)
 * - TRIANGLE_SHAPE: Linear degradation from \f$ \sigma_c \f$ to 0 over \f$ \delta_c \f$
 *
 * \par Internal Variables Storage
 * The model uses Interaction::internalVariables() to store:
 * - Current damage state (beta)
 * - Cohesive force vector
 * - Displacement jump history
 * - Contact geometry (points, normals, tangents)
 *
 * These variables are automatically initialized when the interaction is created
 * and updated at each time step by updateInternalVariables().
 */
class BinaryCohesiveNSL : public siconos::modeling::CohesiveZoneModelNIFNSL {
 public:
  /** \brief Shape type for the cohesive traction-separation law
   *
   * Defines the functional form of the cohesive traction as a function
   * of the displacement jump.
   */
  enum class ShapeType {
    DOOR_SHAPE,     ///<< Constant traction until failure (step function)
    TRIANGLE_SHAPE  ///<< Linear softening from sigma_c to 0
  };

  /** \brief Indices for internal variables in the SharedVector storage
   *
   * These indices map to positions in the internal variables vector
   * stored in the Interaction. Each entry corresponds to a
   * std::shared_ptr<SiconosVector> containing the specific variable data.
   */
  enum InternalVariables {
    BETA_SURFACE = 0,               ///<< [0]=damage beta (1=intact,0=broken), [1]=surface area
    COHESION = 1,                   ///<< Cohesive force intensity (3D)
    DISPLACEMENT_JUMP = 2,          ///<< Current displacement jump vector (3D)
    COHESIVE_POINT_1 = 3,           ///<< Contact point on body 1 in global frame (3D)
    COHESIVE_POINT_2 = 4,           ///<< Contact point on body 2 in global frame (3D)
    NORMAL = 5,                     ///<< Normal vector at contact point (3D)
    TANGENT_1 = 6,                  ///<< First tangent vector (3D)
    TANGENT_2 = 7,                  ///<< Second tangent vector (3D)
    INITIAL_DISPLACEMENT_JUMP = 8,  ///<< Initial displacement jump for reference (3D)
    INITIAL_RELATIVE_COHESIVE_POINT_1 = 9,   ///<< Initial relative contact point 1 (3D)
    INITIAL_RELATIVE_COHESIVE_POINT_2 = 10,  ///<< Initial relative contact point 2 (3D)
    INITIAL_RELATIVE_NORMAL = 11,            ///<< Initial relative normal vector (3D)
    INITIAL_RELATIVE_TANGENT_1 = 12,         ///<< Initial relative tangent 1 (3D)
    INITIAL_RELATIVE_TANGENT_2 = 13,         ///<< Initial relative tangent 2 (3D)
    INTERNAL_VARIABLE_LENGTH = 14            ///<< Total number of internal variables
  };

 private:
  /** \cond DEVEL */
  /** serialization hooks */
  ACCEPT_SERIALIZATION(BinaryCohesiveNSL);
  /** \endcond */

  /** Critical traction stress (Pa).
   * Maximum cohesive traction the interface can sustain when intact.
   * When the traction exceeds this value, damage initiates.
   */
  double _sigma_c{0.0};

  /** Critical displacement (m).
   * Displacement jump at which the interface is considered fully broken.
   * For DOOR_SHAPE: instant failure at \f$ \delta_c \f$
   * For TRIANGLE_SHAPE: linear softening from 0 to \f$ \delta_c \f$
   */
  double _delta_c{0.0};

  /** ratio normal tangent.
   */
  double _gamma{0.0};

  /** Slope for TRIANGLE_SHAPE law.
   * Computed as \f$ -1/\delta_c \f$ during construction.
   */
  double _slope{0.0};

  /** Shape type of the cohesive law.
   * Determines the traction-separation curve shape.
   */
  ShapeType _shape_type{ShapeType::DOOR_SHAPE};

  /** Fallback law when the interface is broken.
   * Inherited from base class, redeclared here for clarity.
   */
  std::shared_ptr<siconos::modeling::NonSmoothLaw> _nslaw_broken;

 protected:
  /** Default constructor is deleted.
   * Use the parameterized constructors instead.
   */
  BinaryCohesiveNSL() = delete;

 public:
  /** \brief Constructor with size only
   *
   * Creates a cohesive law with zero parameters. Useful for deserialization
   * or when parameters will be set via setters.
   *
   * \param size dimension of the non-smooth law (2 or 3)
   */
  explicit BinaryCohesiveNSL(siconos::algebra::Index size);

  /** \brief Constructor with cohesive parameters (default DOOR_SHAPE)
   *
   * \param en normal restitution coefficient [0,1]
   * \param et tangent restitution coefficient [0,1]
   * \param mu friction coefficient (>= 0)
   * \param sigma_c critical traction stress (Pa), must be >= 0
   * \param delta_c critical displacement (m), must be > 0
   * \param size dimension of the law (2 for 2D, 3 for 3D contact)
   * \param gamma ratio of normal to tangent cohesion (default 0.0)
   */
  BinaryCohesiveNSL(double en, double et, double mu, double sigma_c, double delta_c,
                    siconos::algebra::Index size, double gamma = 0.0);

  /** \brief Constructor with shape type selection
   *
   * \param en normal restitution coefficient [0,1]
   * \param et tangent restitution coefficient [0,1]
   * \param mu friction coefficient (>= 0)
   * \param sigma_c critical traction stress (Pa), must be >= 0
   * \param delta_c critical displacement (m), must be > 0
   * \param size dimension of the law (2 for 2D, 3 for 3D contact)
   * \param shape_type shape of the traction-separation law
   * \param gamma ratio of normal to tangent cohesion (default 0.0)
   */
  BinaryCohesiveNSL(double en, double et, double mu, double sigma_c, double delta_c,
                    siconos::algebra::Index size, ShapeType shape_type, double gamma = 0.0);

  /** \brief Destructor */
  ~BinaryCohesiveNSL() noexcept override = default;

  // GETTERS/SETTERS

  /** \brief Get the critical traction stress
   * \return critical traction \f$ \sigma_c \f$ in Pascals
   */
  inline double sigma_c() const { return _sigma_c; };

  /** \brief Set the critical traction stress
   * \param newVal new critical traction value (Pa), must be >= 0
   */
  inline void setSigma_c(double newVal) { _sigma_c = newVal; };

  /** \brief Get the critical displacement
   * \return critical displacement \f$ \delta_c \f$ in meters
   */
  inline double delta_c() const { return _delta_c; };

  /** \brief Set the critical displacement
   * \param newVal new critical displacement value (m), must be > 0
   */
  inline void setDelta_c(double newVal) { _delta_c = newVal; };

  /** \brief Get the shape type of the cohesive law
   * \return the shape type (DOOR_SHAPE or TRIANGLE_SHAPE)
   */
  inline ShapeType shape_type() const { return _shape_type; };

  /** \brief Get the normal/tangent cohesion ratio
   * \return gamma ratio (normal cohesion / tangent cohesion)
   */
  inline double gamma() const { return _gamma; };

  /** \brief Set the normal/tangent cohesion ratio
   * \param newVal new gamma value, must be >= 0
   */
  inline void setGamma(double newVal) { _gamma = newVal; };

  // CZM INTERFACE IMPLEMENTATION

  /** \brief Initialize internal variables for this cohesive law
   *
   * Called automatically when an Interaction using this law is initialized.
   * Allocates storage for all internal variables and sets initial values:
   * - beta = 1.0 (fully intact)
   * - surface area = 1.0 (should be set correctly based on geometry)
   * - cohesive force = zero vector
   * - displacement jump = zero vector
   * - stores initial contact geometry for reference
   *
   * \param inter the Interaction to initialize variables for
   * \return shared pointer to the initialized internal variables vector
   * \note This method requires the Interaction to have a NewtonEuler1DR relation
   */
  std::shared_ptr<siconos::algebra::blocks::SharedVector3> initializeInternalVariables(
      siconos::modeling::Interaction& inter) override;

  /** \brief Update internal variables after each time step
   *
   * This method is called by the simulation after each successful time step
   * to update the cohesive state:
   * - Computes current displacement jump from contact geometry
   * - Updates damage parameter beta based on failure criterion
   * - Computes cohesive traction force
   * - Stores current state for next time step
   *
   * \param inter the Interaction containing internal variables
   * \note The damage evolution is irreversible (beta only decreases)
   */
  void updateInternalVariables(siconos::modeling::Interaction& inter) override;

  /** \brief Check if the cohesive law is active at a given level
   *
   * Cohesive forces are computed at level 1 (predictor step) to influence
   * the contact detection and reaction computation at level 0.
   *
   * \param inter the Interaction
   * \param level the level to check (0 = corrector, 1 = predictor)
   * \return true if the law should be applied at this level
   * \note Returns true only for level == 1 and when beta > 0
   */
  bool isActiveAtLevel(siconos::modeling::Interaction& inter, unsigned int level) override;

  /** \brief Get the cohesive force vector
   *
   * Returns a pointer to the cohesive force vector stored in the internal
   * variables. This is used by the OSNS solver to add cohesive contributions.
   *
   * \param inter the Interaction containing internal variables
   * \return pointer to the 3D cohesive force vector data
   * \note The force is computed during updateInternalVariables()
   */
  double* cohesion(siconos::modeling::Interaction& inter) const override;

  /** \brief Get the damage parameter beta
   *
   * \param inter the Interaction containing internal variables
   * \return damage parameter beta (1.0 = fully intact, 0.0 = fully broken)
   * \note Beta decreases monotonically during the simulation
   */
  double beta(siconos::modeling::Interaction& inter) const;

  /** \brief Print the law's parameters to stdout
   *
   * Displays sigma_c, delta_c, shape_type, and other parameters.
   */
  void display() const override;

  /** \brief Display internal variables for debugging
   *
   * Prints the current values of all internal variables including:
   * - beta (damage parameter)
   * - cohesive force vector
   * - displacement jump
   * - contact geometry
   *
   * \param internalVariables the internal variables vector to display
   */
  void displayInternalVariables(
      siconos::algebra::blocks::SharedVector3& internalVariables) override;

  /** \cond DEVEL */
  /** Visitors hook for type dispatch */
  siconos::modeling::Type acceptType(siconos::types::FindType& ft) const override {
    return ft.visit(*this);
  }
  /** \endcond */
};

}  // namespace siconos::mechanics::czm

#endif  // BINARYCOHESIVENSL_H
