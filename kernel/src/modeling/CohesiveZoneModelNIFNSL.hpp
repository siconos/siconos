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
/*! \file CohesiveZoneModelNIFNSL.hpp
  \brief Base class for cohesive zone models based on NewtonImpactFrictionNSL

  \section sec_czm_overview Overview

  A cohesive zone model (CZM) describes the progressive damage and failure of an 
  interface between two materials. Unlike standard contact laws that only act when 
  bodies are in contact (gap <= 0), cohesive laws can sustain traction forces 
  when the interface is intact (gap > 0 but small).

  \section sec_czm_physics Physical Behavior

  The cohesive interface evolves through several stages:
  - <b>Intact state</b>: The interface can sustain traction up to a critical value.
    The traction force increases as the gap opens (elastic or rigid behavior).
  - <b>Damage initiation</b>: When traction reaches the critical value, damage 
    begins to accumulate.
  - <b>Softening</b>: The traction capacity decreases as the gap continues to open.
  - <b>Complete failure</b>: When the gap exceeds a critical displacement, the 
    interface is fully broken and the standard contact law takes over.

  \section sec_czm_mathematical Mathematical Formulation

  The cohesive traction \f$ t_c \f$ is typically a function of the displacement 
  jump \f$ \delta \f$ across the interface:
  \f[
  t_c = f(\delta, \text{damage parameters})
  \f]

  The total traction is the sum of:
  - Cohesive traction (when interface is not fully broken)
  - Contact reaction (when gap <= 0)

  \section sec_czm_implementation Implementation Notes

  This base class extends NewtonImpactFrictionNSL to add:
  - Internal variables storage (via Interaction) for damage state
  - Pure virtual methods for cohesive force computation
  - A fallback law (_nslaw_broken) used when the interface is fully damaged

  Derived classes must implement:
  - initializeInternalVariables(): Set up storage for damage parameters
  - updateInternalVariables(): Update damage state after each time step
  - r_cohesion(): Compute the cohesive force vector

  \see BinaryCohesiveNSL for a concrete implementation
  \see Interaction::internalVariables() for storage mechanism
  \see CohesiveFrictionContact for the OSNS problem formulation
*/

#ifndef COHESIVEZONEMODELNIFNSL_H
#define COHESIVEZONEMODELNIFNSL_H

#include "NewtonImpactFrictionNSL.hpp"
#include "SiconosVector.hpp"

namespace siconos::modeling {

// Forward declaration
class Interaction;

/** \class CohesiveZoneModelNIFNSL
 * \brief Abstract base class for cohesive zone models extending NewtonImpactFrictionNSL
 * 
 * This class extends the Newton impact-friction law with cohesive behavior,
 * where the interface can sustain traction up to a critical value before
 * progressive damage occurs. It serves as the base class for all cohesive
 * zone models in Siconos.
 *
 * \par Internal Variables
 * Cohesive models require internal state variables (damage parameters) that
 * persist across time steps. These are stored in the Interaction object using
 * Interaction::internalVariables() and Interaction::internalVariables_k().
 *
 * \par Fallback Law
 * When the interface is fully damaged, the cohesive contribution vanishes and
 * the simulation falls back to a standard NewtonImpactFrictionNSL (_nslaw_broken).
 *
 * \par Usage Example
 * \code
 * // Create a cohesive law
 * auto nslaw = std::make_shared<BinaryCohesiveNSL>(
 *     0.0,      // en: normal restitution
 *     0.0,      // et: tangent restitution
 *     0.3,      // mu: friction coefficient
 *     1e6,      // sigma_c: critical traction (Pa)
 *     1e-4,     // delta_c: critical displacement (m)
 *     3         // size: 3D problem
 * );
 * 
 * // Create interaction with cohesive law
 * auto inter = std::make_shared<Interaction>(nslaw, relation);
 * \endcode
 */
class CohesiveZoneModelNIFNSL : public NewtonImpactFrictionNSL {
 private:
  /** serialization hooks
   * \cond DEVEL
   */
  ACCEPT_SERIALIZATION(CohesiveZoneModelNIFNSL);
  /** \endcond */

  /** Fallback non-smooth law used when the interface is fully broken.
   * This is typically a NewtonImpactFrictionNSL with the same restitution
   * and friction parameters as the original cohesive law.
   */
  std::shared_ptr<NonSmoothLaw> _nslaw_broken;

 protected:
  /** Default constructor is deleted.
   * Use the parameterized constructors instead.
   */
  CohesiveZoneModelNIFNSL() = delete;

 public:
  /** \brief Constructor with size only
   * 
   * Creates a cohesive law with default parameters (zero restitution, 
   * zero friction). Useful when parameters will be set later.
   * 
   * \param size dimension of the non-smooth law (2 for 2D, 3 for 3D)
   */
  explicit CohesiveZoneModelNIFNSL(siconos::algebra::Index size);

  /** \brief Constructor with full parameters
   * 
   * \param en normal restitution coefficient (0 = perfectly inelastic, 1 = elastic)
   * \param et tangent restitution coefficient
   * \param mu friction coefficient
   * \param size dimension of the non-smooth law (2 for 2D, 3 for 3D)
   */
  CohesiveZoneModelNIFNSL(double en, double et, double mu, siconos::algebra::Index size);

  /** \brief Virtual destructor */
  ~CohesiveZoneModelNIFNSL() noexcept override = default;

  /** \brief Compute the cohesive force vector
   * 
   * This pure virtual method must be implemented by derived classes to
   * compute the cohesive traction force based on the current state of
   * the interface (displacement jump, damage parameters, etc.).
   * 
   * \param inter the Interaction containing internal variables (damage state)
   * \return pointer to the cohesive force vector data (typically size 3)
   * \note The returned pointer points to internal storage that is valid 
   *       until the next call to updateInternalVariables()
   */
  virtual double* cohesion(Interaction& inter) const = 0;

  /** \brief Get the fallback law for broken interfaces
   * 
   * When the interface is fully damaged, this law replaces the cohesive
   * behavior with standard contact mechanics.
   * 
   * \return shared pointer to the fallback non-smooth law
   */
  std::shared_ptr<NonSmoothLaw> nslawBroken() const { return _nslaw_broken; };

  /** \brief Check if the NS law is active at a given level
   * 
   * Cohesive laws are typically active at level 1 (predictor step) to
   * compute cohesive forces before the contact detection at level 0.
   * 
   * \param inter the Interaction
   * \param level the level to check (0 = corrector, 1 = predictor)
   * \return true if the law should be applied at this level
   * \note The default implementation returns true for level == 1
   */
  bool isActiveAtLevel(Interaction& inter, unsigned int level) const override;

  /** \brief Print the law's data to the screen
   * 
   * Displays the parameters of the cohesive law for debugging purposes.
   */
  void display() const override;

  /** \brief Visitor pattern support
   * \cond DEVEL
   */
  void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
  /** \endcond */
};

}  // namespace siconos::modeling

#endif  // COHESIVEZONEMODELNIFNSL_H
