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
  Binary cohesive model with impact and friction

  The model implements a naive binary (broken or intact) cohesive zone model.
  When the interface is intact, the limit traction threshold is constant and equal to sigma_c.
  The interface is broken if the normal displacement is larger than delta_c.
*/

#ifndef BINARYCOHESIVENSL_H
#define BINARYCOHESIVENSL_H

#include "CohesiveZoneModelNIFNSL.hpp"
#include "SiconosVector.hpp"

namespace siconos::mechanics::czm {

// Forward declaration
class Interaction;

/** BinaryCohesiveNSL
 * Binary cohesive model with impact and friction
 *
 * Implements a binary (broken or intact) cohesive zone model where:
 * - When intact: constant traction threshold sigma_c
 * - When broken: interface fails if normal displacement > delta_c
 */
class BinaryCohesiveNSL : public siconos::modeling::CohesiveZoneModelNIFNSL {
 public:
  /** Shape type for the cohesive law */
  enum class ShapeType { DOOR_SHAPE, TRIANGLE_SHAPE };

  /** Indices for internal variables */
  enum InternalVariables {
    BETA_SURFACE = 0,              ///< Damage parameter (0=broken, 1=intact) and surface area
    R_COHESION = 1,                ///< Cohesion force vector
    DISPLACEMENT_JUMP = 2,         ///< Current displacement jump
    COHESIVE_POINT_1 = 3,          ///< Contact point on body 1
    COHESIVE_POINT_2 = 4,          ///< Contact point on body 2
    NORMAL = 5,                    ///< Normal vector
    TANGENT_1 = 6,                 ///< First tangent vector
    TANGENT_2 = 7,                 ///< Second tangent vector
    INITIAL_DISPLACEMENT_JUMP = 8, ///< Initial displacement jump
    INITIAL_RELATIVE_COHESIVE_POINT_1 = 9,  ///< Initial relative contact point 1
    INITIAL_RELATIVE_COHESIVE_POINT_2 = 10, ///< Initial relative contact point 2
    INITIAL_RELATIVE_NORMAL = 11,           ///< Initial relative normal
    INITIAL_RELATIVE_TANGENT_1 = 12,        ///< Initial relative tangent 1
    INITIAL_RELATIVE_TANGENT_2 = 13,        ///< Initial relative tangent 2
    INTERNAL_VARIABLE_LENGTH = 14           ///< Total number of internal variables
  };

 private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(BinaryCohesiveNSL);

  /** cohesive resistance to traction */
  double _sigma_c{0.0};

  /** critical displacement */
  double _delta_c{0.0};

  /** slope for triangle law */
  double _slope{0.0};

  /** shape of the cohesive law */
  ShapeType _shape_type{ShapeType::DOOR_SHAPE};

  /** fallback law when the interface is broken */
  std::shared_ptr<siconos::modeling::NonSmoothLaw> _nslaw_broken;

 protected:
  /** default constructor */
  BinaryCohesiveNSL() = delete;  // Base class has deleted default constructor

 public:
  /** basic constructor
   *  \param size size of the ns law
   */
  explicit BinaryCohesiveNSL(siconos::algebra::Index size);

  /** constructor with the value of the BinaryCohesiveNSL attributes
   *  \param en double : normal restitution coefficient
   *  \param et double : tangent restitution coefficient
   *  \param mu double : friction coefficient
   *  \param sigma_c double : cohesive resistance to traction
   *  \param delta_c double : critical displacement
   *  \param size unsigned int: size of the ns law
   */
  BinaryCohesiveNSL(double en, double et, double mu, double sigma_c, double delta_c,
                    siconos::algebra::Index size);

  /** constructor with shape type
   *  \param en double : normal restitution coefficient
   *  \param et double : tangent restitution coefficient
   *  \param mu double : friction coefficient
   *  \param sigma_c double : cohesive resistance to traction
   *  \param delta_c double : critical displacement
   *  \param size unsigned int: size of the ns law
   *  \param shape_type shape of the cohesive law
   */
  BinaryCohesiveNSL(double en, double et, double mu, double sigma_c, double delta_c,
                    siconos::algebra::Index size, ShapeType shape_type);

  /** Destructor */
  ~BinaryCohesiveNSL() noexcept override = default;

  // GETTERS/SETTERS

  /** getter of sigma_c
   * \return the value of sigma_c
   */
  inline double sigma_c() const { return _sigma_c; };

  /** setter of sigma_c
   * \param newVal a double to set sigma_c
   */
  inline void setSigma_c(double newVal) { _sigma_c = newVal; };

  /** getter of delta_c
   * \return the value of delta_c
   */
  inline double delta_c() const { return _delta_c; };

  /** setter of delta_c
   * \param newVal a double to set delta_c
   */
  inline void setDelta_c(double newVal) { _delta_c = newVal; };

  /** getter of shape_type
   * \return the shape type
   */
  inline ShapeType shape_type() const { return _shape_type; };

  // CZM INTERFACE IMPLEMENTATION

  /** Initialize internal variables for this cohesive law
   * \param inter the Interaction to initialize variables for
   * \return shared pointer to vector of internal variables
   */
  std::shared_ptr<siconos::algebra::blocks::SharedVector> initializeInternalVariables(
      siconos::modeling::Interaction& inter) override;

  /** Update internal variables after each time step
   * \param inter the Interaction containing internal variables
   */
  void updateInternalVariables(siconos::modeling::Interaction& inter) override;

  /** Check if the NS law is active at a given level
   * \param inter the Interaction
   * \param level the level to check
   * \return true if active
   */
  bool isActiveAtLevel(siconos::modeling::Interaction& inter, unsigned int level) override;

  /** getter of r_cohesion
   * \param inter the Interaction containing internal variables
   * \return pointer to the cohesion force vector data
   */
  double* r_cohesion(siconos::modeling::Interaction& inter) const override;

  /** getter of beta (damage parameter)
   * \param inter the Interaction containing internal variables
   * \return the value of beta (1.0 = intact, 0.0 = broken)
   */
  double beta(siconos::modeling::Interaction& inter) const;

  /** print the data to the screen
   */
  void display() const override;

  /** Display internal variables for debugging
   * \param inter the Interaction containing internal variables
   */
  void displayInternalVariables(siconos::algebra::blocks::SharedVector & internalVariables)  override;

  /** Visitors hook
   */
  // void accept(siconos::modeling::nonsmooth_laws::Visitor& tourist) const override {
  //   tourist.visit(*this);
  // }

  siconos::modeling::Type acceptType(siconos::types::FindType& ft) const override {
    return ft.visit(*this);
  }
};

}  // namespace siconos::mechanics::czm

#endif  // BINARYCOHESIVENSL_H
