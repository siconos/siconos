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
#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP
#include <memory>
#include <vector>

#include "SiconosVector.hpp"
#include "SiconosSerialization.hpp"

namespace siconos::plugins {
class PluggedObject;
}

namespace siconos::algebra {
// class SiconosVector;
}

namespace siconos::modeling {

typedef void (*FPtrPrescribedVelocity)(double, unsigned int, double*);

/** This class models simple boundary conditions for
 *  prescribing the velocities in a Dynamical System. A simple
 *  boundary condition is considered to fix a component \f$ j \f$ of
 *  the velocity vector, i.e., \f$ v_j(t) = bc(t) \f$ where \f$ bc(t) \f$
 *  is a given function of time.
 *
 */
class BoundaryCondition {
 public:
  using Indices = std::vector<unsigned int>;

  /** Basic constructor
   *
   *  \param newVelocityIndices the indices of the velocity subjected to
   * prescribed velocities
   */
  BoundaryCondition(Indices&& newVelocityIndices);
  BoundaryCondition(const Indices& newVelocityIndices);

  /** Constructor with constant prescribed values
   *
   *  \param newVelocityIndices the indices of the velocity subjected to
   * prescribed velocities \param newVelocityValues the values of the prescribed
   * velocities
   */
  BoundaryCondition(
      Indices&& newVelocityIndices,
      std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues);
  BoundaryCondition(
      const Indices& newVelocityIndices,
      std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues);

  /** destructor */
  virtual ~BoundaryCondition() noexcept = default;

  // Rule of five
  BoundaryCondition(const BoundaryCondition&) = delete;
  BoundaryCondition& operator=(const BoundaryCondition&) = delete;
  BoundaryCondition(BoundaryCondition&&) = delete;
  BoundaryCondition& operator=(BoundaryCondition&&) = delete;

  /** \return indices on which constraints are applied. */
  inline const Indices& velocityIndices() { return _velocityIndices; };

  /** \return values applied for each constrained index */
  inline auto prescribedVelocity() { return _prescribedVelocity; };

  /** \return former values applied for each constrained index */
  inline auto prescribedVelocityOld() { return _prescribedVelocityOld; };

  /** allow to set a specified function to compute prescribedVelocity
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputePrescribedVelocityFunction(const std::string& pluginPath,
                                            const std::string& functionName);

  /** default function to compute the precribed velocities
   *
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

  /** \return the number of boundary conditions
   */
  auto size() { return _velocityIndices.size(); }

  /** display */
  void display();

  /** Add a new index on which bc are to be appied
      \param ind new index value
  */
  void appendIndex(unsigned int ind);
  
 protected:
  ACCEPT_SERIALIZATION(BoundaryCondition);

  /** protected default constructor */
  BoundaryCondition() = default;

  /* Indices of the prescribed component of the velocity vector */
  Indices _velocityIndices = {};

  /* Values of the prescribed component of the velocity vector */
  std::shared_ptr<siconos::algebra::SiconosVector> _prescribedVelocity{nullptr};

  /* Old values of the prescribed component of the velocity vector */
  std::shared_ptr<siconos::algebra::SiconosVector> _prescribedVelocityOld{
      nullptr};

  /*plugin defining the function V(t)*/
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginPrescribedVelocity{
      nullptr};
};
}  // namespace siconos::modeling

#endif  // BOUNDARYCONDITION_HPP
