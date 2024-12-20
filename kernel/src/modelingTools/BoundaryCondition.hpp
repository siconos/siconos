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

#include "FunctionTypes.hpp"
#include "SiconosSerialization.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {
// class SiconosVector;
}

namespace siconos::modeling {

/** This class models simple boundary conditions for
 *  prescribing the velocities in a Dynamical System.
 *
 * A simple boundary condition is considered to fix a component
 *  \f$ j \f$ of the velocity vector, i.e.,
 *
 * \f$ v_j(t) = bc(t) \f$
 *
 *  where \f$ bc(t) \f$  is a given function of time.
 *
 */
class BoundaryCondition {
 public:
  /** Type of list of indices used in boundary conditions */
  using Indices = std::vector<unsigned int>;

 protected:
  ACCEPT_SERIALIZATION(BoundaryCondition);

  /** Indices of the prescribed component of the velocity vector */
  Indices velocityIndices_ = {};

  /** Values of the prescribed component of the velocity vector */
  std::shared_ptr<siconos::algebra::SiconosVector> prescribedVelocity_{nullptr};

  /** function wrapper used to compute external forces */
  siconos::modeling::func_prototypes::FunctionS_V computePrescribedVelocity_{nullptr};

  /** Old values of the prescribed component of the velocity vector */
  std::shared_ptr<siconos::algebra::SiconosVector> prescribedVelocityOld_{nullptr};

  /** protected default constructor */
  BoundaryCondition() = default;
  // Rule of five
  BoundaryCondition(const BoundaryCondition& ds) = delete;
  BoundaryCondition& operator=(const BoundaryCondition&) = delete;
  BoundaryCondition(BoundaryCondition&&) = delete;
  BoundaryCondition& operator=(BoundaryCondition&&) = delete;

 public:
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
  BoundaryCondition(Indices&& newVelocityIndices,
                    std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues);
  BoundaryCondition(const Indices& newVelocityIndices,
                    std::shared_ptr<siconos::algebra::SiconosVector> newVelocityValues);

  /** destructor */
  virtual ~BoundaryCondition() noexcept = default;

  /** \return indices on which constraints are applied. */
  inline const Indices& velocityIndices() { return velocityIndices_; };

  /** \return  a read-only view on the prescribed velocities vector  */
  inline const auto prescribedVelocity() const {
    return siconos::algebra::ConstMapVectorType(prescribedVelocity_->data(),
                                                prescribedVelocity_->size());
  }

  /** \return  a read-only view on the former prescribed velocities vector  */
  inline const auto prescribedVelocityOld() const {
    return siconos::algebra::ConstMapVectorType(prescribedVelocityOld_->data(),
                                                prescribedVelocityOld_->size());
  }

  /** set a user-defined function to compute the prescribed velocities vector
   *
   *   \param func the user-defined function (std::function, lambda ...)
   */
  void setComputePrescribedVelocityFunction(
      const siconos::modeling::func_prototypes::FunctionS_V& func);

  /** compute the prescribed velocities
   *
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

  /** \return the number of boundary conditions */
  auto size() { return velocityIndices_.size(); }

  /** display */
  void display() const;

  /** Add a new index on which bc are to be appied
   *
   * \param ind new index value
   */
  void appendIndex(unsigned int ind);
};
}  // namespace siconos::modeling

#endif  // BOUNDARYCONDITION_HPP
