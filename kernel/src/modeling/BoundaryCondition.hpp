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
// #include <memory>
#include <span>
#include <vector>

#include "FunctionTypes.hpp"
// #include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"
#include "SiconosVector.hpp"
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
  using Indices = std::vector<siconos::algebra::SiconosVector::Index>;

 protected:
  ACCEPT_SERIALIZATION(BoundaryCondition);

  /** Indices of the prescribed component of the velocity vector */
  Indices velocityIndices_ = {};

  /** Values of the prescribed component of the velocity vector */
  siconos::algebra::SiconosVector prescribedVelocity_{};

  /** function wrapper used to compute external forces */
  siconos::modeling::func_prototypes::FunctionS_V computePrescribedVelocity_{nullptr};

  /** Old values of the prescribed component of the velocity vector */
  siconos::algebra::SiconosVector prescribedVelocityOld_{};

  /** protected default constructor */
  BoundaryCondition() = default;
  // Rule of five
  BoundaryCondition(const BoundaryCondition& ds) = delete;
  BoundaryCondition& operator=(const BoundaryCondition&) = delete;
  BoundaryCondition(BoundaryCondition&&) = delete;
  BoundaryCondition& operator=(BoundaryCondition&&) = delete;

 public:
  /** @brief constructor with the list of prescribed indices
   *
   *  Rq: std::move is always called for attributes init, so
   *  to optimize memory use, the best usage is
   *  BoundaryCondition bc{BoundaryCondition::Indices{v1, v2,...}}
   *  See examples and comments in BoundaryConditionTests.cpp
   *
   *  Prescribed velocity values are set to zero.
   *
   *  @param newVelocityIndices the indices of the velocity subjected to prescribed velocity
   *   values
   *
   */
  BoundaryCondition(Indices newVelocityIndices)
      : velocityIndices_{std::move(newVelocityIndices)},
        prescribedVelocity_{siconos::algebra::SiconosVector::Zero(
            siconos::algebra::to_index(velocityIndices_.size()))},
        prescribedVelocityOld_{siconos::algebra::SiconosVector::Zero(
            siconos::algebra::to_index(velocityIndices_.size()))} {};
  // Notice the std::move in constructor  + pass by value
  // BoundaryCondition{BoundaryCondition::Indices{v1, v2,...}}
  //  --> rvalue + move : only one memory print
  // BoundaryCondition{existing_indices}
  // --> lvalue + move : one copy (pass by value), as expected

  /** @brief Constructor with constant prescribed values
   *
   *  @param newVelocityIndices the indices of the velocity subjected to
   *    prescribed velocities
   *  @param newVelocityValues the values of the prescribed velocities
   */
  BoundaryCondition(
      Indices newVelocityIndices,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& newVelocityValues);

  /** destructor */
  virtual ~BoundaryCondition() noexcept = default;

  /** @return indices on which constraints are applied (read only) */
  inline std::span<const siconos::algebra::Index> velocityIndices() const {
    return velocityIndices_;
  };

  /** @return  a read-only view on the prescribed velocities vector  */
  inline auto prescribedVelocity() const {
    return siconos::algebra::ConstMapVectorType(prescribedVelocity_.data(),
                                                prescribedVelocity_.size());
  }

  /** @return  a read-only view on the former prescribed velocities vector  */
  inline auto prescribedVelocityOld() const {
    return siconos::algebra::ConstMapVectorType(prescribedVelocityOld_.data(),
                                                prescribedVelocityOld_.size());
  }

  /** set a user-defined function to compute the prescribed velocities vector
   *
   *   @param func the user-defined function (std::function, lambda ...)
   */
  void setComputePrescribedVelocityFunction(
      const siconos::modeling::func_prototypes::FunctionS_V& func);

  /** compute the prescribed velocities
   *
   *  @param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

  /** @return the number of boundary conditions */
  auto size() const { return velocityIndices_.size(); }

  /** display */
  void display() const;

  /** @brief Add a new index on which bc are to be applied
   *
   * @param ind new index value
   * @param imposedVelocity velocity value for the given index
   */
  void appendIndex(siconos::algebra::Index ind, double imposedVelocity = 0.);
};

}  // namespace siconos::modeling

#endif  // BOUNDARYCONDITION_HPP
