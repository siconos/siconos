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
/*! \file NewtonEulerJointR.hpp

*/
#ifndef NewtonEulerJointRELATION_H
#define NewtonEulerJointRELATION_H

#include <optional>
#include <vector>

#include "NewtonEulerR.hpp"
namespace siconos::joints {

enum class DofType {
  INVALID = 0,
  LINEAR = 1,
  ANGULAR = 2,
};

/**
    This class implements an abstract Joint relation (articulation) between one or two
   Newton/Euler dynamical systems.
*/
class NewtonEulerJointR : public siconos::modeling::NewtonEulerR {
 protected:
  ACCEPT_SERIALIZATION(NewtonEulerJointR);

  /** A flag determining whether this joint should block
   *  "self-collision", i.e., if true, bodies connected by this joint
   *  will not enter into unilateral contact. */
  bool allowSelfCollide_{false};

  /** Vector of points used to defined the joint constraint. Default size = 0. */
  std::vector<siconos::algebra::SiconosVector3> points_ = {};

  /** Vector of axes used to defined the joint constraint. Default size = 0 */
  std::vector<siconos::algebra::SiconosVector3> axes_ = {};

  /** Defines whether points and axes are specified in absolute or
   * relative frame. */
  bool absoluteRef_{true};

 public:
  // No constructors : default is implemented in base class NewtonEulerR

  virtual ~NewtonEulerJointR() noexcept = default;

  /** Set a point for this joint. The role of each point is specific
   *  to the joint subclass. Won't take effect until
   *  setBasePositions is called.
   *
   *  \param index The index of the point to be set
   *  \param point The point coordinates as a vector of size 3
   */
  void setPoint(size_t index, const Eigen::Ref<siconos::algebra::SiconosVector3>& point);

  /** Set an axis for this joint. The role of each axis is specific to
   *  the joint subclass. Won't take effect until setBasePositions
   *  is called.
   *
   *  \param index The index of the axis to be set
   *  \param axis The axis coordinates as a vector of size 3
   */
  void setAxis(size_t index, const Eigen::Ref<siconos::algebra::SiconosVector3>& axis);

  /** \return a read-only view on the axis at position index
   *
   *  \param index required position in axes vector
   */
  inline auto axis(unsigned int index) {
    return siconos::algebra::ConstMapVectorType(axes_[index].data(), axes_[index].size());
  }

  /** Set whether points and axes should be interpreted in absolute or
   *  relative frame. Won't take effect until setBasePositions is
   *  called.
   *
   *  \param absoluteRef true for absolute frame, false for relative frame.
   */
  void setAbsolute(bool absoluteRef) { absoluteRef_ = absoluteRef; }

  /** Get whether points and axes are interpreted in absolute or
   *  relative frame.
   *
   *  \return True for absolute frame, false for relative frame.
   */
  bool absolute() { return absoluteRef_; }

  /** Initialize the joint constants based on the provided base positions.
   *
   *  \param[in] q1 a vector of size 7 indicating translation and orientation in inertial
   * coordinates.
   *  \param[in] q2 an optional vector of size 7 indicating translation and orientation; if
   * null, the inertial frame will be considered as the second base.
   */
  virtual void setBasePositions(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2 =
          std::nullopt) = 0;


  /** \return the projection of a vector onto the given 0-indexed free axis. Useful for
   *  calculating velocities in the axis, or for calculating
   *  axis-aligned forces applied to connected bodies.  If axis is of
   *  angular type (see typeOfDoF), then the projection is onto the
   *  axis of rotation.
   *
   *  \param[in] v The vector to project
   *  \param[in] q0 The state q of the first NewtonEulerDS
   *  \param[in] q1 The state q of the second NewtonEulerDS (optional)
   *  \param[in] absoluteRef If true, v and the result are in the inertial frame,
   *  otherwise the q1 frame is assumed.
   */
  siconos::algebra::SiconosVector3 projectVectorDoF(const siconos::algebra::SiconosVector& v,
                                                    const siconos::algebra::SiconosVector& q0,
                                                    const std::optional<Eigen::Ref<siconos::algebra::SiconosVector>>& q1,
                                                    int axis = 0, bool absoluteRef = true);

  

  /** \return the axis of rotation.
   *  Retrieve a normal in the direction of a 0-indexed free
   *  axis. Useful for calculating velocities in the axis, or for
   *  calculating axis-aligned forces applied to connected bodies.  If
   *  axis is of angular type (see typeOfDoF), then the returned normal
   *  is the axis of rotation.
   *
   *  \param[in] q0 The state q of one NewtonEulerDS
   *  \param[in] q1 The state q of the second NewtonEulerDS (optional)
   *  \param[in] axis
   *  \param[in] absoluteRef If true, ans is in the inertial frame,
   *  otherwise the q1 frame is assumed.
   */
  virtual siconos::algebra::SiconosVector3 normalDoF(
      const siconos::algebra::SiconosVector& q0,
      const std::optional<Eigen::Ref<siconos::algebra::SiconosVector>>& q1 = std::nullopt, int axis = 0,
      bool absoluteRef = true) {
    throw std::logic_error("normalDof  not implemented for this kind of joint");
  }

  /** Return the value of the allowSelfCollide_ flag. */
  bool allowSelfCollide() { return allowSelfCollide_; }

  /** Set the value of the allowSelfCollide_ flag. */
  void setAllowSelfCollide(bool x) { allowSelfCollide_ = x; }

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() const = 0;

  /**
     Return the number of degrees of freedom of this joint.

     \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() const = 0;

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
   */
  virtual DofType typeOfDoF(unsigned int axis) const = 0;

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y, unsigned int axis = 0) {}

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(siconos::modeling::Interaction& inter,
                               const siconos::algebra::BlockVector& q0,
                               Eigen::Ref<siconos::algebra::SiconosMatrix> jachq,
                               unsigned int axis = 0) {}
  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::joints
#endif  // NewtonEulerJointRELATION_H
