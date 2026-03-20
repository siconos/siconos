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
/*! \file PivotJointR.hpp
 */
#ifndef PivotJointRELATION_H
#define PivotJointRELATION_H

#include "KneeJointR.hpp"

namespace siconos::joints {
/**
   This class implements a pivots joint between one or two Newton/Euler Dynamical system. -
   Inherits from KneeJointR
*/
class PivotJointR : public KneeJointR {
 protected:
  ACCEPT_SERIALIZATION(PivotJointR);

  /*Axis coordonates*/
  std::shared_ptr<siconos::algebra::SiconosVector3> axis_coords_{nullptr};

  siconos::algebra::SiconosVector3 axis1_{siconos::algebra::SiconosVector3::Zero()};
  siconos::algebra::SiconosVector3 axis2_{siconos::algebra::SiconosVector3::Zero()};

  /*Initial conditions*/
  double _initial_AscalA{0.}, _initial_AscalA1{0.}, _initial_AscalA2{0.};
  siconos::algebra::SiconosVector4 cq2q_ = siconos::algebra::SiconosVector4::Zero();

  /** Cumulative number of twists around the joint relative to initial
   * angular difference. */
  int _twistCount{0};        // TODO: Should be in a graph work vector?
  double _previousAngle{0};  // Needed to track _twistCount, TODO: work vector?

  void buildOrthonormalBase();

  /** compute H_NE matrix
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                            const siconos::algebra::BlockVector& q0) override;

 public:
  /** Default constructor */
  PivotJointR();

  /** Constructor based on one or two dynamical systems, a point and an axis.
   *
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, (null for absolute frame)
   *  \param P a vector that defines the point around which rotation is allowed.
   *  \param A a vector that defines the cylindrical axis.
   *  \param absoluteRef if true, P and A are in the absolute frame,g
   *                     otherwise P and A are in d1 frame.
   */
  PivotJointR(const Eigen::Ref<siconos::algebra::SiconosVector3>& P,
              const Eigen::Ref<siconos::algebra::SiconosVector3>& A, bool absoluteRef,
              std::shared_ptr<siconos::modeling::NewtonEulerDS> d1 = nullptr,
              std::shared_ptr<siconos::modeling::NewtonEulerDS> d2 = nullptr);

  virtual ~PivotJointR() noexcept = default;

  /** Initialize the joint constants based on the provided base positions.
   *
   *  \param[in] q1 a vector of size 7 indicating translation and orientation in inertial
   * coordinates.
   *  \param[in] q2 an optional vector of size 7 indicating translation and orientation; if
   * null, the inertial frame will be considered as the second base.
   */
  virtual void setBasePositions(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2 =
          std::nullopt) override;

  /** \return the axis of rotation.
   *  Retrieve a normal in the direction of a 0-indexed free
   *  axis. Useful for calculating velocities in the axis, or for
   *  calculating axis-aligned forces applied to connected bodies.  If
   *  axis is of angular type (see typeOfDoF), then the returned normal
   *  is the axis of rotation.
   *
   *  \param[in] q0 The state q of the fiest NewtonEulerDS
   *  \param[in] q1 The state q of the second NewtonEulerDS
   *  \param[in] axis
   *  \param[in] absoluteRef If true, ans is in the inertial frame,
   *  otherwise the q1 frame is assumed.
   */
  virtual siconos::algebra::SiconosVector3 normalDoF(
      const siconos::algebra::SiconosVector7& q0,
      const std::optional<Eigen::Ref<siconos::algebra::SiconosVector7>>& q1 = std::nullopt,
      int axis = 0, bool absoluteRef = true) override;

  /**
     to compute the output y = h(q) of the Relation

      \param[in] q1 generalized coordinates vector of the fist dynamical system involved
      in the relation
      \param[in] q2 generalized coordinates vector of the second dynamical system
      involved in the relation
      \param[in,out] y the resulting vector
  */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y) override;
  /**
       Get the number of constraints defined in the joint

       \return the number of constraints
     */
  virtual siconos::algebra::Index numberOfConstraints() const override { return 5; }

  /**
     Return the number of degrees of freedom of this joint.

     \return the number of degrees of freedom (DoF)
   */
  virtual siconos::algebra::Index numberOfDoF() const override { return 1; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(siconos::algebra::Index axis) const override {
    if (axis == 0)
      return DofType::ANGULAR;
    else
      return DofType::INVALID;
  };
  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y,
      siconos::algebra::Index axis = 0) override;

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(
      siconos::modeling::Interaction& inter,
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosDenseMatrix> jachq,
      siconos::algebra::Index axis = 0) override;
  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};

// Free functions for pivot joints
namespace pivot {
void computeH_for_2DS(const Eigen::Ref<const siconos::algebra::SiconosVector4>& qp1,
                      const siconos::algebra::SiconosVector3& coords1,
                      const Eigen::Ref<const siconos::algebra::SiconosVector4>& qp2,
                      const siconos::algebra::SiconosVector3& coords2,
                      const siconos::algebra::SiconosVector3& A1,
                      const siconos::algebra::SiconosVector3& A2,
                      const siconos::algebra::SiconosVector4& c2q,
                      Eigen::Ref<siconos::algebra::MapType> result);

void computeH_for_1DS(const Eigen::Ref<const siconos::algebra::SiconosVector4>& qp1,
                      const siconos::algebra::SiconosVector3& coords1,
                      const siconos::algebra::SiconosVector3& A1,
                      const siconos::algebra::SiconosVector3& A2,
                      const siconos::algebra::SiconosVector4& c2q,
                      Eigen::Ref<siconos::algebra::MapType> result);

void rot2to1(const Eigen::Ref<const siconos::algebra::SiconosVector4>& qp1,
             const Eigen::Ref<const siconos::algebra::SiconosVector4>& qp2,
             const siconos::algebra::SiconosVector4& cq2q,
             Eigen::Ref<siconos::algebra::SiconosVector4> result);

}  // namespace pivot

}  // namespace siconos::joints
#endif  // PivotJointRELATION_H
