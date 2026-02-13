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
/** \file KneeJointR.hpp
 */
#ifndef KneeJointRELATION_H
#define KneeJointRELATION_H

#include <NewtonEulerJointR.hpp>

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::joints {

/**
   This class implements a knee joint between one or two Newton/Euler Dynamical system
*/
class KneeJointR : public NewtonEulerJointR {
 protected:
  ACCEPT_SERIALIZATION(KneeJointR);

  // For this joint, points_[0] corresponds to the coordinate of the knee point
  // in the body frame of the first dynamical system _d1

  /**Absolute coodinates of the vector  G1P0 when d1 is located in q=(0,0,0,1,0,0,0)
   * i.e. points_[0] in the body frame of d1.
   * These values are computed when the constructor is called.
   */
  siconos::algebra::SiconosVector3 G1P0_ = siconos::algebra::SiconosVector3::Zero();

  /** Absolute coodinates of the vector G2P0 when d2 is located in q=(0,0,0,1,0,0,0)
   *  i.e. points_[0] in the body frame of d2.
   * These values are computed when the constructor is called.
   */
  siconos::algebra::SiconosVector3 G2P0_ = siconos::algebra::SiconosVector3::Zero();

  /** compute H_NE matrix
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                            const siconos::algebra::BlockVector& q0) override;

 public:
  /** Default constructor
   */
  KneeJointR();

  /** Constructor based on one or two dynamical systems and a point.
   *
   *  \param P SiconosVector of size 3 that defines the point around
   *   which rotation is allowed.
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, or NULL
   *  for absolute frame.
   *  \param absoluteRef if true, P is in the absolute frame,
   *   otherwise P is in d1 frame.
   */
  KneeJointR(const Eigen::Ref<siconos::algebra::SiconosVector3>& P, bool absoluteRef,
             std::shared_ptr<siconos::modeling::NewtonEulerDS> d1 = nullptr,
             std::shared_ptr<siconos::modeling::NewtonEulerDS> d2 = nullptr);

  /** destructor */
  virtual ~KneeJointR() noexcept = default;

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
          std::nullopt) override;

  /** Perform some checks on the initial conditions. */
  void checkInitPos(const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
                    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>&
                        q2 = std::nullopt);

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual siconos::algebra::Index numberOfConstraints() const override { return 3; }

  /**
     Get the number of degrees of freedom defined in the joint

     \return the number of degrees of freedom (DoF)
   */
  virtual siconos::algebra::Index numberOfDoF() const override { return 3; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(unsigned int axis) const override {
    if (axis < 3)
      return DofType::ANGULAR;
    else
      return DofType::INVALID;
  };

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
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};

// Free functions for knee joints

namespace knee {
void hfunction(const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
               const siconos::algebra::SiconosVector3& coords1,
               const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2opt,
               const siconos::algebra::SiconosVector3& coords2,
               Eigen::Ref<siconos::algebra::SiconosVector> result);

void computeH_for_2DS(const Eigen::Ref<const siconos::algebra::SiconosVector>& qp1,
                      const siconos::algebra::SiconosVector3& coords1,
                      const Eigen::Ref<const siconos::algebra::SiconosVector>& qp2,
                      const siconos::algebra::SiconosVector3& coords2,
                      Eigen::Ref<siconos::algebra::MapType> result);

void computeH_for_1DS(const Eigen::Ref<const siconos::algebra::SiconosVector>& qp1,
                      const siconos::algebra::SiconosVector3& coords1,
                      Eigen::Ref<siconos::algebra::MapType> result);

void computeH_dot_for1DS(const Eigen::Ref<const siconos::algebra::SiconosVector>& qpdot1,
                         const siconos::algebra::SiconosVector3& coords1,
                         Eigen::Ref<siconos::algebra::MapType> result);

void computeH_dot_for2DS(const Eigen::Ref<const siconos::algebra::SiconosVector>& qpdot1,
                         const siconos::algebra::SiconosVector3& coords1,
                         const Eigen::Ref<const siconos::algebra::SiconosVector>& qpdot2,
                         const siconos::algebra::SiconosVector3& coords2,
                         Eigen::Ref<siconos::algebra::MapType> result);

}  // namespace knee

}  // namespace siconos::joints
#endif  // KneeJointRELATION_H
