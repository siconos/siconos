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
/*! \file CouplerJointR.hpp
 */
#ifndef CouplerJointRELATION_H
#define CouplerJointRELATION_H

#include <NewtonEulerJointR.hpp>

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::joints {

/**
 This class implements a coupling (equality) between two
 DoFs of any NewtonEulerJointR.  Can be used e.g. to implement a
 screw relation (coupled rotation and translation) based on
 CylindricalJointR.
*/
class CouplerJointR : public NewtonEulerJointR {
 protected:
  ACCEPT_SERIALIZATION(CouplerJointR);

  std::shared_ptr<NewtonEulerJointR> joint1_{nullptr};
  std::shared_ptr<NewtonEulerJointR> joint2_{nullptr};

  std::shared_ptr<const siconos::algebra::SiconosVector> ref1_{nullptr};
  std::shared_ptr<const siconos::algebra::SiconosVector> ref2_{nullptr};

  siconos::algebra::Index dof1_{0}, dof2_{0};
  siconos::algebra::Index ref1_index_{0}, ref2_index_{0};
  double ratio_{0.};
  double offset_{0.};

  // Types to handle the different cases (q2 or not, ref or not ...)
  // when we need to call computehDoF or similar functions
  using OptRef7 = std::optional<Eigen::Map<const siconos::algebra::SiconosVector7>>;
  // std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>;
  using VArray = std::array<OptRef7, 2>;

  VArray v1_;  // for joint1
  VArray v2_;  // for joint2

  /** An internal helper function to assign reference vectors during
   * computeh and computeJachq. */
  void build_vectors(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2);

  // // void resolveVectors(const siconos::algebra::SiconosVector7* q1,
  // //                     const siconos::algebra::SiconosVector7* q2,
  // //                     const siconos::algebra::SiconosVector7*& v1_1,
  // //                     const siconos::algebra::SiconosVector7*& v1_2,
  // //                     const siconos::algebra::SiconosVector7*& v2_1,
  // //                     const siconos::algebra::SiconosVector7*& v2_2) const;

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                            const siconos::algebra::BlockVector& q0) override;

 public:
  /** Set reference joints and vectors.  This constraint maintains the
   *  equality theta2=theta1*ratio; theta1 is measured by joint1 with
   *  reference to some vector ref1 which must replace either the
   *  first or second DS of the current relation being defined.  If
   *  ref1 and ref2 are null, then no reference is used.  Typically
   *  ref1 and ref2 will be equal in order to implement gear ratios
   *  for example.  joint1 must be between ref1 and ds1 specified in
   *  setBasePositions(), while joint2 must be between ref2 and ds2.
   */

  /** Initialize a coupler*/
  CouplerJointR(std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
                std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2,
                double ratio, siconos::algebra::Index ref1_index = 0,
                std::shared_ptr<siconos::algebra::SiconosVector> ref1 = nullptr,
                siconos::algebra::Index ref2_index = 0,
                std::shared_ptr<siconos::algebra::SiconosVector> ref2 = nullptr);

  /** Initialize a coupler */
  CouplerJointR(std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
                std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2,
                double ratio,
                std::shared_ptr<siconos::modeling::NewtonEulerDS> refds1 = nullptr,
                siconos::algebra::Index ref1_index = 0,
                std::shared_ptr<siconos::modeling::NewtonEulerDS> refds2 = nullptr,
                siconos::algebra::Index ref2_index = 0);

  ~CouplerJointR() noexcept = default;

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

  /* Return the joint DoF index assigned to the first DoF. */
  siconos::algebra::Index dof1() { return dof1_; }

  /* Return the joint DoF index assigned to the second DoF. */
  siconos::algebra::Index dof2() { return dof2_; }

  /* Return the first reference joint assigned to this friction relation. */
  auto joint1() { return joint1_; }

  /* Return the second reference joint assigned to this friction relation. */
  auto joint2() { return joint2_; }

  /* Set the gear ratio. */
  void setRatio(double ratio) { ratio_ = ratio; };

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

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual siconos::algebra::Index numberOfConstraints() const override { return 1; }

  /**
     Get the number of degrees of freedom defined in the joint

     \return the number of degrees of freedom (DoF)
   */
  virtual siconos::algebra::Index numberOfDoF() const override { return 1; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(siconos::algebra::Index axis) const override {
    if (axis < 1)
      return DofType::LINEAR;
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
}  // namespace siconos::joints
#endif  // CouplerJointRELATION_H
