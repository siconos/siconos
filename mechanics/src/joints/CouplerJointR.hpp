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

  std::shared_ptr<NewtonEulerJointR> _joint1{nullptr};
  std::shared_ptr<NewtonEulerJointR> _joint2{nullptr};

  std::shared_ptr<siconos::algebra::SiconosVector> _ref1{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _ref2{nullptr};

  unsigned int _dof1{0}, _dof2{0};
  unsigned int _ref1_index{0}, _ref2_index{0};
  double _ratio{0.}, _offset{0.};

  /** Return the normal of the linear DoF axis.  \param axis must be 0 */
  virtual void _normalDoF(siconos::algebra::SiconosVector& ans,
                          const siconos::algebra::BlockVector& q0, int axis,
                          bool absoluteRef = true) override;

  /** An internal helper function to assign reference vectors during
   * computeh and computeJachq. */
  void makeBlockVectors(std::shared_ptr<siconos::algebra::SiconosVector> q1,
                        std::shared_ptr<siconos::algebra::SiconosVector> q2,
                        siconos::algebra::BlockVector& q01,
                        siconos::algebra::BlockVector& q02);

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeJacobianhOver_q_(double time, siconos::modeling::Interaction& inter,
                                       const siconos::algebra::BlockVector& q0) override;

 public:
  /** Initialize a coupler. See setReferences() for an explanation of
   *  the parameters. */
  CouplerJointR(std::shared_ptr<NewtonEulerJointR> joint1, unsigned int dof1,
                std::shared_ptr<NewtonEulerJointR> joint2, unsigned int dof2, double ratio,
                std::shared_ptr<siconos::algebra::SiconosVector> ref1 = nullptr,
                unsigned int ref1_index = 0,
                std::shared_ptr<siconos::algebra::SiconosVector> ref2 = nullptr,
                unsigned int ref2_index = 0);

  /** Initialize a coupler. See setReferences() for an explanation of
   *  the parameters. */
  CouplerJointR(std::shared_ptr<NewtonEulerJointR> joint1, unsigned int dof1,
                std::shared_ptr<NewtonEulerJointR> joint2, unsigned int dof2, double ratio,
                std::shared_ptr<siconos::modeling::NewtonEulerDS> refds1,
                unsigned int ref1_index,
                std::shared_ptr<siconos::modeling::NewtonEulerDS> refds2,
                unsigned int ref2_index);

  ~CouplerJointR() noexcept = default;

  /**
      to compute the output y = h(q) of the Relation

      \param[in] q generalized coordinates vector of the concerned dynamical systems
      \param[in,out] y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /* Return the joint DoF index assigned to the first DoF. */
  unsigned int dof1() { return _dof1; }

  /* Return the joint DoF index assigned to the second DoF. */
  unsigned int dof2() { return _dof2; }

  /* Return the first reference joint assigned to this friction relation. */
  auto joint1() { return _joint1; }

  /* Return the second reference joint assigned to this friction relation. */
  auto joint2() { return _joint2; }

  /** Set reference joints and vectors.  This constraint maintains the
   *  equality theta2=theta1*ratio; theta1 is measured by joint1 with
   *  reference to some vector ref1 which must replace either the
   *  first or second DS of the current relation being defined.  If
   *  ref1 and ref2 are null, then no reference is used.  Typically
   *  ref1 and ref2 will be equal in order to implement gear ratios
   *  for example.  joint1 must be between ref1 and ds1 specified in
   *  setBasePositions(), while joint2 must be between ref2 and ds2.
   *
   *  \param joint1 The joint for the first reference measurement theta1.
   *  \param dof1 The degree of freedom index of joint1 to use for
   *              the first reference measurement theta1.
   *  \param ref1 The optional reference position for the first
   *              reference measurement theta1.
   *  \param ref1_index Must be 0 or 1, depending on where ref1
   *                    appears in joint1.
   *  \param joint2 The joint for the second reference measurement theta2.
   *  \param dof2 The degree of freedom index of joint2 to use for
   *              the second reference measurement theta2.
   *  \param ref2 The optional reference position for the second
   *              reference measurement theta2.
   *  \param ref2_index Must be 0 or 1, depending on where ref2
   *                    appears in joint2.
   */
  void setReferences(std::shared_ptr<NewtonEulerJointR> joint1, unsigned int dof1,
                     std::shared_ptr<NewtonEulerJointR> joint2, unsigned int dof2,
                     std::shared_ptr<siconos::algebra::SiconosVector> ref1,
                     unsigned int ref1_index,
                     std::shared_ptr<siconos::algebra::SiconosVector> ref2,
                     unsigned int ref2_index);

  /** Set reference joints and vectors.  This constraint maintains the
   *  equality theta2=theta1*ratio; theta1 is measured by joint1 with
   *  reference to some vector ref1 which must replace either the
   *  first or second DS of the current relation being defined.  If
   *  ref1 and ref2 are null, then no reference is used.  Typically
   *  ref1 and ref2 will be equal in order to implement gear ratios
   *  for example.  joint1 must be between ref1 and ds1 specified in
   *  setBasePositions(), while joint2 must be between ref2 and ds2.
   *  This alternative setReferences() takes NewtonEulerDS as
   *  parameters, but the reference vectors will be taken as
   *  refds1->q() and refds2->q() respectively.
   *
   *  \param joint1 The joint for the first reference measurement theta1.
   *  \param dof1 The degree of freedom index of joint1 to use for
   *              the first reference measurement theta1.
   *  \param refds1 The optional reference vector for the first
   *                reference measurement theta1.
   *  \param ref1_index Must be 0 or 1, depending on where ref1
   *                    appears in joint1.
   *  \param joint2 The joint for the second reference measurement theta2.
   *  \param dof2 The degree of freedom index of joint2 to use for
   *              the second reference measurement theta2.
   *  \param refds2 The optional reference vector for the second
   *                reference measurement theta2.
   *  \param ref2_index Must be 0 or 1, depending on where ref2
   *                    appears in joint2.
   */
  void setReferences(std::shared_ptr<NewtonEulerJointR> joint1, unsigned int dof1,
                     std::shared_ptr<NewtonEulerJointR> joint2, unsigned int dof2,
                     std::shared_ptr<siconos::modeling::NewtonEulerDS> refds1,
                     unsigned int ref1_index,
                     std::shared_ptr<siconos::modeling::NewtonEulerDS> refds2,
                     unsigned int ref2_index);

  /* Set the gear ratio. */
  void setRatio(double ratio);

  /* From provided position vectors, calculate initial offsets to be
   * maintained in the relation. */
  void setBasePositions(std::shared_ptr<siconos::algebra::SiconosVector> q1,
                        std::shared_ptr<siconos::algebra::SiconosVector> q2 =
                            std::shared_ptr<siconos::algebra::SiconosVector>()) override;

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() override { return 1; }

  /**
     Get the number of degrees of freedom defined in the joint

     \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() override { return 1; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(unsigned int axis) override {
    if (axis < 1)
      return DofType::LINEAR;
    else
      return DofType::INVALID;
  };

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(const siconos::algebra::BlockVector& q0,
                           Eigen::Ref<siconos::algebra::SiconosVector> y,
                           unsigned int axis = 0) override;

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(siconos::modeling::Interaction& inter,
                               const siconos::algebra::BlockVector& q0,
                               Eigen::Ref<siconos::algebra::SiconosMatrix> jachq,
                               unsigned int axis = 0) override;
};
}  // namespace siconos::joints
#endif  // CouplerJointRELATION_H
