/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file JointFrictionR.hpp

*/
#ifndef JointFrictionRELATION_H
#define JointFrictionRELATION_H

#include <vector>

#include "NewtonEulerR.hpp"

namespace siconos::joints {

class NewtonEulerJointR;

/**
   This class implements a friction on a DoF for any NewtonEulerJointR.
*/
class JointFrictionR : public siconos::modeling::NewtonEulerR {
 protected:
  ACCEPT_SERIALIZATION(JointFrictionR);

  std::shared_ptr<NewtonEulerJointR> _joint{nullptr};

  std::shared_ptr<std::vector<siconos::algebra::Index>> _axis{nullptr};

  siconos::algebra::Index _axisMin{0}, _axisMax{0};
  std::shared_ptr<siconos::algebra::SiconosDenseMatrix> jacobianhOver_q_Tmp{nullptr};

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeConstraintJacobian(double time, siconos::modeling::Interaction& inter,
                                         const siconos::algebra::BlockVector& q0) override;

 public:
  /** Initialize a joint friction for a common case: a single axis with a
   *  single friction, either positive or negative. For use with
   *  NewtonImpactNSL. */
  JointFrictionR(std::shared_ptr<NewtonEulerJointR> joint, siconos::algebra::Index axis);

  /** Initialize a multidimensional joint friction, e.g. the cone friction on
   *  a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
  JointFrictionR(std::shared_ptr<NewtonEulerJointR> joint,
                 std::shared_ptr<std::vector<siconos::algebra::Index>> axes = nullptr);

  virtual ~JointFrictionR() noexcept = default;

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

  virtual size_t numberOfConstraints() const;

  /** Return the joint axis number assigned to a friction axis. */
  siconos::algebra::Index axis(size_t _index);

  /** Return the joint assigned to this friction relation. */
  std::shared_ptr<NewtonEulerJointR> joint() { return _joint; }

  size_t numberOfAxes();

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::joints
#endif  // JointFrictionRELATION_H
