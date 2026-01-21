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

  std::shared_ptr<std::vector<unsigned int>> _axis{nullptr};

  unsigned int _axisMin{0}, _axisMax{0};
  std::shared_ptr<siconos::algebra::SiconosMatrix> jacobianhOver_q_Tmp{nullptr};

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                            const siconos::algebra::BlockVector& q0) override;

 public:
  /** Initialize a joint friction for a common case: a single axis with a
   *  single friction, either positive or negative. For use with
   *  NewtonImpactNSL. */
  JointFrictionR(std::shared_ptr<NewtonEulerJointR> joint, unsigned int axis);

  /** Initialize a multidimensional joint friction, e.g. the cone friction on
   *  a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
  JointFrictionR(std::shared_ptr<NewtonEulerJointR> joint,
                 std::shared_ptr<std::vector<unsigned int>> axes = nullptr);

  virtual ~JointFrictionR() noexcept = default;

  /**
     to compute the output y = h(q) of the Relation

     \param[in] q generalized coordinates vector of the dynamical systems (at most 2) involved
    in the relation \param[in,out] y the resulting vector
 */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /**
    to compute the output y = h(q) of the Relation

   \param[in] q1 generalized coordinates vector of the first dynamical system involved
   in the relation
   \param[in] q2 generalized coordinates vector of the second dynamical system
   involved in the relation
   \param[in,out] y the resulting vector
 */
  void computeh(const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
                const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  virtual unsigned int numberOfConstraints() const;

  /** Return the joint axis number assigned to a friction axis. */
  unsigned int axis(unsigned int _index);

  /** Return the joint assigned to this friction relation. */
  std::shared_ptr<NewtonEulerJointR> joint() { return _joint; }

  unsigned int numberOfAxes();

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::joints
#endif  // JointFrictionRELATION_H
