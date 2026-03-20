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
/*! \file JointStopR.hpp

*/
#ifndef JointStopRELATION_H
#define JointStopRELATION_H

#include <vector>

#include "NewtonEulerR.hpp"

namespace siconos::joints {

class NewtonEulerJointR;

/**
   This class implements a stop on a DoF for any NewtonEulerJointR.
*/
class JointStopR : public siconos::modeling::NewtonEulerR {
 protected:
  ACCEPT_SERIALIZATION(JointStopR);

  std::shared_ptr<NewtonEulerJointR> joint_{nullptr};

  std::vector<siconos::algebra::Index> axis_ = {};
  siconos::algebra::SiconosVector pos_;
  siconos::algebra::SiconosVector dir_;

  siconos::algebra::Index axisMin_{0}, axisMax_{0};
  siconos::algebra::SiconosMatrix jacobian_buffer_;

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                            const siconos::algebra::BlockVector& q0) override;

 public:
  /** Initialize a joint stop for a common case: a single axis with a
   *  single stop, either positive or negative. For use with
   *  NewtonImpactNSL. */
  JointStopR(std::shared_ptr<NewtonEulerJointR> joint, double pos, bool dir,
             siconos::algebra::Index axis = 0);

  /** Initialize a multidimensional joint stop, e.g. the cone stop on
   *  a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
  JointStopR(std::shared_ptr<NewtonEulerJointR> joint,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& pos,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& dir,
             const std::vector<siconos::algebra::Index>& axes);

#if 0
  /* The following constructor is disabled for now.  In fact
   * JointStopR is designed to support multiple stops, but it does not
   * work in practice since only the first value of y is examined to
   * determine whether to put the Interaction into indexSet1,
   * c.f. MoreauJeanOSI::addInteractionInIndexSet(). Use an individual
   * instance of JointStopR for each stop instead.*/

  /** Initialize a joint stop for a common case: a single axis with a
   * double stop, one positive and one negative. */
  JointStopR(std::shared_ptr<NewtonEulerJointR> joint, double pos, double neg,
             siconos::algebra::Index axis=0);
#endif
  /** destructor */
  virtual ~JointStopR() noexcept = default;

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

  virtual std::size_t numberOfConstraints() const;

  /** Return the joint axis number assigned to a stop index. */
  auto axis(size_t _index);

  /** Return the joint position assigned to a stop index. */
  double position(siconos::algebra::Index _index);

  /** Return the direction (1 or -1) assigned to a stop index. */
  double direction(siconos::algebra::Index _index);

  /** Return the joint assigned to this joint stop relation. */
  std::shared_ptr<NewtonEulerJointR> joint() { return joint_; }

  /** Return the number of joint axes indexed by this relation. */
  auto numberOfAxes();

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::joints
#endif  // JointStopRELATION_H
