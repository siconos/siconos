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

  std::shared_ptr<NewtonEulerJointR> _joint{nullptr};

  std::shared_ptr<std::vector<unsigned int>> _axis{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _pos{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _dir{nullptr};

  unsigned int _axisMin{0}, _axisMax{0};
  std::shared_ptr<siconos::algebra::SiconosMatrix> jacobianhOver_q_Tmp{nullptr};

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeJacobianhOver_q_(double time, siconos::modeling::Interaction& inter,
                                       const siconos::algebra::BlockVector& q0) override;

 public:
  /** Initialize a joint stop for a common case: a single axis with a
   *  single stop, either positive or negative. For use with
   *  NewtonImpactNSL. */
  JointStopR(std::shared_ptr<NewtonEulerJointR> joint, double pos, bool dir,
             unsigned int axis = 0);

  /** Initialize a multidimensional joint stop, e.g. the cone stop on
   *  a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
  JointStopR(std::shared_ptr<NewtonEulerJointR> joint,
             std::shared_ptr<siconos::algebra::SiconosVector> pos,
             std::shared_ptr<siconos::algebra::SiconosVector> dir,
             std::shared_ptr<std::vector<unsigned int>> axes);

  /* The following constructor is disabled for now.  In fact
   * JointStopR is designed to support multiple stops, but it does not
   * work in practice since only the first value of y is examined to
   * determine whether to put the Interaction into indexSet1,
   * c.f. MoreauJeanOSI::addInteractionInIndexSet(). Use an individual
   * instance of JointStopR for each stop instead.*/
#if 0
  /** Initialize a joint stop for a common case: a single axis with a
   * double stop, one positive and one negative. */
  JointStopR(std::shared_ptr<NewtonEulerJointR> joint, double pos, double neg,
             unsigned int axis=0);
#endif
  /** destructor */
  virtual ~JointStopR() noexcept = default;

  /**
      to compute the output y = h(q) of the Relation

      \param[in] q generalized coordinates vector of the concerned dynamical systems
      \param[in,out] y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  virtual unsigned int numberOfConstraints();

  /** Return the joint axis number assigned to a stop index. */
  unsigned int axis(unsigned int _index);

  /** Return the joint position assigned to a stop index. */
  double position(unsigned int _index);

  /** Return the direction (1 or -1) assigned to a stop index. */
  double direction(unsigned int _index);

  /** Return the joint assigned to this joint stop relation. */
  std::shared_ptr<NewtonEulerJointR> joint() { return _joint; }

  /** Return the number of joint axes indexed by this relation. */
  unsigned int numberOfAxes();
};
}  // namespace siconos::joints
#endif  // JointStopRELATION_H
