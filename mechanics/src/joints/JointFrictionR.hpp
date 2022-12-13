/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "NewtonEulerR.hpp"
#include <vector>

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
  std::shared_ptr<siconos::algebra::SimpleMatrix> _jachqTmp{nullptr};

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
     to compute the output y = h(t,q,z) of the Relation

     \param time current time value
     \param q coordinates of the dynamical systems involved in the relation
     \param y the resulting vector
  */
  virtual void computeh(double time, const siconos::algebra::BlockVector& q0,
                        siconos::algebra::SiconosVector& y) override;

  virtual void computeJachq(double time, siconos::modeling::Interaction& inter,
                            std::shared_ptr<siconos::algebra::BlockVector> q0) override;

  virtual unsigned int numberOfConstraints();

  /** Return the joint axis number assigned to a friction axis. */
  unsigned int axis(unsigned int _index);

  /** Return the joint assigned to this friction relation. */
  std::shared_ptr<NewtonEulerJointR> joint() { return _joint; }

  unsigned int numberOfAxes();
};
}  // namespace siconos::joints
#endif  // JointFrictionRELATION_H
