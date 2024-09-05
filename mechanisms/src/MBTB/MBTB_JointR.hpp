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

#ifndef MBTB_JOINTR
#define MBTB_JOINTR

#include <memory>

namespace siconos::algebra {
class SiconosVector;
class SimpleMatrix;
}  // namespace siconos::algebra

namespace siconos::modeling {
class NewtonEulerDS;
class Interaction;
}  // namespace siconos::modeling

namespace siconos::joints {
class PivotJointR;
class NewtonEulerJointR;
}  // namespace siconos::joints

namespace siconos::mechanisms {
/**
 * \brief This class implements a joint in a multi-bodies system.
 * It is an aggregation to the class siconos::NewtonEulerR. Mainly, it consists
 * in adding members needed for the computation of joint forces.
 */
class MBTB_JointR {
  friend class siconos::joints::PivotJointR;

 protected:
 public:
  MBTB_JointR();
  //! it is assumed that _interaction
  std::shared_ptr<siconos::modeling::Interaction> _interaction{nullptr};
  //! it is assumed that _joinrR is a pivot.
  std::shared_ptr<siconos::joints::NewtonEulerJointR> _jointR{nullptr};
  //! The first dynamical systems of the joint.
  std::shared_ptr<siconos::modeling::NewtonEulerDS> _ds1{nullptr};

  //! vector G0Ci, where Ci Contact points where the forces of joint must be
  //! computed.
  std::shared_ptr<siconos::algebra::SiconosVector> _G0C1{nullptr};
  //! vector G0Ci, where Ci Contact points where the forces of joint must be
  //! computed.
  std::shared_ptr<siconos::algebra::SiconosVector> _G0C2{nullptr};
  //! Joint forces. F1: F[0,1,2]. F2: F[3,4,5].
  std::shared_ptr<siconos::algebra::SiconosVector> _F{nullptr};
  //! A matrix such that  _M * F = BLambda
  std::shared_ptr<siconos::algebra::SimpleMatrix> _M{nullptr};

  //! It consists in building the system  _M * F = BLambda  and solving it.
  void computeEquivalentForces();
};
}  // namespace siconos::mechanisms
#endif
