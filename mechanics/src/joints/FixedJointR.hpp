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
/** \file FixedJointR.hpp
 */
#ifndef FixedJointRELATION_H
#define FixedJointRELATION_H

#include <NewtonEulerJointR.hpp>

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::joints {

/**
   This class implements a fixed joint between one or two Newton/Euler Dynamical system
*/
class FixedJointR : public NewtonEulerJointR {
 protected:
  ACCEPT_SERIALIZATION(FixedJointR);

  /*Initial conditions*/
  double _G10G20d1x{0.}, _G10G20d1y{0.}, _G10G20d1z{0.};
  double _cq2q101{0.}, _cq2q102{0.}, _cq2q103{0.}, _cq2q104{0.};

 public:
  /** constructor,
   *
   *  \param a std::shared_ptr<siconos::modeling::NewtonEulerDS> d1, a dynamical system
   * containing the initial position \param a std::shared_ptr<siconos::modeling::NewtonEulerDS>
   * d2, a dynamical system containing the initial position
   */
  FixedJointR(std::shared_ptr<siconos::modeling::NewtonEulerDS> d1,
              std::shared_ptr<siconos::modeling::NewtonEulerDS> d2 = nullptr);

  /** destructor
   */
  virtual ~FixedJointR() noexcept = default;

  /** Initialize the joint constants based on the provided base positions.
   *
   *  \param q1 A siconos::algebra::SiconosVector of size 7 indicating translation and
   *  orientation in inertial coordinates.
   *  \param q2 An optional siconos::algebra::SiconosVector of size 7 indicating
   *  translation and orientation; if null, the inertial
   *  frame will be considered as the second base. */
  virtual void setBasePositions(
      std::shared_ptr<siconos::algebra::SiconosVector> q1,
      std::shared_ptr<siconos::algebra::SiconosVector> q2 = nullptr) override;

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() override { return 6; }

  virtual void computeJachq(double time, siconos::modeling::Interaction& inter,
                            std::shared_ptr<siconos::algebra::BlockVector> q0) override;

  /**
     to compute the output y = h(t,q,z) of the Relation

     \param time current time value
     \param q coordinates of the dynamical systems involved in the relation
     \param y the resulting vector
  */
  virtual void computeh(double time, const siconos::algebra::BlockVector& q0,
                        siconos::algebra::SiconosVector& y) override;

  virtual unsigned int numberOfDoF() override { return 0; }

  virtual DofType typeOfDoF(unsigned int axis) override { return DofType::INVALID; }

 protected:
  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12,
                     double q13, double X2, double Y2, double Z2, double q20, double q21,
                     double q22, double q23);

  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12,
                   double q13);
};
}  // namespace siconos::joints
#endif  // FixedJointRELATION_H
