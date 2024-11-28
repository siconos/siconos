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
/*! \file PivotJointR.hpp
 */
#ifndef PivotJointRELATION_H
#define PivotJointRELATION_H

#include "KneeJointR.hpp"

namespace siconos::joints {
/**
   This class implements a pivots joint between one or two Newton/Euler Dynamical system. -
   Inherits from KneeJointR
*/
class PivotJointR : public KneeJointR {
 protected:
  ACCEPT_SERIALIZATION(PivotJointR);

  /*Axis coordonates*/
  std::shared_ptr<siconos::algebra::SiconosVector> _A{nullptr};
  double _A1x{0.}, _A1y{0.}, _A1z{0.};
  double _A2x{0.}, _A2y{0.}, _A2z{0.};

  /*Initial conditions*/
  double _cq2q101{0.}, _cq2q102{0.}, _cq2q103{0.}, _cq2q104{0.};
  double _initial_AscalA{0.}, _initial_AscalA1{0.}, _initial_AscalA2{0.};

  /** Cumulative number of twists around the joint relative to initial
   * angular difference. */
  int _twistCount{0};     // TODO: Should be in a graph work vector?
  int _previousAngle{0};  // Needed to track _twistCount, TODO: work vector?

  void buildA1A2();

  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12,
                     double q13, double X2, double Y2, double Z2, double q20, double q21,
                     double q22, double q23) override;

  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12,
                   double q13) override;

  void rot2to1(double q10, double q11, double q12, double q13, double q20, double q21,
               double q22, double q23, double* q2to1w, double* q2to1x, double* q2to1y,
               double* q2to1z);

  double AscalA1(double q2to1x, double q2to1y, double q2to1z);
  double AscalA2(double q2to1x, double q2to1y, double q2to1z);
  double AscalA(double q2to1x, double q2to1y, double q2to1z);

  /** Return the normal of the angular DoF axis of rotation.
   *
   *  \param axis must be 0 */
  virtual void _normalDoF(siconos::algebra::SiconosVector& ans,
                          const siconos::algebra::BlockVector& q0, int axis,
                          bool absoluteRef = true) override;

 public:
  /** Empty constructor. The relation may be initialized later by
   * setPoint, setAxis, setAbsolute, and setBasePositions. */
  PivotJointR();

  /** Constructor based on one or two dynamical systems, a point and an axis.
   *
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, or NULL
   *            for absolute frame.
   *  \param P siconos::algebra::SiconosVector of size 3 that defines the point around
   *           which rotation is allowed.
   *  \param A siconos::algebra::SiconosVector of size 3 that defines the cylindrical axis.
   *  \param absoluteRef if true, P and A are in the absolute frame,
   *                     otherwise P and A are in d1 frame.
   */
  PivotJointR(std::shared_ptr<siconos::algebra::SiconosVector> P,
              std::shared_ptr<siconos::algebra::SiconosVector> A, bool absoluteRef,
              std::shared_ptr<siconos::modeling::NewtonEulerDS> d1 = nullptr,
              std::shared_ptr<siconos::modeling::NewtonEulerDS> d2 = nullptr);

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

  virtual ~PivotJointR() noexcept = default;

  std::shared_ptr<siconos::algebra::SiconosVector> A() { return _A; }

  /**
     to compute the output y = h(q) of the Relation

     \param[in] q generalized coordinates vector of the dynamical systems (at most 2) involved
    in the relation \param[in,out] y the resulting vector
 */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;
  /**
       Get the number of constraints defined in the joint

       \return the number of constraints
     */
  virtual unsigned int numberOfConstraints() override { return 5; }

  /**
     Return the number of degrees of freedom of this joint.

     \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() override { return 1; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(unsigned int axis) override {
    if (axis == 0)
      return DofType::ANGULAR;
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
#endif  // PivotJointRELATION_H
