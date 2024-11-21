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
/** \file KneeJointR.hpp
 */
#ifndef KneeJointRELATION_H
#define KneeJointRELATION_H

#include <NewtonEulerJointR.hpp>

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::joints {

/**
   This class implements a knee joint between one or two Newton/Euler Dynamical system
*/
class KneeJointR : public NewtonEulerJointR {
 protected:
  ACCEPT_SERIALIZATION(KneeJointR);

  /** Coordinate of the knee point in the body frame of the first dynamical system _d1
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _P0{nullptr};

  /**Absolute coodinates of the vector  G1P0 when d1 is located in q=(0,0,0,1,0,0,0)
   * i.e. P0 in the body frame of d1.
   * These values are computed when the constructor is called.
   */
  double _G1P0x{0.};
  double _G1P0y{0.};
  double _G1P0z{0.};

  /** Absolute coodinates of the vector G2P0 when d2 is located in q=(0,0,0,1,0,0,0)
   *  i.e. P0 in the body frame of d2.
   * These values are computed when the constructor is called.
   */
  double _G2P0x{0.};
  double _G2P0y{0.};
  double _G2P0z{0.};

 public:
  /** Empty constructor. The relation may be initialized later by
   *  setPoint, setAbsolute, and setBasePositions. */
  KneeJointR();

  /** Constructor based on one or two dynamical systems and a point.
   *
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, or NULL
   *  for absolute frame.
   *  \param P SiconosVector of size 3 that defines the point around
   *  which rotation is allowed.
   *  \param absoluteRef if true, P is in the absolute frame,
   *  otherwise P is in d1 frame.
   */
  KneeJointR(std::shared_ptr<siconos::algebra::SiconosVector> P, bool absoluteRef,
             std::shared_ptr<siconos::modeling::NewtonEulerDS> d1 = nullptr,
             std::shared_ptr<siconos::modeling::NewtonEulerDS> d2 = nullptr);

  /** destructor
   */
  virtual ~KneeJointR() noexcept = default;

  virtual void initialize(siconos::modeling::Interaction& inter) override;

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

  /** Perform some checks on the initial conditions. */
  void checkInitPos(std::shared_ptr<siconos::algebra::SiconosVector> q1,
                    std::shared_ptr<siconos::algebra::SiconosVector> q2);

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() override { return 3; }

  /**
     Get the number of degrees of freedom defined in the joint

     \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() override { return 3; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(unsigned int axis) override {
    if (axis < 3)
      return DofType::ANGULAR;
    else
      return DofType::INVALID;
  };

  virtual void computeJacobianhOver_q(double time, siconos::modeling::Interaction& inter,
                            std::shared_ptr<siconos::algebra::BlockVector> q0) override;

  /**
     to compute the output y = h(t,q,z) of the Relation

     \param time current time value
     \param q coordinates of the dynamical systems involved in the relation
     \param y the resulting vector
  */
  virtual void computeh(double time, const siconos::algebra::BlockVector& q0,
                        siconos::algebra::SiconosVector& y) override;

  virtual void computeDotJachq(double time, const siconos::algebra::BlockVector& workQ,
                               siconos::algebra::BlockVector& workZ,
                               const siconos::algebra::BlockVector& workQdot) override;

  virtual void computeDotJachq(
      double time, std::shared_ptr<siconos::algebra::SiconosVector> qdot1,
      std::shared_ptr<siconos::algebra::SiconosVector> qdot2 = nullptr);

  std::shared_ptr<siconos::algebra::SiconosVector> P() { return _P0; }

 protected:
  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12,
                     double q13, double X2, double Y2, double Z2, double q20, double q21,
                     double q22, double q23);
  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12,
                   double q13);

  /* \warning, the following function should also depend on q */
  virtual void DotJd1d2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11,
                        double qdot12, double qdot13, double Xdot2, double Ydot2, double Zdot2,
                        double qdot20, double qdot21, double qdot22, double qdot23);
  virtual void DotJd1(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11,
                      double qdot12, double qdot13);

  // public:
  double Hx(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  double Hy(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  double Hz(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
};
}  // namespace siconos::joints
#endif  // KneeJointRELATION_H
