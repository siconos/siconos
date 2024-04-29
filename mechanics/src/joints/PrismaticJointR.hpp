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
/*! \file PrismaticJointR.hpp

*/
#ifndef PrismaticJointRELATION_H
#define PrismaticJointRELATION_H

#include "NewtonEulerJointR.hpp"

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::joints {
/**
   This class implements a prismatic joint between one or two Newton/Euler Dynamical system

   From a given axis, we construct two unit othorgonal vectors to the axis V1 and V2 such that
   (axis,V1,V2) is an orthogonal frame

*/
class PrismaticJointR : public NewtonEulerJointR {
 protected:
  ACCEPT_SERIALIZATION(PrismaticJointR);

  /** Axis of the prismatic point in the q1 frame of reference
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _axis0{nullptr};

  /** _V1 is an unit vector that is orthogonal to the prismatic axis _axis0.
   * It forms with _V2 and _axis0 a base such that (_axis0,_V1,_v2) is an orthogonal
   * frame
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _V1{nullptr};

  /** _V2 is an unit vector that is orthogonal to the prismatic axis _axis0.
   * It forms with _V2 and _axis0 a base such that (_axis0,_V1,_v2) is an orthogonal
   * frame
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _V2{nullptr};

  /** Convenient storage of the components of _V1 and _V2
   */
  double _V1x{0.};
  double _V1y{0.};
  double _V1z{0.};
  double _V2x{0.};
  double _V2y{0.};
  double _V2z{0.};

  double _G10G20d1x{0.};
  double _G10G20d1y{0.};
  double _G10G20d1z{0.};

  double _cq2q101{0.};
  double _cq2q102{0.};
  double _cq2q103{0.};
  double _cq2q104{0.};

  /** Return the normal of the linear DoF axis.
   *
   *  \param ans
   *  \param q0
   *  \param axis must be 0
   *  \param absoluteRef
   */
  virtual void _normalDoF(siconos::algebra::SiconosVector& ans,
                          const siconos::algebra::BlockVector& q0, int axis,
                          bool absoluteRef = true) override;

 public:
  /** Empty constructor. The relation may be initialized later by
   *  setPoint, setAbsolute, and setBasePositions. */
  PrismaticJointR();

  /** Constructor based on one or two dynamical systems and an axis.
   *
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, or NULL
   *  for absolute frame.
   *  \param axis siconos::algebra::SiconosVector of size 3 that defines the prismatic axis.
   *  \param absoluteRef if true, A is in the absolute frame,
   *  otherwise A is in d1 frame.
   */
  PrismaticJointR(std::shared_ptr<siconos::algebra::SiconosVector> axis, bool absoluteRef,
                  std::shared_ptr<siconos::modeling::NewtonEulerDS> d1 = nullptr,
                  std::shared_ptr<siconos::modeling::NewtonEulerDS> d2 = nullptr);

  /** destructor */
  virtual ~PrismaticJointR() noexcept = default;

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

  void displayInitialPosition();

  void computeV1V2FromAxis();

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, const siconos::algebra::BlockVector& q0,
                           siconos::algebra::SiconosVector& y, unsigned int axis) override;

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, siconos::modeling::Interaction& inter,
                               std::shared_ptr<siconos::algebra::BlockVector> q0,
                               siconos::algebra::SiconosMatrix& jachq,
                               unsigned int axis) override;

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

  virtual void computeDotJachq(double time, const siconos::algebra::BlockVector& workQ,
                               siconos::algebra::BlockVector& workZ,
                               const siconos::algebra::BlockVector& workQdot) override;

  /* The options were    : operatorarrow */
  double H1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H3(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  /* The options were    : operatorarrow */
  double H5(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H4(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
             double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);

  void DotJd1d2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11,
                double qdot12, double qdot13, double Xdot2, double Ydot2, double Zdot2,
                double qdot20, double qdot21, double qdot22, double qdot23);

  void DotJd2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11,
              double qdot12, double qdot13, double X2, double Y2, double Z2, double qdot20,
              double qdot21, double qdot22, double qdot23);

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
  */
  virtual unsigned int numberOfConstraints() override { return 5; }

  /**
     Return the number of degrees of freedom of this joint.

     \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() override { return 2; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(unsigned int axis) override {
    if (axis == 0)
      return DofType::LINEAR;
    else
      return DofType::INVALID;
  };
};
}  // namespace siconos::joints
#endif  // PrismaticJointRELATION_H
