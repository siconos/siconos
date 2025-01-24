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

/*! \file CylindricalJointR.hpp
 */

#ifndef CylindricalJointRELATION_H
#define CylindricalJointRELATION_H
#include <NewtonEulerJointR.hpp>

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::joints {

/**
    This class implements a cylindrical joint between one or
    two Newton/Euler Dynamical system.  It is similar to a
    PrismaticJointR but allows for rotation around the axis.

    From a given axis, we construct two unit othorgonal vectors to the
    axis V1 and V2 such that (axis,V1,V2) is an orthogonal frame
*/
class CylindricalJointR : public NewtonEulerJointR {
 protected:
  ACCEPT_SERIALIZATION(CylindricalJointR);

  /** Axis of the cylindrical point in the q1 frame of reference
   */
  std::shared_ptr<siconos::algebra::SiconosVector3> _axis0{nullptr};

  /** _V1 is an unit vector that is orthogonal to the cylindrical axis
   * _axis0.  It forms with _V2 and _axis0 a base such that
   * (_axis0,_V1,_v2) is an orthogonal frame
   */
  std::shared_ptr<siconos::algebra::SiconosVector3> _V1{nullptr};

  /** _V2 is an unit vector that is orthogonal to the cylindrical axis
   * _axis0.  It forms with _V2 and _axis0 a base such that
   * (_axis0,_V1,_v2) is an orthogonal frame
   */
  std::shared_ptr<siconos::algebra::SiconosVector3> _V2{nullptr};

  double _cq2q101{0.};
  double _cq2q102{0.};
  double _cq2q103{0.};
  double _cq2q104{0.};

  /** P is the point defining the location of the line created by
   * _axis0.  It is stored in the q1 frame, i.e. the vector from
   * initial G1 to P, called _G1P0. */
  std::shared_ptr<siconos::algebra::SiconosVector3> _G1P0{nullptr};

  /** _G2P0 is the vector from initial G1 to P */
  std::shared_ptr<siconos::algebra::SiconosVector3> _G2P0{nullptr};

  /** Cumulative number of twists around the joint relative to initial
   * angular difference. */
  int _twistCount{0};         // TODO: Should be in a graph work vector?
  double _previousAngle{0.};  // Needed to track _twistCount, TODO: work vector?
  double _initialAngle{0.};

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                                       const siconos::algebra::BlockVector& q0) override;

 public:
  /** Constructor based on one or two dynamical systems, a point and an axis.
   *
   *  \param P a vector to define the point around which rotation is allowed.
   *  \param A a vector to define the cylindrical axis.
   *  \param absoluteRef if true, P and A are in the absolute frame, otherwise P and A are in
   *  d1 frame.
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint (nullptr for absolute frame)
   */
  CylindricalJointR(const Eigen::Ref<siconos::algebra::SiconosVector3>& P,
                    const Eigen::Ref<siconos::algebra::SiconosVector3>& A, bool absoluteRef,
                    std::shared_ptr<siconos::modeling::NewtonEulerDS> d1 = nullptr,
                    std::shared_ptr<siconos::modeling::NewtonEulerDS> d2 = nullptr);

  /** destructor */
  virtual ~CylindricalJointR() noexcept = default;

  /** Initialize the joint constants based on the provided base positions.
   *
   *  \param[in] q1 a vector of size 7 indicating translation and orientation in inertial
   * coordinates.
   *  \param[in] q2 an optional vector of size 7 indicating translation and orientation; if
   * null, the inertial frame will be considered as the second base.
   */
  virtual void setBasePositions(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2 =
          std::nullopt) override;

  /** \return the axis of rotation.
   *  Retrieve a normal in the direction of a 0-indexed free
   *  axis. Useful for calculating velocities in the axis, or for
   *  calculating axis-aligned forces applied to connected bodies.  If
   *  axis is of angular type (see typeOfDoF), then the returned normal
   *  is the axis of rotation.
   *
   *  \param[in] q0 The state q of one or more NewtonEulerDS
   *  \param[in] axis
   *  \param[in] absoluteRef If true, ans is in the inertial frame,
   *  otherwise the q1 frame is assumed.
   */
  virtual siconos::algebra::SiconosVector3 normalDoF(const siconos::algebra::BlockVector& q0,
                                                     int axis,
                                                     bool absoluteRef = true) override;

  void computeV1V2FromAxis();

  int twistCount() { return _twistCount; }

  /**
     to compute the output y = h(q) of the Relation

     \param[in] q generalized coordinates vector of the dynamical systems (at most 2) involved
    in the relation \param[in,out] y the resulting vector
 */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
             double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);

  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() const override { return 4; }

  /**
     Return the number of degrees of freedom of this joint.

     \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() const override { return 2; }

  /**
     Return the type of a degree of freedom of this joint.

     \return the type of the degree of freedom (DoF)
  */
  virtual DofType typeOfDoF(unsigned int axis) const override {
    if (axis == 0)
      return DofType::LINEAR;
    else if (axis == 1)
      return DofType::ANGULAR;
    else
      return DofType::INVALID;
  }
  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y, unsigned int axis = 0) override;

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(siconos::modeling::Interaction& inter,
                               const siconos::algebra::BlockVector& q0,
                               Eigen::Ref<siconos::algebra::SiconosMatrix> jachq,
                               unsigned int axis = 0) override;
};
}  // namespace siconos::joints
#endif  // CylindricalJointRELATION_H
