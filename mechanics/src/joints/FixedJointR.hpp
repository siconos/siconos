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

  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, siconos::modeling::Interaction& inter,
                            const siconos::algebra::BlockVector& q0) override;

 public:
  /** default constructor */
  FixedJointR() = default;

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
   *  \param[in] q1 a vector of size 7 indicating translation and orientation in inertial
   * coordinates.
   *  \param[in] q2 an optional vector of size 7 indicating translation and orientation; if
   * null, the inertial frame will be considered as the second base.
   */
  virtual void setBasePositions(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2 =
          std::nullopt) override;
  /**
     Get the number of constraints defined in the joint

     \return the number of constraints
   */
  virtual siconos::algebra::Index numberOfConstraints() const override { return 6; }

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

  virtual siconos::algebra::Index numberOfDoF() const override { return 0; }

  virtual DofType typeOfDoF(unsigned int axis) const override { return DofType::INVALID; }

 protected:
  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12,
                     double q13, double X2, double Y2, double Z2, double q20, double q21,
                     double q22, double q23);

  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12,
                   double q13);

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::joints
#endif  // FixedJointRELATION_H
